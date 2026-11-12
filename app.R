# app.R  -  App Shiny: arricchimento di siti di legame di TF nei promotori umani.
#
# Come e' fatta un'app Shiny: tre pezzi.
#   1. cio' che si fa UNA volta all'avvio (qui: caricare i dati);
#   2. `ui`: come appare la pagina;
#   3. `server`: cosa succede quando l'utente fa qualcosa.
# Le funzioni di R/ (analyze_list, MAX_LIST, ...) non vanno caricate con source():
# Shiny carica da solo i file della cartella R/ prima di leggere questo file.

library(shiny)
library(bslib)
library(DT)                               # tabelle interattive (DTOutput / renderDT)
library(ggplot2)                          # il grafico di confronto (compare_plot, in R/compare.R)

# --- 1. Dati ------------------------------------------------------------------
# Fuori da server(): vengono letti una volta sola quando parte l'app e condivisi
# da tutti gli utenti. Dentro server() verrebbero riletti per ogni visitatore.
genes <- readRDS("data/genes.rds")        # una riga per gene
hits <- readRDS("data/motif_hits.rds")    # una voce per motivo, con i geni che hanno il sito

# --- 2. UI --------------------------------------------------------------------
# page_navbar() fa una pagina con una barra in alto e una scheda (nav_panel) per voce.
ui <- page_navbar(
  title = "tfbalance",
  id = "tabs",
  fillable = FALSE,                       # la pagina scorre normalmente
  nav_panel(
    "Analysis",
    layout_sidebar(                       # colonna a sinistra + contenuto a destra
      sidebar = sidebar(
        textAreaInput("genes", "Gene list", rows = 10,
                      placeholder = "One gene per line, or separated by spaces or commas"),
        helpText("Symbols, Entrez IDs or Ensembl IDs. At most ", MAX_LIST,
                 " recognised genes."),
        actionButton("go", "Analyse", class = "btn-primary")
      ),
      card(
        card_header("Results"),
        uiOutput("summary")               # qui il server scrive il riquadro dei conteggi
      ),
      card(
        card_header("Comparison"),
        plotOutput("plot"),               # grafico: baseline contro corretto, un punto per motivo
        uiOutput("plot_note")             # riga sotto il grafico: quanti motivi sono disegnati
      ),
      card(
        card_header("Motifs"),
        uiOutput("caption"),              # didascalia sopra la tabella
        DTOutput("tbl")                   # la tabella (DT = DataTables)
      )
    )
  ),
  nav_panel("Unrecognised genes", uiOutput("notfound")),
  nav_panel(
    "How it works",                       # testo di spiegazione (inglese, ASCII): titoli h4, paragrafi p, elenchi ul/li
    h4("What this app does"),
    p(paste0(
      'Paste a list of human genes (symbols, Entrez IDs or Ensembl IDs). For each transcription ',
      'factor (TF) motif, the app asks whether the promoters of your genes contain predicted ',
      'binding sites more often than expected by chance. A promoter is the region from 1000 bp ',
      'upstream to 100 bp downstream of the transcription start site (GRCh38). Sites come from ',
      'scanning 1019 JASPAR 2026 vertebrate motifs at 90% of the maximum score; 852 motifs are ',
      'tested (those found in 0.5% to 90% of promoters).')),
    h4("Two tests, side by side"),
    tags$ul(
      tags$li(paste0(
        'Baseline: a one-sided hypergeometric test on "genes with at least one site". It does ',
        'not account for the fact that promoters differ in GC content and in how many sites they ',
        'contain overall, which can inflate some results.')),
      tags$li(paste0(
        'Corrected: a score test of the same question in a logistic model that also includes GC ',
        'content and the number of predicted sites per promoter (log scale, leaving out the 109 ',
        'motifs found in more than 90% of promoters).'))
    ),
    p(paste0('Both tests look only for enrichment, not depletion. P-values are adjusted with ',
             'Benjamini-Hochberg (FDR).')),
    h4("Reading the results"),
    tags$ul(
      tags$li(paste0(
        'Table: k = genes of your list with a site; K = genes in the background with a site; ',
        'Expected (E) = number of your genes expected to have a site by chance; p baseline and ',
        'FDR baseline; z, p corrected and FDR corrected. NA in the corrected test means it was ',
        'not computed (for example because E is below 10).')),
      tags$li(paste0(
        'Counts box: number of motifs with FDR < 0.05 for the baseline on all motifs, the ',
        'baseline on the same motifs that have a corrected test (the fair comparison), and the ',
        'corrected test.')),
      tags$li(paste0(
        'Plot: one point per motif that has a corrected test, baseline on the x axis and ',
        'corrected on the y axis (-log10 p). Points below the dashed diagonal are less ',
        'significant after correction. Red lines mark FDR 0.05 for each test that has ',
        'significant motifs. Triangles show p = 0 (underflow), drawn at the plot limit.'))
    ),
    h4("Limits"),
    p(paste0(
      'Lists with fewer than 12 genes have no corrected test, and the maximum is 5000 genes. ',
      'The corrected test also cannot be fitted for lists made only of extreme genes (for ',
      'example the highest GC content). Motifs share genes, so the tests are not independent. ',
      'The corrected test was checked mainly on random lists of 300 genes. Predicted sites are ',
      'not proof of binding: use the results to generate hypotheses.'))
  )
)

# --- 3. Server ----------------------------------------------------------------
server <- function(input, output, session) {

  # `res` e' l'UNICO calcolo dell'app. Un "reactive" e' un calcolo che Shiny rifa
  # da solo quando cambia qualcosa da cui dipende. bindEvent(input$go) dice: rifallo
  # SOLO quando si preme "Analizza" (input$go conta i click). Scrivere nel campo di
  # testo non lo rilancia, e prima del primo click non parte.
  res <- reactive({
    req(!is.null(input$genes))            # se il testo non esiste ancora, aspetta
    withProgress(message = "Analysing...", {   # messaggio mentre calcola
      tryCatch(
        # suppressWarnings: l'avviso "the null model cannot be fitted" non serve, c'e' gia'
        # la nota corrected_note nel risultato.
        suppressWarnings(analyze_list(input$genes, genes, hits)),
        # Se analyze_list da' un errore (messaggio italiano) lo salviamo invece di
        # far cadere l'app: res() sara' list(error = "...").
        error = function(e) list(error = conditionMessage(e))
      )
    })
  }) |> bindEvent(input$go)

  # Il riquadro: un errore in rosso, oppure i conteggi.
  output$summary <- renderUI({
    r <- res()
    if (!is.null(r$error))
      return(div(class = "alert alert-danger", r$error))

    # Conteggi dei motivi con FDR < 0.05 (R/compare.R). B e C sono NA se nessun motivo ha il
    # test corretto: si scrive "not available", non 0.
    cnt <- fdr_counts(r$table)
    show <- function(x) if (is.na(x)) "not available" else x
    tagList(
      tags$ul(
        tags$li(paste0("Recognised genes: ", r$n_list)),
        tags$li(paste0("Repeated tokens or synonyms ignored: ", r$n_duplicates)),
        tags$li(paste0("Unrecognised: ", length(r$not_found))),
        tags$li(paste0("Motifs tested: ", r$n_motifs, " (excluded: ", r$n_excluded_low,
                       " below the minimum, ", r$n_excluded_high, " above the maximum)")),
        tags$li(paste0("Corrected test available for ", r$n_corrected, " of ",
                       r$n_motifs, " motifs"))
      ),
      p("Motifs with FDR < 0.05"),
      tags$ul(
        tags$li(paste0("baseline, all ", cnt$n_all, " motifs: ", cnt$base_all)),
        tags$li(paste0("baseline, same ", cnt$n_sub, " motifs as the corrected test: ",
                       show(cnt$base_sub))),
        tags$li(paste0("corrected test: ", show(cnt$corrected)))
      ),
      if (!is.null(r$corrected_note))     # la nota c'e' solo se nessun motivo ha il test corretto
        div(class = "alert alert-warning", r$corrected_note)
    )
  })

  # Il grafico. compare_plot (R/compare.R) disegna un punto per ogni motivo con il test corretto.
  # - se analyze_list e' fallita: req() ferma l'output in silenzio (l'errore e' gia' nel riquadro);
  # - se nessun motivo ha il test corretto: validate(need(...)) mostra un messaggio al posto del grafico;
  # - alt e' il testo per chi non vede il grafico (screen reader); e' una funzione, cosi' Shiny la
  #   rilegge a ogni analisi e il numero di motivi e' quello della tabella.
  output$plot <- renderPlot({
    r <- res()
    req(is.null(r$error))
    validate(need(any(!is.na(r$table$z)), "No motif has a corrected test: plot not available"))
    compare_plot(r$table)
  }, alt = function() {
    sprintf(paste0("Scatter plot of %d motifs: -log10 p baseline (x) against -log10 p corrected (y), ",
                   "with the diagonal y = x. A line marks FDR 0.05 for each test that has ",
                   "significant motifs; triangles mark p = 0 (underflow)."),
            sum(!is.na(res()$table$z)))
  })

  # Riga sotto il grafico: N = motivi con z non NA (quelli disegnati), M = motivi testati.
  output$plot_note <- renderUI({
    r <- res()
    req(is.null(r$error))
    p(class = "text-muted small",
      sprintf("%d motifs drawn out of %d (the others have no corrected test)",
              sum(!is.na(r$table$z)), nrow(r$table)))
  })

  # Didascalia: spiega gli NA, ma solo con frasi VERE. "NA = atteso < soglia" vale per le
  # liste normali; gli NA possono pero' nascere anche da altre cause (modello non
  # stimabile, motivo presente in tutti i geni, ...). Quindi la frase precisa si usa solo se
  # controllando la tabella e' vera: nessuna nota e nessun NA con atteso >= soglia.
  # La soglia non e' scritta a mano: la leggo dal valore di default di min_expected nella
  # firma di analyze_list (che l'app non cambia), cosi' se quel default cambia, cambia
  # anche il testo.
  output$caption <- renderUI({
    r <- res()
    req(is.null(r$error))                 # con un errore non c'e' tabella: niente didascalia
    soglia <- formals(analyze_list)$min_expected
    na_con_atteso_alto <- any(is.na(r$table$z) & r$table$expected >= soglia)
    testo <- if (is.null(r$corrected_note) && !na_con_atteso_alto)
      paste0("NA in the corrected test = expected < ", format(soglia), " genes: test not reliable.")
    else
      paste0("NA in the corrected test = test not computed (the cause is in the note above; ",
             "with no note, the motif cannot be analysed, for example because it is present in all genes).")
    p(class = "text-muted small", testo)
  })

  # La tabella. Il risultato di analyze_list e' gia' ordinato: order = list() dice a DT di
  # NON riordinarlo. (DT 0.34.0 lo fa gia' da solo se manca; lo scrivo lo stesso perche'
  # l'ordine e' una scelta dell'app e non deve dipendere dal default di DT.)
  # formatSignif/formatRound cambiano solo come si VEDONO i numeri (3 cifre
  # significative per p e FDR, 2 decimali per atteso e z): l'ordinamento resta numerico.
  output$tbl <- renderDT({
    r <- res()
    req(is.null(r$error))
    datatable(
      r$table,
      rownames = FALSE,                   # niente colonna con i numeri di riga
      colnames = c("Motif", "TF", "k", "K", "Expected (E)", "p baseline", "FDR baseline",
                   "z", "p corrected", "FDR corrected"),
      options = list(order = list(), pageLength = 25)
    ) |>
      formatSignif(c("p_baseline", "fdr_baseline", "p_corrected", "fdr_corrected"), 3) |>
      formatRound(c("expected", "z"), 2)
  })

  # Scheda "Unrecognised genes": il numero e l'elenco, un nome per riga.
  # tags$pre mostra il testo cosi' com'e' (anche un nome con "<" non viene letto come HTML).
  output$notfound <- renderUI({
    r <- res()
    if (!is.null(r$error))
      return(p("No analysis available: see the error in the Analysis tab."))
    if (length(r$not_found) == 0)
      return(p("No unrecognised tokens."))
    tagList(
      p(paste0("Unrecognised tokens: ", length(r$not_found))),
      tags$pre(paste(r$not_found, collapse = "\n"))
    )
  })
}

shinyApp(ui, server)
