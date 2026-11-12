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
    "Analisi",
    layout_sidebar(                       # colonna a sinistra + contenuto a destra
      sidebar = sidebar(
        textAreaInput("genes", "Lista di geni", rows = 10,
                      placeholder = "Un gene per riga, o separati da spazi o virgole"),
        helpText("Simboli, ID Entrez o ID Ensembl. Massimo ", MAX_LIST,
                 " geni riconosciuti."),
        actionButton("go", "Analizza", class = "btn-primary")
      ),
      card(
        card_header("Risultato"),
        uiOutput("summary")               # qui il server scrive il riquadro dei conteggi
      ),
      card(
        card_header("Motivi"),
        uiOutput("caption"),              # didascalia sopra la tabella
        DTOutput("tbl")                   # la tabella (DT = DataTables)
      )
    )
  ),
  nav_panel("Geni non riconosciuti", uiOutput("notfound")),
  nav_panel("Come funziona", p("(testo da scrivere)"))     # segnaposto
)

# --- 3. Server ----------------------------------------------------------------
server <- function(input, output, session) {

  # `res` e' l'UNICO calcolo dell'app. Un "reactive" e' un calcolo che Shiny rifa
  # da solo quando cambia qualcosa da cui dipende. bindEvent(input$go) dice: rifallo
  # SOLO quando si preme "Analizza" (input$go conta i click). Scrivere nel campo di
  # testo non lo rilancia, e prima del primo click non parte.
  res <- reactive({
    req(!is.null(input$genes))            # se il testo non esiste ancora, aspetta
    withProgress(message = "Analisi in corso...", {   # messaggio mentre calcola
      tryCatch(
        # suppressWarnings: l'avviso "modello non stimabile" non serve, c'e' gia'
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

    tagList(
      tags$ul(
        tags$li(paste0("Geni riconosciuti: ", r$n_list)),
        tags$li(paste0("Token ripetuti o sinonimi ignorati: ", r$n_duplicates)),
        tags$li(paste0("Non riconosciuti: ", length(r$not_found))),
        tags$li(paste0("Motivi testati: ", r$n_motifs, " (esclusi: ", r$n_excluded_low,
                       " sotto il minimo, ", r$n_excluded_high, " sopra il massimo)")),
        tags$li(paste0("Test corretto disponibile per ", r$n_corrected, " di ",
                       r$n_motifs, " motivi"))
      ),
      if (!is.null(r$corrected_note))     # la nota c'e' solo se nessun motivo ha il test corretto
        div(class = "alert alert-warning", r$corrected_note)
    )
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
      paste0("NA nel test corretto = atteso < ", format(soglia), " geni: test non affidabile.")
    else
      paste0("NA nel test corretto = test non calcolato (la causa e' nella nota sopra; ",
             "senza nota, il motivo non e' analizzabile, per esempio e' presente in tutti i geni).")
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
      colnames = c("Motivo", "TF", "k", "K", "Attesi (E)", "p baseline", "FDR baseline",
                   "z", "p corretto", "FDR corretto"),
      options = list(order = list(), pageLength = 25)
    ) |>
      formatSignif(c("p_baseline", "fdr_baseline", "p_corrected", "fdr_corrected"), 3) |>
      formatRound(c("expected", "z"), 2)
  })

  # Scheda "Geni non riconosciuti": il numero e l'elenco, un nome per riga.
  # tags$pre mostra il testo cosi' com'e' (anche un nome con "<" non viene letto come HTML).
  output$notfound <- renderUI({
    r <- res()
    if (!is.null(r$error))
      return(p("Nessuna analisi disponibile: guarda l'errore nella scheda Analisi."))
    if (length(r$not_found) == 0)
      return(p("Nessun token non riconosciuto."))
    tagList(
      p(paste0("Token non riconosciuti: ", length(r$not_found))),
      tags$pre(paste(r$not_found, collapse = "\n"))
    )
  })
}

shinyApp(ui, server)
