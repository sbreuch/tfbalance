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

# --- 1. Dati ------------------------------------------------------------------
# Fuori da server(): vengono letti una volta sola quando parte l'app e condivisi
# da tutti gli utenti. Dentro server() verrebbero riletti per ogni visitatore.
genes <- readRDS("data/genes.rds")        # una riga per gene
hits <- readRDS("data/motif_hits.rds")    # una voce per motivo, con i geni che hanno il sito

# --- 2. UI --------------------------------------------------------------------
# page_navbar() fa una pagina con una barra in alto e una scheda (nav_panel) per
# voce. Per ora c'e' solo "Analisi"; le altre arrivano nel passo successivo.
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
      )
    )
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
}

shinyApp(ui, server)
