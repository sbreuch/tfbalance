# tests/testthat/test-app.R
# Test di app.R scritti PRIMA del codice (Fase B, passo B2a), senza browser:
# shiny::testServer() fa girare il server dell'app con input finti.
# Usano i dati veri (data/*.rds): se mancano, i test vengono saltati.
#
# Cosa si aspetta da app.R (interfaccia):
#   input$genes (testo), input$go (bottone)
#   res()            reactive agganciato a input$go: il risultato di analyze_list(),
#                    oppure list(error = <messaggio italiano>) se analyze_list fallisce
#   output$summary   riquadro dei conteggi (o dell'errore)

# Radice dell'app vista da tests/testthat: due cartelle sopra. test_path() funziona
# comunque lo si lanci (da tests/testthat o dalla radice del progetto).
# normalizePath() lo rende ASSOLUTO: mentre l'app gira Shiny cambia la cartella di
# lavoro (la sposta nella radice dell'app), e un percorso relativo si romperebbe.
app_dir <- normalizePath(testthat::test_path("..", ".."))
data_ok <- file.exists(file.path(app_dir, "data", "genes.rds")) &&
           file.exists(file.path(app_dir, "data", "motif_hits.rds"))

# n simboli veri e unici (scelti a caso, seme fisso): ognuno e' un gene diverso
real_symbols <- function(n, seed = 1) {
  g <- readRDS(file.path(app_dir, "data", "genes.rds"))
  s <- toupper(g$symbol)
  ok <- !is.na(s) & s != "" & !duplicated(s) & !duplicated(s, fromLast = TRUE)
  set.seed(seed)
  sample(g$symbol[ok], n)
}

# testo semplice del riquadro: dall'HTML tolgo i tag, rimetto i caratteri speciali
# (nell'HTML ">" e' scritto "&gt;", "'" come "&#39;") e riduco gli spazi
box_text <- function(x) {
  s <- gsub("<[^>]+>", " ", as.character(x$html))
  s <- gsub("&gt;", ">", s, fixed = TRUE)
  s <- gsub("&lt;", "<", s, fixed = TRUE)
  s <- gsub("&#39;", "'", s, fixed = TRUE)
  s <- gsub("&amp;", "&", s, fixed = TRUE)
  trimws(gsub("\\s+", " ", s))
}

# Il codice dentro testServer() vede solo gli oggetti del server (input, res, ...),
# non quelli di R/ ne' i dati. Si raggiungono dall'ambiente padre del server,
# cioe' quello di app.R, che a sua volta eredita quello dei file di R/.
app_get <- function(name, session) get(name, envir = parent.env(session$env))

COLS <- c("id", "tf", "k", "K", "expected", "p_baseline", "fdr_baseline",
          "z", "p_corrected", "fdr_corrected")

# --- (0) il percorso e' giusto ----------------------------------------------
# Serve a distinguere, se gli altri test sono rossi, "percorso sbagliato" da
# "app.R non c'e' ancora" (questo test e' verde da subito).
test_that("app_dir e' la radice del progetto (ci sono R/analyze.R e R/enrich.R)", {
  expect_true(file.exists(file.path(app_dir, "R", "analyze.R")))
  expect_true(file.exists(file.path(app_dir, "R", "enrich.R")))
})

# --- (a) R/ viene caricato da Shiny, senza source() -----------------------------
test_that("(a) le funzioni di R/ sono visibili al server e vengono dal caricamento automatico", {
  skip_if_not(data_ok)
  shiny::testServer(app_dir, {
    f <- app_get("analyze_list", session)
    expect_true(is.function(f))
    expect_equal(app_get("MAX_LIST", session), 5000)
    # il caricamento automatico le valuta in un ambiente figlio di globalenv, non in
    # globalenv: cosi' non si confonde con un source() fatto da un altro file di test
    expect_false(identical(environment(f), globalenv()))
    expect_identical(parent.env(environment(f)), globalenv())
    # i dati sono stati letti fuori da server(), una volta
    expect_gt(nrow(app_get("genes", session)), 1000)
  })
})

# --- (b) lista valida ------------------------------------------------------------
test_that("(b) lista valida: res() e' il risultato di analyze_list", {
  skip_if_not(data_ok)
  shiny::testServer(app_dir, {
    session$setInputs(genes = paste(real_symbols(60), collapse = " "), go = 1)
    r <- res()
    expect_null(r$error)
    expect_equal(r$n_list, 60)
    expect_identical(names(r$table), COLS)
    expect_equal(nrow(r$table), r$n_motifs)
    expect_false(is.unsorted(r$table$p_corrected[!is.na(r$table$p_corrected)]))
  })
})

test_that("(b) testo del riquadro: conteggi e 'test corretto disponibile per X di Y motivi'", {
  skip_if_not(data_ok)
  shiny::testServer(app_dir, {
    # 60 simboli veri + 2 inventati, uno scritto due volte
    txt <- paste(c(real_symbols(60), "FOOBAR", "FOOBAR"), collapse = ", ")
    session$setInputs(genes = txt, go = 1)
    r <- res()
    s <- box_text(output$summary)
    expect_match(s, "Geni riconosciuti: 60", fixed = TRUE)
    expect_match(s, "Token ripetuti o sinonimi ignorati: 1", fixed = TRUE)
    expect_match(s, "Non riconosciuti: 1", fixed = TRUE)
    expect_match(s, sprintf("Motivi testati: %d (esclusi: %d sotto il minimo, %d sopra il massimo)",
                            r$n_motifs, r$n_excluded_low, r$n_excluded_high), fixed = TRUE)
    # X e Y sono n_corrected e n_motifs di analyze_list, non numeri scritti a mano
    expect_match(tolower(s), sprintf("test corretto disponibile per %d di %d motivi",
                                     r$n_corrected, r$n_motifs), fixed = TRUE)
    expect_gt(r$n_corrected, 0)                    # con 60 geni qualche motivo ha E >= 10
    expect_null(r$corrected_note)
    expect_false(grepl("alert-danger", as.character(output$summary$html)))
  })
})

test_that("(b) lista piccola: 'disponibile per 0 di Y motivi' e compare corrected_note", {
  skip_if_not(data_ok)
  shiny::testServer(app_dir, {
    session$setInputs(genes = paste(real_symbols(10), collapse = " "), go = 1)
    r <- res()
    expect_equal(r$n_corrected, 0)                 # 10 geni: E <= 9 per ogni motivo
    s <- box_text(output$summary)
    expect_match(tolower(s), sprintf("test corretto disponibile per 0 di %d motivi", r$n_motifs),
                 fixed = TRUE)
    expect_type(r$corrected_note, "character")
    expect_match(s, r$corrected_note, fixed = TRUE)   # la nota e' mostrata nel riquadro
  })
})

# --- (c) errori di analyze_list ---------------------------------------------------
test_that("(c) nessun gene riconosciuto: res()$error e riquadro rosso con il messaggio", {
  skip_if_not(data_ok)
  shiny::testServer(app_dir, {
    session$setInputs(genes = "FOO BAR", go = 1)
    r <- res()
    expect_match(r$error, "nessun gene riconosciuto")
    expect_null(r$table)
    expect_match(as.character(output$summary$html), "alert-danger", fixed = TRUE)
    s <- box_text(output$summary)
    expect_match(s, "nessun gene riconosciuto", fixed = TRUE)
    expect_false(grepl("Geni riconosciuti", s))    # niente conteggi accanto all'errore
  })
})

test_that("(c) testo vuoto: arriva ad analyze_list e l'utente vede l'errore (non resta muto)", {
  skip_if_not(data_ok)
  shiny::testServer(app_dir, {
    session$setInputs(genes = "", go = 1)
    expect_match(res()$error, "nessun gene riconosciuto")
    expect_match(box_text(output$summary), "nessun gene riconosciuto", fixed = TRUE)
  })
})

test_that("(c) un errore sui limiti arriva com'e': piu' di MAX_TOKENS token", {
  skip_if_not(data_ok)
  shiny::testServer(app_dir, {
    max_tokens <- app_get("MAX_TOKENS", session)
    session$setInputs(genes = paste(rep("G1", max_tokens + 1), collapse = " "), go = 1)
    expect_match(res()$error, "50000|50\\.000")
    expect_match(box_text(output$summary), "50000|50\\.000")
  })
})

# --- (d) l'analisi parte solo dal bottone -------------------------------------------
test_that("(d) senza click non si calcola; scrivere nel campo non rilancia l'analisi", {
  skip_if_not(data_ok)
  # Il bottone vero conta i click (1, 2, 3...), quindi simulo due click con go = 1 e
  # go = 2. NB: non serve per far funzionare il test: in Shiny 1.14.0 il server non
  # scarta un input con lo stesso valore di prima (ReactiveValues con dedupe = FALSE,
  # sia in ShinySession sia in MockShinySession), e con go = 1 due volte il test
  # darebbe lo stesso verde. Il contatore crescente rispecchia il browser.
  # (Negli altri test c'e' un solo click per sessione, quindi go = 1 basta.)
  shiny::testServer(app_dir, {
    # niente click: res() e' "in attesa" (errore silenzioso di Shiny), non un risultato
    expect_error(res(), class = "shiny.silent.error")
    session$setInputs(genes = paste(real_symbols(60), collapse = " "))
    expect_error(res(), class = "shiny.silent.error")   # testo scritto, ma nessun click

    session$setInputs(go = 1)
    expect_equal(res()$n_list, 60)

    # cambio il testo senza ricliccare: il risultato resta quello di prima
    session$setInputs(genes = paste(real_symbols(30, seed = 2), collapse = " "))
    expect_equal(res()$n_list, 60)

    # nuovo click: ora si ricalcola sul testo nuovo
    session$setInputs(go = 2)
    expect_equal(res()$n_list, 30)
  })
})
