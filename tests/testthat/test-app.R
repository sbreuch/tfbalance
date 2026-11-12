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

# nell'HTML ">" e' scritto "&gt;", "'" come "&#39;": li rimetto come testo normale
unescape_html <- function(s) {
  s <- gsub("&gt;", ">", s, fixed = TRUE)
  s <- gsub("&lt;", "<", s, fixed = TRUE)
  s <- gsub("&#39;", "'", s, fixed = TRUE)
  s <- gsub("&quot;", "\"", s, fixed = TRUE)
  gsub("&amp;", "&", s, fixed = TRUE)          # per ultimo
}

# testo semplice di un output HTML: tolgo i tag, rimetto i caratteri speciali, riduco gli spazi
box_text <- function(x) {
  trimws(gsub("\\s+", " ", unescape_html(gsub("<[^>]+>", " ", as.character(x$html)))))
}

# Il codice dentro testServer() vede solo gli oggetti del server (input, res, ...),
# non quelli di R/ ne' i dati. Si raggiungono dall'ambiente padre del server,
# cioe' quello di app.R, che a sua volta eredita quello dei file di R/.
app_get <- function(name, session) get(name, envir = parent.env(session$env))

# i n geni con GC piu' alto (simboli veri e unici): la lista e' separabile dalle covariate
# del modello, quindi il test corretto non e' stimabile
top_gc_symbols <- function(n) {
  g <- readRDS(file.path(app_dir, "data", "genes.rds"))
  s <- toupper(g$symbol)
  cand <- which(!is.na(s) & s != "" & !duplicated(s) & !duplicated(s, fromLast = TRUE))
  g$symbol[cand[order(g$gc[cand], decreasing = TRUE)][seq_len(n)]]
}

# output$tbl dentro testServer e' una stringa JSON: la leggo come lista
tbl_json <- function(x) jsonlite::fromJSON(unclass(x), simplifyVector = FALSE)

# intestazioni italiane della tabella, nell'ordine delle colonne di res()$table
TABLE_HEADERS <- c("Motivo", "TF", "k", "K", "Attesi (E)", "p baseline", "FDR baseline",
                   "z", "p corretto", "FDR corretto")

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

# === B2b ======================================================================
# Interfaccia aggiunta: output$tbl (tabella DT), output$caption (didascalia sopra la
# tabella), output$notfound (tab "Geni non riconosciuti"), oggetto `ui` in app.R.

# --- (e) lista separabile ----------------------------------------------------------
test_that("(e) lista separabile: nessun avviso arriva al test e la nota 'non e' stimabile' e' mostrata", {
  skip_if_not(data_ok)
  shiny::testServer(app_dir, {
    # enrich_corrected avverte "modello non stimabile" (verificato sui dati veri); l'app
    # lo silenzia con suppressWarnings(), la nota nel riquadro basta
    expect_no_warning({
      session$setInputs(genes = paste(top_gc_symbols(60), collapse = " "), go = 1)
      r <- res()
      s <- box_text(output$summary)
    })
    expect_equal(r$n_list, 60)
    expect_equal(r$n_corrected, 0)
    expect_match(r$corrected_note, "non e' stimabile", fixed = TRUE)
    expect_match(s, "non e' stimabile", fixed = TRUE)
    expect_match(s, "Test corretto disponibile per 0 di", fixed = TRUE)
    expect_false(grepl("nessun motivo", s))        # non e' un problema di E bassa
  })
})

# --- (f) tabella DT ----------------------------------------------------------------
test_that("(f) output$tbl: intestazioni italiane, formatSignif(3), formatRound(2), interi, order vuoto, rownames = FALSE", {
  skip_if_not(data_ok)
  shiny::testServer(app_dir, {
    session$setInputs(genes = paste(real_symbols(60), collapse = " "), go = 1)
    r <- res()
    # ATTENZIONE: questo test legge il JSON interno che DT produce (columnDefs, le funzioni
    # JavaScript DTWidget.formatSignif/formatRound). Se un aggiornamento di DT lo rompe
    # senza che l'app sia sbagliata, si adatta il test (non l'app).
    j <- tbl_json(output$tbl)
    opts <- j$x$options

    # rownames = FALSE: 10 intestazioni, nessuna colonna in piu' per i nomi di riga
    th <- regmatches(j$x$container, gregexpr("<th[^>]*>.*?</th>", j$x$container))[[1]]
    expect_identical(trimws(gsub("<[^>]+>", "", th)), TABLE_HEADERS)

    # formattazione: DT scrive una funzione JavaScript per colonna (targets = numero di
    # colonna da 0). Guardo solo le voci con `render` (le altre sono classi e nomi).
    cd <- Filter(function(x) !is.null(x$render), opts$columnDefs)
    target <- vapply(cd, function(x) as.numeric(x$targets), numeric(1))
    js <- vapply(cd, function(x) x$render, "")
    col0 <- function(nm) match(nm, names(r$table)) - 1      # posizione 0-based (senza rownames)
    sig <- c("p_baseline", "fdr_baseline", "p_corrected", "fdr_corrected")
    rnd <- c("expected", "z")
    expect_setequal(target[grepl("formatSignif(data, 3,", js, fixed = TRUE)], col0(sig))
    expect_setequal(target[grepl("formatRound(data, 2,", js, fixed = TRUE)], col0(rnd))
    expect_length(cd, length(sig) + length(rnd))             # nient'altro e' formattato
    # k e K: interi, senza formattazione (nessuna voce per le loro colonne)
    expect_type(r$table$k, "integer")
    expect_type(r$table$K, "integer")
    expect_false(any(target %in% col0(c("k", "K"))))

    # order vuoto: DT non riordina le righe, resta l'ordine di analyze_list. Controllo il
    # risultato (JSON con order vuoto), non come e' scritto in app.R: DT 0.34.0 mette
    # order = list() da solo se manca, quindi questo test non distingue le due cose.
    # Serve che l'opzione ESISTA (se mancasse, length() darebbe 0 lo stesso)
    expect_true("order" %in% names(opts))
    expect_length(opts$order, 0)
  })
})

test_that("(f) con un errore ne' tabella ne' didascalia (req ferma l'output)", {
  skip_if_not(data_ok)
  shiny::testServer(app_dir, {
    session$setInputs(genes = "FOO BAR", go = 1)
    expect_error(output$tbl, class = "shiny.silent.error")
    expect_error(output$caption, class = "shiny.silent.error")
  })
})

# --- (g) didascalia ----------------------------------------------------------------
test_that("(g) didascalia: NA nel test corretto = atteso < N geni, con N letto dalle formals", {
  skip_if_not(data_ok)
  shiny::testServer(app_dir, {
    # N e' il default di min_expected nella firma di analyze_list: se cambia li', cambia
    # anche la didascalia, e questo test continua a valere
    n <- formals(app_get("analyze_list", session))$min_expected
    expect_true(is.numeric(n))
    session$setInputs(genes = paste(real_symbols(60), collapse = " "), go = 1)
    s <- box_text(output$caption)
    expect_match(s, "NA nel test corretto", fixed = TRUE)
    expect_match(s, sprintf("atteso < %s geni", format(n)), fixed = TRUE)
    expect_match(s, "test non affidabile", fixed = TRUE)
    # e la frase e' vera: ogni NA della tabella ha E sotto la soglia
    r <- res()
    expect_true(any(is.na(r$table$z)))
    expect_true(all(r$table$expected[is.na(r$table$z)] < n))
  })
})

# La didascalia deve essere VERA in ogni caso. "NA = atteso < N" non lo e' sempre: gli NA nascono
# anche da altre cause (modello non stimabile, K = N, guardia SCORE_VAR_MIN: vedi R/enrich.R).
# Regola: se c'e' una corrected_note, o se qualche NA ha E >= soglia, la didascalia NON deve
# affermare che NA = atteso < N (la frase vale solo per la lista "normale").
test_that("(g2) con corrected_note la didascalia non afferma 'atteso < N'", {
  skip_if_not(data_ok)
  shiny::testServer(app_dir, {
    # due liste che producono una nota: 10 geni (nessun motivo con E >= soglia) e 60 geni con
    # GC piu' alto (modello non stimabile: tutti NA anche con E >= soglia)
    session$setInputs(genes = paste(real_symbols(10), collapse = " "), go = 1)
    expect_type(res()$corrected_note, "character")
    s <- box_text(output$caption)
    expect_match(s, "NA nel test corretto", fixed = TRUE)
    expect_false(grepl("atteso <", s, fixed = TRUE))

    session$setInputs(genes = paste(top_gc_symbols(60), collapse = " "), go = 2)
    expect_match(res()$corrected_note, "non e' stimabile", fixed = TRUE)
    expect_true(any(res()$table$expected >= formals(app_get("analyze_list", session))$min_expected))
    s <- box_text(output$caption)
    expect_match(s, "NA nel test corretto", fixed = TRUE)
    expect_false(grepl("atteso <", s, fixed = TRUE))     # qui NA anche con E alto: sarebbe falso
  })
})

test_that("(g3) un NA con E >= soglia senza nota: la didascalia non dice '= atteso < N'", {
  skip_if_not(data_ok)
  shiny::testServer(app_dir, {
    # Caso reale ma non raggiungibile dall'app con il filtro di default (3 motivi hanno K = N,
    # NA con E = n; verificato in data-raw/validation/step8_na_with_E.R): lo simulo sostituendo
    # analyze_list con una funzione che chiama la vera e poi mette NA al primo motivo, con E = 50.
    orig <- app_get("analyze_list", session)
    env <- environment(orig)                              # dove R/ ha definito analyze_list
    fake <- function(text, genes, hits, frac_range = c(0.005, 0.9), min_expected = 10) {
      r <- orig(text, genes, hits, frac_range, min_expected)
      r$table$expected[1] <- 50
      r$table[1, c("z", "p_corrected", "fdr_corrected")] <- NA_real_
      r
    }
    assign("analyze_list", fake, envir = env)
    on.exit(assign("analyze_list", orig, envir = env), add = TRUE)

    session$setInputs(genes = paste(real_symbols(60), collapse = " "), go = 1)
    r <- res()
    expect_null(r$corrected_note)                         # nessuna nota...
    expect_true(any(is.na(r$table$z) & r$table$expected >= 10))   # ...ma un NA con E alto
    s <- box_text(output$caption)
    expect_match(s, "NA nel test corretto", fixed = TRUE)
    expect_false(grepl("= atteso <", s, fixed = TRUE))
  })
})

test_that("(g4) la soglia nella didascalia e' quella delle formals di analyze_list (7, non 10)", {
  skip_if_not(data_ok)
  shiny::testServer(app_dir, {
    # analyze_list sostituita da una versione con min_expected = 7 di default: l'app non passa
    # min_expected, quindi vale il default della firma. Un 10 scritto a mano nel testo darebbe
    # lo stesso risultato di leggere le formals finche' il default e' 10: qui no.
    orig <- app_get("analyze_list", session)
    env <- environment(orig)
    fake <- function(text, genes, hits, frac_range = c(0.005, 0.9), min_expected = 7) {
      orig(text, genes, hits, frac_range, min_expected)
    }
    assign("analyze_list", fake, envir = env)
    on.exit(assign("analyze_list", orig, envir = env), add = TRUE)

    expect_equal(formals(fake)$min_expected, 7)
    session$setInputs(genes = paste(real_symbols(60), collapse = " "), go = 1)
    r <- res()
    expect_null(r$corrected_note)                       # lista normale: la frase precisa
    expect_false(any(is.na(r$table$z) & r$table$expected >= 7))   # e' vera con soglia 7
    s <- box_text(output$caption)
    expect_match(s, "atteso < 7 geni", fixed = TRUE)
    expect_false(grepl("atteso < 10", s, fixed = TRUE))
  })
})

# --- (h) tab "Geni non riconosciuti" e "Come funziona" ----------------------------
test_that("(h) tab Geni non riconosciuti: numero e un token per riga (testo del tutto inerte)", {
  skip_if_not(data_ok)
  shiny::testServer(app_dir, {
    # tre ignoti diversi (uno scritto due volte, uno che sembra HTML)
    txt <- paste(c(real_symbols(60), "FOOBAR", "BAZ", "FOOBAR", "<b>X</b>"), collapse = " ")
    session$setInputs(genes = txt, go = 1)
    r <- res()
    expect_identical(r$not_found, c("FOOBAR", "BAZ", "<b>X</b>"))
    html <- as.character(output$notfound$html)
    expect_match(box_text(output$notfound), "Token non riconosciuti: 3", fixed = TRUE)
    # elenco dentro <pre>: un token per riga, nell'ordine di not_found
    pre <- regmatches(html, regexec("(?s)<pre[^>]*>(.*?)</pre>", html, perl = TRUE))[[1]][2]
    expect_identical(strsplit(unescape_html(pre), "\n")[[1]], r$not_found)
    # il token che sembra HTML e' mostrato come testo, non interpretato
    expect_false(grepl("<b>X</b>", html, fixed = TRUE))
  })
})

test_that("(h) tab Geni non riconosciuti: se non ce ne sono, un messaggio (niente elenco)", {
  skip_if_not(data_ok)
  shiny::testServer(app_dir, {
    session$setInputs(genes = paste(real_symbols(60), collapse = " "), go = 1)
    expect_length(res()$not_found, 0)
    expect_match(box_text(output$notfound), "Nessun token non riconosciuto", fixed = TRUE)
    expect_false(grepl("<pre", as.character(output$notfound$html), fixed = TRUE))
  })
})

test_that("(h) tab Geni non riconosciuti: se l'analisi da' errore, lo dice (non resta vuoto)", {
  skip_if_not(data_ok)
  shiny::testServer(app_dir, {
    session$setInputs(genes = "FOO BAR", go = 1)
    expect_match(box_text(output$notfound), "Nessuna analisi disponibile", fixed = TRUE)
  })
})

test_that("(h) i tre tab ci sono, nell'ordine, e 'Come funziona' e' solo un segnaposto", {
  skip_if_not(data_ok)
  shiny::testServer(app_dir, {
    html <- htmltools::renderTags(app_get("ui", session))$html
    tabs <- unique(gsub('data-value="|"', "", regmatches(html, gregexpr('data-value="[^"]+"', html))[[1]]))
    expect_identical(tabs, c("Analisi", "Geni non riconosciuti", "Come funziona"))
    expect_match(html, "(testo da scrivere)", fixed = TRUE)
  })
})
