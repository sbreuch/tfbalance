# tests/testthat/test-analyze.R
# Test scritti PRIMA del codice (Fase B, passo B1). Fixture sintetica, niente data/*.rds.
# Interfaccia:
#   analyze_list(text, genes, hits, frac_range = c(0.005, 0.9), min_expected = 10)
#     -> list(table, not_found, n_list, n_duplicates, n_motifs,
#             n_excluded_low, n_excluded_high, n_corrected, corrected_note)
#   MAX_LIST = 5000 (costante: massimo di geni riconosciuti, misurato nel Passo 7)

source(testthat::test_path("..", "..", "R", "enrich.R"))
source(testthat::test_path("..", "..", "R", "analyze.R"))

# --- fixture --------------------------------------------------------------
# N geni (G1..GN, Entrez 1001.., Ensembl ENSG00000000001..); un motivo per ogni
# valore di `K` (K = geni con sito). Con N = 200 e le frazioni K/N:
#   K = 0 e' sotto 0.005 (escluso "low"), K = 190 sopra 0.9 (escluso "high"),
#   gli altri sei (10, 20, 30, 40, 60, 100) passano il filtro di default.
# Con una lista di n = 60 geni, E = n * K / N vale 3, 6, 9, 12, 18, 30:
#   solo gli ultimi tre hanno E >= 10 (min_expected di default).
make_toy <- function(N = 200, K = c(0, 10, 20, 30, 40, 60, 100, 190), seed = 1) {
  set.seed(seed)
  genes <- data.frame(
    entrez = as.character(1000L + seq_len(N)),
    symbol = paste0("G", seq_len(N)),
    ensembl = sprintf("ENSG%011d", seq_len(N)),
    gc = runif(N, 0.3, 0.8),
    total_sites = rpois(N, 40),
    total_sites_noext = rpois(N, 30),
    stringsAsFactors = FALSE
  )
  hits <- lapply(seq_along(K), function(j) {
    list(id = sprintf("MA%04d.1", j), tf = paste0("TF", j), L = 8L,
         genes = sort(sample(N, K[j])))
  })
  list(genes = genes, hits = hits)
}

# testo con i simboli G<i> separati da `sep`
sym_text <- function(i, sep = " ") paste0("G", i, collapse = sep)

COLS <- c("id", "tf", "k", "K", "expected", "p_baseline", "fdr_baseline",
          "z", "p_corrected", "fdr_corrected")

# === parsing ================================================================

test_that("separatori: spazio, virgola, punto e virgola e a capo danno lo stesso risultato", {
  toy <- make_toy()
  ref <- analyze_list(sym_text(1:60), toy$genes, toy$hits)
  expect_equal(ref$n_list, 60)
  expect_equal(analyze_list(sym_text(1:60, ","),  toy$genes, toy$hits), ref)
  expect_equal(analyze_list(sym_text(1:60, ";"),  toy$genes, toy$hits), ref)
  expect_equal(analyze_list(sym_text(1:60, "\n"), toy$genes, toy$hits), ref)
  # separatori misti
  mixed <- paste0(sym_text(1:20, " "), ",", sym_text(21:40, ";"), "\n", sym_text(41:60, ", "))
  expect_equal(analyze_list(mixed, toy$genes, toy$hits), ref)
})

test_that("separatore tab (anche misto a spazi) come lo spazio", {
  toy <- make_toy()
  ref <- analyze_list(sym_text(1:60), toy$genes, toy$hits)
  expect_equal(analyze_list(sym_text(1:60, "\t"), toy$genes, toy$hits), ref)
  expect_equal(analyze_list(sym_text(1:60, " \t "), toy$genes, toy$hits), ref)
})

test_that("a capo in stile Windows (\\r\\n) come \\n: nessun \\r resta attaccato ai token", {
  toy <- make_toy()
  ref <- analyze_list(sym_text(1:60), toy$genes, toy$hits)
  r <- analyze_list(paste0(sym_text(1:60, "\r\n"), "\r\n"), toy$genes, toy$hits)
  expect_equal(r, ref)
  expect_length(r$not_found, 0)    # "G5\r" non deve finire in not_found
})

test_that("token vuoti e spazi attorno ai token vengono scartati", {
  toy <- make_toy()
  ref <- analyze_list(sym_text(1:60), toy$genes, toy$hits)
  messy <- paste0("  ,; \n", sym_text(1:60, " ,, ;\n\n "), " ;,\n  ")
  r <- analyze_list(messy, toy$genes, toy$hits)
  expect_equal(r, ref)
  expect_length(r$not_found, 0)    # i vuoti non finiscono in not_found
})

# === token non riconosciuti =================================================

test_that("token non riconosciuti non sono un errore: vanno in not_found, senza ripetizioni", {
  toy <- make_toy()
  txt <- paste(sym_text(1:60), "FOO BAR FOO")
  r <- analyze_list(txt, toy$genes, toy$hits)
  expect_equal(r$n_list, 60)
  expect_length(r$not_found, 2)
  expect_setequal(r$not_found, c("FOO", "BAR"))
  # e il resto del risultato e' quello della lista senza i token sbagliati
  expect_equal(r$table, analyze_list(sym_text(1:60), toy$genes, toy$hits)$table)
})

test_that("senza token non riconosciuti not_found e' character(0)", {
  toy <- make_toy()
  expect_identical(analyze_list(sym_text(1:60), toy$genes, toy$hits)$not_found, character(0))
})

test_that("nessun token riconosciuto: errore italiano chiaro", {
  toy <- make_toy()
  expect_error(analyze_list("FOO BAR", toy$genes, toy$hits), "No gene recognised")
  expect_error(analyze_list("", toy$genes, toy$hits), "No gene recognised")
  expect_error(analyze_list("  ,; \n ", toy$genes, toy$hits), "No gene recognised")
})

# === duplicati dopo la mappatura ============================================

test_that("lo stesso gene scritto in modi diversi conta una volta; n_duplicates conta le ripetizioni", {
  toy <- make_toy()
  ref <- analyze_list(sym_text(1:60), toy$genes, toy$hits)
  expect_equal(ref$n_duplicates, 0)
  # G5 anche come Entrez (1005), G10 anche come Ensembl, G20 scritto due volte: 3 ripetizioni
  txt <- paste(sym_text(1:60), "1005 ENSG00000000010 G20")
  r <- analyze_list(txt, toy$genes, toy$hits)
  expect_equal(r$n_list, 60)
  expect_equal(r$n_duplicates, 3)
  expect_equal(r$table, ref$table)
  expect_length(r$not_found, 0)
})

test_that("un ignoto ripetuto conta 1 in n_duplicates e compare una volta in not_found", {
  toy <- make_toy()
  r <- analyze_list(paste(sym_text(1:60), "FOO FOO"), toy$genes, toy$hits)
  expect_equal(r$n_duplicates, 1)
  expect_identical(r$not_found, "FOO")
  # tre volte = 2 ripetizioni; anche scritto con maiuscole diverse (stessa chiave di map_ids)
  r3 <- analyze_list(paste(sym_text(1:60), "FOO FOO foo"), toy$genes, toy$hits)
  expect_equal(r3$n_duplicates, 2)
  expect_length(r3$not_found, 1)
  # ignoti diversi non sono ripetizioni
  r2 <- analyze_list(paste(sym_text(1:60), "FOO BAR"), toy$genes, toy$hits)
  expect_equal(r2$n_duplicates, 0)
  expect_length(r2$not_found, 2)
})

test_that("un token corrisponde sempre a un solo gene: simbolo, poi Ensembl, poi Entrez; a parita' la prima riga", {
  toy <- make_toy()
  g <- toy$genes
  g$symbol[70] <- "1010"   # il simbolo del gene 70 coincide con l'Entrez del gene 10
  g$symbol[80] <- "G79"    # due righe con lo stesso simbolo
  # "1010" vale il gene 70 (il simbolo vince sull'Entrez), non il 10: lista di 60 geni distinti
  # G1..G59 piu' il gene 70. Riferimento: lo stesso gene scritto come Ensembl.
  ref <- analyze_list(paste(sym_text(1:59), "ENSG00000000070"), g, toy$hits)
  r <- analyze_list(paste(sym_text(1:59), "1010"), g, toy$hits)
  expect_equal(r$n_list, 60)            # se valesse il gene 10 sarebbe 59 (gia' in lista)
  expect_equal(r$n_duplicates, 0)
  expect_equal(r, ref)
  # "G79" con due righe: vale la prima (riga 79), una sola volta
  r79 <- analyze_list(paste(sym_text(1:58), "G79"), g, toy$hits)
  expect_equal(r79, analyze_list(paste(sym_text(1:58), "ENSG00000000079"), g, toy$hits))
  expect_equal(r79$n_list, 59)
})

# === errori di check_list e limite MAX_LIST =================================

test_that("gli errori di check_list si propagano con il loro messaggio", {
  toy <- make_toy()
  expect_error(analyze_list(sym_text(1:200), toy$genes, toy$hits), "all genes")
})

test_that("MAX_LIST e' 5000", {
  expect_equal(MAX_LIST, 5000)
})

test_that("fino a MAX_LIST geni riconosciuti va bene; non contano ne' i non riconosciuti ne' i duplicati", {
  toy <- make_toy(N = 6000, K = c(40, 40, 40))
  txt <- paste(sym_text(seq_len(MAX_LIST)), "NOTAGENE1 NOTAGENE2 G1 G2 G3")
  r <- analyze_list(txt, toy$genes, toy$hits)
  expect_equal(r$n_list, MAX_LIST)
  expect_equal(r$n_duplicates, 3)
  expect_length(r$not_found, 2)
})

test_that("oltre MAX_LIST geni riconosciuti: errore chiaro che cita il limite", {
  toy <- make_toy(N = 6000, K = c(40, 40, 40))
  expect_error(analyze_list(sym_text(seq_len(MAX_LIST + 1)), toy$genes, toy$hits),
               "maximum is 5000")
})

test_that("MAX_TOKENS e' 50000", {
  expect_equal(MAX_TOKENS, 50000)
})

test_that("fino a MAX_TOKENS token grezzi va bene (anche con molte ripetizioni)", {
  toy <- make_toy()
  txt <- paste(rep_len(paste0("G", 1:60), MAX_TOKENS), collapse = " ")
  r <- analyze_list(txt, toy$genes, toy$hits)
  expect_equal(r$n_list, 60)
  expect_equal(r$n_duplicates, MAX_TOKENS - 60)
})

test_that("oltre MAX_TOKENS token grezzi: errore chiaro che cita il limite, PRIMA di map_ids", {
  toy <- make_toy()
  txt <- paste(rep("G1", MAX_TOKENS + 1), collapse = " ")
  expect_error(analyze_list(txt, toy$genes, toy$hits), "maximum is 50000")
  # map_ids sostituita da una funzione che fallisce con un altro messaggio: se venisse
  # chiamata prima del controllo, "50000" non comparirebbe nell'errore
  env <- environment(analyze_list)
  orig <- get("map_ids", envir = env)
  assign("map_ids", function(...) stop("map_ids chiamata"), envir = env)
  on.exit(assign("map_ids", orig, envir = env), add = TRUE)
  expect_error(analyze_list(txt, toy$genes, toy$hits), "maximum is 50000")
})

# === risultato ==============================================================

test_that("il risultato e' una lista con i campi previsti", {
  toy <- make_toy()
  r <- analyze_list(sym_text(1:60), toy$genes, toy$hits)
  expect_type(r, "list")
  expect_setequal(names(r), c("table", "not_found", "n_list", "n_duplicates", "n_motifs",
                              "n_excluded_low", "n_excluded_high", "n_corrected",
                              "corrected_note"))
  expect_s3_class(r$table, "data.frame")
  expect_equal(r$n_list, 60)
})

test_that("table ha le colonne previste, una riga per motivo tenuto dal filtro", {
  toy <- make_toy()
  r <- analyze_list(sym_text(1:60), toy$genes, toy$hits)
  expect_identical(names(r$table), COLS)
  expect_equal(nrow(r$table), r$n_motifs)
  expect_equal(r$n_motifs, 6)
  # i motivi esclusi (K = 0 e K = 190, cioe' il primo e l'ultimo) non compaiono
  expect_false(any(c("MA0001.1", "MA0008.1") %in% r$table$id))
  expect_setequal(r$table$id, sprintf("MA%04d.1", 2:7))
})

test_that("i valori coincidono con enrich_baseline e enrich_corrected sui motivi tenuti", {
  toy <- make_toy()
  N <- nrow(toy$genes)
  r <- analyze_list(sym_text(1:60), toy$genes, toy$hits)
  keep <- select_motifs(toy$hits, N)
  b <- enrich_baseline(1:60, toy$hits[keep], N)
  cc <- enrich_corrected(1:60, toy$hits[keep], toy$genes)
  tab <- r$table[match(b$id, r$table$id), ]     # stesso ordine dei riferimenti
  expect_equal(tab$tf, b$tf)
  expect_equal(tab$k, b$k)
  expect_equal(tab$K, b$K)
  expect_equal(tab$expected, b$expected)
  expect_equal(tab$p_baseline, b$p)
  expect_equal(tab$fdr_baseline, b$fdr)     # BH solo sui motivi tenuti
  expect_equal(tab$z, cc$z)
  expect_equal(tab$p_corrected, cc$p)
  expect_equal(tab$fdr_corrected, cc$fdr)
})

test_that("table e' ordinata per p_corrected crescente, NA in fondo, spareggio su p_baseline", {
  toy <- make_toy()
  tab <- analyze_list(sym_text(1:60), toy$genes, toy$hits)$table
  na <- is.na(tab$p_corrected)
  expect_equal(sum(na), 3)                         # E = 3, 6, 9 < 10
  expect_equal(which(na), 4:6)                     # gli NA stanno tutti in fondo
  expect_false(is.unsorted(tab$p_corrected[!na]))
  expect_false(is.unsorted(tab$p_baseline[na]))    # tra gli NA decide p_baseline
  expect_gt(length(unique(tab$p_baseline[na])), 1) # (il test non e' vuoto)
})

# === filtro e conteggi ======================================================

test_that("conteggi del filtro: n_motifs, n_excluded_low, n_excluded_high", {
  toy <- make_toy()
  r <- analyze_list(sym_text(1:60), toy$genes, toy$hits)
  expect_equal(r$n_excluded_low, 1)
  expect_equal(r$n_excluded_high, 1)
  expect_equal(r$n_motifs, 6)
  expect_equal(r$n_motifs + r$n_excluded_low + r$n_excluded_high, length(toy$hits))
})

test_that("frac_range = c(0, 1): nessun motivo escluso", {
  toy <- make_toy()
  r <- analyze_list(sym_text(1:60), toy$genes, toy$hits, frac_range = c(0, 1))
  expect_equal(r$n_excluded_low, 0)
  expect_equal(r$n_excluded_high, 0)
  expect_equal(r$n_motifs, length(toy$hits))
  expect_equal(nrow(r$table), length(toy$hits))
})

test_that("il filtro esclude tutti i motivi: errore italiano che nasce in analyze_list", {
  toy <- make_toy()
  # frazioni dei motivi: 0, .05, .1, .15, .2, .3, .5, .95; nessuna in [0.6, 0.9]
  expect_error(analyze_list(sym_text(1:60), toy$genes, toy$hits, frac_range = c(0.6, 0.9)),
               "excludes all motifs")
  # l'errore deve venire PRIMA di enrich_baseline: la sostituisco con una funzione che
  # fallisce con un altro messaggio (se fosse chiamata, "filtro" non comparirebbe)
  env <- environment(analyze_list)
  orig <- get("enrich_baseline", envir = env)
  assign("enrich_baseline", function(...) stop("enrich_baseline chiamata"), envir = env)
  on.exit(assign("enrich_baseline", orig, envir = env), add = TRUE)
  expect_error(analyze_list(sym_text(1:60), toy$genes, toy$hits, frac_range = c(0.6, 0.9)),
               "excludes all motifs")
})

test_that("frac_range e min_expected non validi: l'errore di enrich.R arriva all'utente", {
  toy <- make_toy()
  expect_error(analyze_list(sym_text(1:60), toy$genes, toy$hits, frac_range = c(0.9, 0.1)),
               "must be c\\(minimum, maximum\\)")
  expect_error(analyze_list(sym_text(1:60), toy$genes, toy$hits, min_expected = -1),
               "min_expected.*single finite number")
})

# === test corretto: n_corrected e corrected_note ============================

test_that("n_corrected conta i motivi con z non NA; corrected_note e' NULL se > 0", {
  toy <- make_toy()
  r <- analyze_list(sym_text(1:60), toy$genes, toy$hits)
  expect_equal(r$n_corrected, 3)
  expect_equal(r$n_corrected, sum(!is.na(r$table$z)))
  expect_null(r$corrected_note)
})

test_that("min_expected = 0: tutti i motivi tenuti hanno il test corretto", {
  toy <- make_toy()
  r <- analyze_list(sym_text(1:60), toy$genes, toy$hits, min_expected = 0)
  expect_equal(r$n_corrected, r$n_motifs)
  expect_null(r$corrected_note)
})

test_that("nessun motivo con E >= min_expected: n_corrected = 0 e corrected_note spiega perche'", {
  toy <- make_toy()
  r <- analyze_list(sym_text(1:10), toy$genes, toy$hits)    # E = 10 * K / 200 <= 5
  expect_equal(r$n_corrected, 0)
  expect_true(all(is.na(r$table$z)))
  expect_true(all(is.na(r$table$p_corrected)))
  expect_false(anyNA(r$table$p_baseline))                   # il baseline c'e' comunque
  expect_false(is.unsorted(r$table$p_baseline))             # ordine = p_baseline
  expect_type(r$corrected_note, "character")
  expect_length(r$corrected_note, 1)
  expect_match(r$corrected_note, "no motif has E")
  # la nota cita la soglia usata
  r12 <- analyze_list(sym_text(1:10), toy$genes, toy$hits, min_expected = 12)
  expect_match(r12$corrected_note, "12")
})

# === modello non stimabile ==================================================

test_that("modello non stimabile (enrich_corrected sostituita: tutti NA): la nota non parla di E", {
  toy <- make_toy()
  env <- environment(analyze_list)
  orig <- get("enrich_corrected", envir = env)
  # finta enrich_corrected: stesso avviso e stesse colonne dell'originale, z/p/fdr tutti NA,
  # ma `expected` vero: con n = 60 tre motivi hanno E >= 10
  assign("enrich_corrected", function(idx, hits, genes, min_expected = 10) {
    warning("the null model cannot be fitted (mock)", call. = FALSE)
    K <- lengths(lapply(hits, `[[`, "genes"))
    data.frame(id = vapply(hits, function(h) h$id, ""),
               tf = vapply(hits, function(h) h$tf, ""),
               expected = length(idx) * K / nrow(genes),
               z = NA_real_, p = NA_real_, fdr = NA_real_, stringsAsFactors = FALSE)
  }, envir = env)
  on.exit(assign("enrich_corrected", orig, envir = env), add = TRUE)

  # l'avviso di enrich_corrected arriva al chiamante: analyze_list non lo intercetta
  expect_warning(r <- analyze_list(sym_text(1:60), toy$genes, toy$hits), "null model")
  expect_equal(r$n_corrected, 0)
  expect_true(any(r$table$expected >= 10))     # quindi NON e' un problema di E bassa
  expect_match(r$corrected_note, "cannot be fitted")
  expect_false(grepl("no motif has E", r$corrected_note))
  expect_false(is.unsorted(r$table$p_baseline))   # tutti NA: ordine = p_baseline
})

test_that("modello non stimabile con input reale (lista = i geni con GC piu' alto: separazione)", {
  toy <- make_toy()
  top <- order(toy$genes$gc, decreasing = TRUE)[1:60]
  txt <- paste(toy$genes$symbol[top], collapse = " ")
  expect_warning(r <- analyze_list(txt, toy$genes, toy$hits), "null model")
  expect_equal(r$n_corrected, 0)
  expect_match(r$corrected_note, "cannot be fitted")
  expect_false(anyNA(r$table$p_baseline))      # il baseline c'e' comunque
})
