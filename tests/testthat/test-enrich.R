# tests/testthat/test-enrich.R
# Test scritti PRIMA del codice. Usano una fixture sintetica piccola (non data/*.rds).
# Interfaccia (piano approvato):
#   map_ids(x, genes)                     -> list(idx = integer, not_found = character)
#   select_motifs(hits, n_genes, frac_range = c(0.005, 0.9)) -> logical
#   enrich_baseline(idx, hits, n_genes)   -> data.frame(id, tf, k, K, expected, p, fdr)
#   enrich_corrected(idx, hits, genes)    -> data.frame(id, tf, z, p, fdr)   # test di punteggio
#   odds_ratio_one(idx, motif_genes, genes) -> data.frame(beta, se, or, or_low, or_high, note)
#   random_list(genes, n, biased = FALSE, by = "gc", power = 10) -> indici di riga

source(testthat::test_path("..", "..", "R", "enrich.R"))

# --- fixture --------------------------------------------------------------
# N geni con colonne come genes.rds; M motivi con "genes" = indici di riga.
make_toy <- function(N = 200, M = 8, seed = 1) {
  set.seed(seed)
  genes <- data.frame(
    entrez = as.character(1000L + seq_len(N)),   # come genes.rds: testo, non numero
    symbol = paste0("G", seq_len(N)),
    ensembl = sprintf("ENSG%011d", seq_len(N)),
    gc = runif(N, 0.3, 0.8),
    total_sites = rpois(N, 40),
    total_sites_noext = rpois(N, 30),
    stringsAsFactors = FALSE
  )
  hits <- lapply(seq_len(M), function(j) {
    list(id = sprintf("MA%04d.1", j), tf = paste0("TF", j), L = 8L,
         genes = sort(sample(N, 40)))
  })
  list(genes = genes, hits = hits)
}

# === map_ids ================================================================

test_that("map_ids riconosce simbolo, Entrez ed Ensembl misti", {
  toy <- make_toy()
  r <- map_ids(c("G5", "1010", "ENSG00000000020"), toy$genes)
  expect_equal(sort(r$idx), c(5, 10, 20))
  expect_length(r$not_found, 0)
})

test_that("map_ids accetta Entrez numerici", {
  toy <- make_toy()
  expect_equal(sort(map_ids(c(1010, 1020), toy$genes)$idx), c(10, 20))
})

test_that("map_ids segnala i nomi inventati e non li include in idx", {
  toy <- make_toy()
  r <- map_ids(c("G5", "NOTAGENE"), toy$genes)
  expect_equal(r$idx, 5)
  expect_equal(r$not_found, "NOTAGENE")
})

test_that("map_ids non duplica lo stesso gene scritto in due modi", {
  toy <- make_toy()
  r <- map_ids(c("G5", "G5", "1005"), toy$genes)
  expect_equal(r$idx, 5)
})

test_that("map_ids ignora spazi, righe vuote e NA", {
  toy <- make_toy()
  r <- map_ids(c(" G5 ", "", NA), toy$genes)
  expect_equal(r$idx, 5)
  expect_length(r$not_found, 0)
})

test_that("map_ids: simboli case-insensitive e Ensembl con suffisso di versione", {
  toy <- make_toy()
  expect_equal(map_ids("g5", toy$genes)$idx, 5)
  expect_equal(map_ids("ENSG00000000020.7", toy$genes)$idx, 20)
})

test_that("map_ids: input vuoto -> nessun indice e nessun non riconosciuto", {
  toy <- make_toy()
  r <- map_ids(character(0), toy$genes)
  expect_length(r$idx, 0)
  expect_length(r$not_found, 0)
})

test_that("map_ids: lo stesso nome sbagliato scritto in modi diversi compare una volta", {
  toy <- make_toy()
  r <- map_ids(c("foo", "FOO", " foo "), toy$genes)
  expect_equal(r$not_found, "foo")
})

test_that("map_ids: Ensembl mancante (NA) in genes non crea falsi match", {
  toy <- make_toy()
  # come i 7.752 geni reali senza Ensembl; conteggio:
  #   sum(is.na(readRDS("data/genes.rds")$ensembl))   # = 7752, nessuna stringa vuota
  toy$genes$ensembl[5] <- NA
  expect_equal(map_ids("G5", toy$genes)$idx, 5)
  expect_equal(map_ids("NA", toy$genes)$not_found, "NA")
})

test_that("map_ids: Entrez numerico grande non diventa notazione scientifica", {
  toy <- make_toy()
  toy$genes$entrez[7] <- "100000"                 # as.character(100000) darebbe "1e+05"
  expect_equal(map_ids(100000, toy$genes)$idx, 7)
})

test_that("map_ids: se mancano colonne in genes e' un errore", {
  toy <- make_toy()
  expect_error(map_ids("G5", toy$genes[, c("symbol", "gc")]), "deve avere le colonne")
})

# === select_motifs ==========================================================

mk_hits <- function(sizes) lapply(sizes, function(K) list(id = "x", tf = "x", genes = seq_len(K)))

test_that("select_motifs: tiene frac tra 0.005 e 0.9, estremi inclusi", {
  # N = 1000: K = 4 -> 0.004 (fuori), 5 -> 0.005 (dentro), 900 -> 0.9 (dentro), 901 (fuori)
  hits <- mk_hits(c(0, 4, 5, 500, 900, 901, 1000))
  expect_equal(select_motifs(hits, 1000), c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE))
})

test_that("select_motifs: frac_range = c(0, 1) tiene tutti i motivi", {
  expect_true(all(select_motifs(mk_hits(c(0, 4, 1000)), 1000, c(0, 1))))
})

test_that("select_motifs: frac_range non valido e' un errore", {
  hits <- mk_hits(c(10, 20))
  expect_error(select_motifs(hits, 1000, 0.5), "frac_range")
  expect_error(select_motifs(hits, 1000, c(0.9, 0.1)), "frac_range")
})

# === enrich_baseline ========================================================

test_that("enrich_baseline coincide con fisher.test su tabella 2x2 fatta a mano", {
  toy <- make_toy()
  N <- nrow(toy$genes)
  set.seed(2)
  idx <- sample(N, 30)
  res <- enrich_baseline(idx, toy$hits, N)
  expect_equal(nrow(res), length(toy$hits))
  for (j in seq_along(toy$hits)) {
    K <- length(toy$hits[[j]]$genes)           # geni con sito nell'universo
    k <- length(intersect(idx, toy$hits[[j]]$genes))  # geni con sito nella lista
    n <- length(idx)
    # righe: in lista / fuori; colonne: con sito / senza sito
    tab <- matrix(c(k, K - k, n - k, N - K - (n - k)), nrow = 2)
    p_ref <- fisher.test(tab, alternative = "greater")$p.value
    expect_equal(res$p[j], p_ref, tolerance = 1e-10)
    expect_equal(res$k[j], k)
    expect_equal(res$K[j], K)
  }
})

test_that("enrich_baseline: attesi = n * K / N", {
  toy <- make_toy()
  N <- nrow(toy$genes)
  idx <- 1:30
  res <- enrich_baseline(idx, toy$hits, N)
  expect_equal(res$expected, vapply(toy$hits, function(h) 30 * length(h$genes) / N, 0))
})

test_that("enrich_baseline: fdr = Benjamini-Hochberg fatto a mano, su tutti i motivi passati", {
  toy <- make_toy()
  N <- nrow(toy$genes)
  toy$hits[[1]]$genes <- integer(0)      # K = 0: p = 1, ma conta tra i motivi testati
  toy$hits[[2]]$genes <- seq_len(N)      # K = N: p = 1
  res <- enrich_baseline(1:30, toy$hits, N)
  expect_equal(res$p[1:2], c(1, 1))
  m <- length(res$p)                     # BH: m = numero di motivi testati
  o <- order(res$p, decreasing = TRUE)   # dal p piu' grande al piu' piccolo
  fdr_hand <- numeric(m)
  fdr_hand[o] <- cummin(res$p[o] * m / (m:1))   # p * m / rango, poi minimo cumulato
  expect_equal(res$fdr, pmin(fdr_hand, 1))
  expect_true(all(res$fdr >= res$p))
})

test_that("enrich_baseline: nessun gene della lista col sito -> p = 1", {
  toy <- make_toy()
  toy$hits[[1]]$genes <- 101:140
  res <- enrich_baseline(1:30, toy$hits, nrow(toy$genes))
  expect_equal(res$p[1], 1)
})

test_that("enrich_baseline: liste non valide sono errori con messaggio chiaro", {
  toy <- make_toy()
  N <- nrow(toy$genes)
  expect_error(enrich_baseline(integer(0), toy$hits, N), "vuota")
  expect_error(enrich_baseline(seq_len(N), toy$hits, N), "tutti i geni")
  expect_error(enrich_baseline(c(1, N + 1), toy$hits, N), "tra 1 e")
  expect_error(enrich_baseline(c(1, NA), toy$hits, N), "tra 1 e")
  expect_error(enrich_baseline(c(1, 1, 2), toy$hits, N), "ripetuti")
})

# === enrich_corrected =======================================================

test_that("enrich_corrected trova un segnale forte e piantato", {
  toy <- make_toy(N = 300)
  toy$hits[[1]]$genes <- 1:80                    # motivo A: geni 1..80
  # lista: 50 geni con A + 10 senza. Se tutti i geni senza A fossero fuori lista
  # ci sarebbe separazione perfetta e il p-value di Wald crollerebbe (vedi sotto).
  idx <- c(1:50, 101:110)
  res <- enrich_corrected(idx, toy$hits, toy$genes)
  expect_equal(nrow(res), length(toy$hits))
  expect_equal(which.min(res$p), 1)
  expect_lt(res$p[1], 1e-6)
  expect_gt(res$or[1], 1)
  expect_equal(res$or[1], exp(res$beta[1]))
})

test_that("enrich_corrected: fdr = Benjamini-Hochberg dei p (NA esclusi)", {
  toy <- make_toy(N = 300)
  res <- enrich_corrected(1:60, toy$hits, toy$genes)
  expect_equal(res$fdr, p.adjust(res$p, "BH"))
})

# Il test che conta: confondimento da GC. Il motivo B non ha alcun effetto
# *dato il GC*, ma e' piu' frequente nei geni ad alto GC, e la lista e'
# sbilanciata verso alto GC. Il baseline lo vede arricchito, il corretto no.
# ATTENZIONE: test statistico con seed; non scritto con R a disposizione.
# Se fallisce, ragionaci sopra: non adattare soglie per farlo passare.
test_that("enrich_corrected riduce il falso positivo da confondimento GC", {
  set.seed(42)
  N <- 2000
  gc <- seq(0.3, 0.8, length.out = N)
  genes <- data.frame(
    entrez = 1000L + seq_len(N), symbol = paste0("G", seq_len(N)),
    ensembl = sprintf("ENSG%011d", seq_len(N)), gc = gc,
    total_sites = rpois(N, 40), total_sites_noext = rpois(N, 30),
    stringsAsFactors = FALSE
  )
  p_gc <- plogis(12 * (gc - 0.55))
  B <- which(rbinom(N, 1, p_gc) == 1)            # presenza di B dipende dal GC
  idx <- which(rbinom(N, 1, p_gc) == 1)          # lista dipende dal GC, non da B
  hits <- list(list(id = "B", tf = "TFB", L = 8L, genes = B))
  base <- enrich_baseline(idx, hits, N)
  corr <- enrich_corrected(idx, hits, genes)
  expect_lt(base$p, 1e-6)
  expect_gt(corr$p, 1e-3)
  expect_gt(corr$p, base$p)
})

test_that("enrich_corrected: motivo senza siti nella lista non va in errore", {
  toy <- make_toy(N = 300)
  toy$hits[[1]]$genes <- 201:240                 # nessun gene della lista 1:60
  res <- enrich_corrected(1:60, toy$hits, toy$genes)
  expect_equal(nrow(res), length(toy$hits))
  expect_true(is.na(res$p[1]) || res$p[1] > 0.5)
})

test_that("enrich_corrected: motivo presente in tutti i geni -> p = NA, nessun errore", {
  toy <- make_toy(N = 300)
  toy$hits[[1]]$genes <- seq_len(300)
  res <- enrich_corrected(1:60, toy$hits, toy$genes)
  expect_equal(nrow(res), length(toy$hits))
  expect_true(is.na(res$p[1]))
})

# === random_list ============================================================

test_that("random_list: n indici validi e distinti", {
  toy <- make_toy(N = 500)
  set.seed(3)
  x <- random_list(toy$genes, 50)
  expect_length(x, 50)
  expect_false(anyDuplicated(x) > 0)
  expect_true(all(x >= 1 & x <= 500))
})

test_that("random_list: stesso seed, stessa lista", {
  toy <- make_toy(N = 500)
  set.seed(4); a <- random_list(toy$genes, 50, biased = TRUE)
  set.seed(4); b <- random_list(toy$genes, 50, biased = TRUE)
  expect_equal(a, b)
})

test_that("random_list biased: GC medio piu' alto dell'uniforme", {
  toy <- make_toy(N = 2000)
  set.seed(5)
  m_u <- mean(replicate(20, mean(toy$genes$gc[random_list(toy$genes, 300)])))
  m_b <- mean(replicate(20, mean(toy$genes$gc[random_list(toy$genes, 300, biased = TRUE)])))
  expect_gt(m_b, m_u + 0.03)
})

# Dichiarato prima di eseguire: con power < 0 un gene a 0 pesa 0^(-10) = Inf se
# non gestito. Politica: lo 0 vale come il piu' piccolo valore positivo, quindi
# e' tra i geni piu' favoriti ma con peso finito. Atteso: scelto in >= 45 liste su 50.
test_that("random_list power < 0: un gene con valore 0 non da' peso infinito", {
  toy <- make_toy(N = 500)
  toy$genes$total_sites_noext[1] <- 0L
  set.seed(7)
  lists <- replicate(50, random_list(toy$genes, 50, biased = TRUE,
                                     by = "total_sites_noext", power = -10),
                     simplify = FALSE)
  valid <- vapply(lists, function(x) length(x) == 50 && !anyDuplicated(x) &&
                    all(x >= 1 & x <= 500), TRUE)
  expect_true(all(valid))
  expect_gte(sum(vapply(lists, function(x) 1 %in% x, TRUE)), 45)
})

test_that("random_list biased by total_sites_noext: densita' media piu' alta", {
  toy <- make_toy(N = 2000)
  set.seed(6)
  d <- toy$genes$total_sites_noext
  m_u <- mean(replicate(20, mean(d[random_list(toy$genes, 300)])))
  m_b <- mean(replicate(20, mean(d[random_list(toy$genes, 300, biased = TRUE,
                                              by = "total_sites_noext")])))
  expect_gt(m_b, m_u)
})
