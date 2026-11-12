# tests/testthat/test-enrich.R
# Test scritti PRIMA del codice. Usano una fixture sintetica piccola (non data/*.rds).
# Interfaccia (piano approvato):
#   map_ids(x, genes)                     -> list(idx = integer, not_found = character)
#   select_motifs(hits, n_genes, frac_range = c(0.005, 0.9)) -> logical
#   enrich_baseline(idx, hits, n_genes)   -> data.frame(id, tf, k, K, expected, p, fdr)
#   enrich_corrected(idx, hits, genes)    -> data.frame(id, tf, z, p, fdr)   # test di punteggio
#   odds_ratio_one(idx, motif_genes, genes) -> data.frame(beta, se, or, or_low, or_high, note)
#   random_list(genes, n, biased = FALSE, by = "gc", power = 2) -> indici di riga, ordinati

source(testthat::test_path("..", "..", "R", "enrich.R"))

# --- fixture --------------------------------------------------------------
# N geni con colonne come genes.rds; M motivi con "genes" = indici di riga.
make_toy <- function(N = 200, M = 8, seed = 1, K = 40) {   # K = geni con sito per motivo
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
         genes = sort(sample(N, K)))
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
# Test di punteggio (score) sul modello in_list ~ has_site + gc + log1p(total_sites_noext).
#
# GUARDIA SULLA VARIANZA DELLO SCORE (dichiarata prima del codice)
# Per un motivo, z = U / sqrt(v), con v = informazione di `has_site` DOPO aver tolto
# quella gia' spiegata da intercetta, gc e densita'. v e' una differenza di due numeri
# grandi: se `has_site` e' (quasi) l'intercetta (K = N) v e' solo rumore di
# arrotondamento e z = U / sqrt(rumore) sarebbe un numero casuale. Il motivo vale NA se
#     K == 0   oppure   K == N   oppure   v / sum(w) < SCORE_VAR_MIN,   SCORE_VAR_MIN = 1e-8
# dove sum(w) e' l'informazione totale del modello nullo (w = mu * (1 - mu)).
# Perche' RELATIVA: v scala con la dimensione della lista (w e' minuscolo per liste
# piccole), una soglia assoluta sarebbe arbitraria. Perche' 1e-8: l'errore di arrotondamento
# e' ~1e-16 * sum(w), quindi 1e-8 sta 8 ordini sopra il rumore; e sta molto sotto
# l'informazione vera anche di un solo gene (K = N - 1 o K = 1: v ~ w di quel gene,
# v / sum(w) ~ 1e-3 .. 1e-5). K = N - 1 NON e' scartato: e' un motivo quasi costante
# ma informativo. Il test "la soglia sta nel divario" lo verifica sull'implementazione.

# BH a mano: gli NA restano NA e non contano in m (numero di test validi)
bh_hand <- function(p) {
  ok <- !is.na(p); m <- sum(ok); out <- rep(NA_real_, length(p))
  q <- p[ok]; o <- order(q, decreasing = TRUE)
  f <- numeric(m); f[o] <- cummin(q[o] * m / (m:1))
  out[ok] <- pmin(f, 1)
  out
}

test_that("enrich_corrected trova un segnale forte e piantato", {
  toy <- make_toy(N = 300)
  toy$hits[[1]]$genes <- 1:80                    # motivo A: geni 1..80
  idx <- c(1:50, 101:110)                        # lista: 50 geni con A + 10 senza
  res <- enrich_corrected(idx, toy$hits, toy$genes, min_expected = 0)
  expect_named(res, c("id", "tf", "expected", "z", "p", "fdr"))
  expect_equal(res$id, vapply(toy$hits, function(h) h$id, ""))
  expect_equal(which.min(res$p), 1)
  expect_lt(res$p[1], 1e-6)
  expect_gt(res$z[1], 0)
})

test_that("enrich_corrected: fdr = BH fatto a mano; gli NA non contano in m", {
  toy <- make_toy(N = 300)
  toy$hits[[1]]$genes <- seq_len(300)            # K = N -> p = NA
  res <- enrich_corrected(1:60, toy$hits, toy$genes, min_expected = 0)
  expect_true(is.na(res$p[1]) && is.na(res$fdr[1]))
  expect_equal(res$fdr, bh_hand(res$p))
  expect_equal(sum(!is.na(res$fdr)), length(toy$hits) - 1)
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

test_that("enrich_corrected: motivo senza geni della lista -> z negativo, nessun errore", {
  toy <- make_toy(N = 300)
  toy$hits[[1]]$genes <- 201:240                 # nessun gene della lista 1:60
  res <- enrich_corrected(1:60, toy$hits, toy$genes, min_expected = 0)
  expect_equal(nrow(res), length(toy$hits))
  expect_lt(res$z[1], 0)
  expect_gt(res$p[1], 0.5)
})

test_that("enrich_corrected: K = 0 e K = N danno NA; K = N - 1 e K = 1 un z finito", {
  toy <- make_toy(N = 300)
  toy$hits[[1]]$genes <- integer(0)              # K = 0
  toy$hits[[2]]$genes <- seq_len(300)            # K = N
  toy$hits[[3]]$genes <- 2:300                   # K = N - 1: manca il gene 1 (in lista)
  toy$hits[[4]]$genes <- 1:299                   # K = N - 1: manca il gene 300 (fuori lista)
  toy$hits[[5]]$genes <- 7L                      # K = 1
  res <- enrich_corrected(1:60, toy$hits, toy$genes, min_expected = 0)
  expect_true(all(is.na(res$z[1:2]) & is.na(res$p[1:2]) & is.na(res$fdr[1:2])))
  expect_true(all(is.finite(res$z[3:5])))
  expect_true(all(is.finite(res$p[3:5])))
})

# Lista interamente dentro il motivo: tutti i 60 geni della lista hanno il sito e
# nessuno ne e' privo (il motivo ne ha anche altri 20 fuori lista, quindi K = 80 < N).
# Atteso, dichiarato prima di eseguire: z finito, z > 0, p < 1e-6.
test_that("enrich_corrected: lista tutta dentro il motivo -> z finito e positivo, p < 1e-6", {
  toy <- make_toy(N = 300)
  toy$hits[[1]]$genes <- 1:80                    # motivo: geni 1..80
  idx <- 1:60                                    # lista: tutti con il sito, nessuno senza
  res <- enrich_corrected(idx, toy$hits, toy$genes)
  expect_true(is.finite(res$z[1]))
  expect_gt(res$z[1], 0)
  expect_lt(res$p[1], 1e-6)
})

# Motivo = esattamente la lista (K = n = 60): has_site coincide con in_list, cioe'
# separazione perfetta. Per il Wald di glm e' un caso patologico; lo score usa il
# modello nullo (senza has_site), quindi deve dare un valore.
# Atteso, dichiarato prima di eseguire: z finito, z > 0, p < 1e-6.
test_that("enrich_corrected: motivo = esattamente la lista (K = n) -> z finito e positivo, p < 1e-6", {
  toy <- make_toy(N = 300)
  toy$hits[[1]]$genes <- 1:60                    # il motivo coincide con la lista
  res <- enrich_corrected(1:60, toy$hits, toy$genes)
  expect_true(is.finite(res$z[1]))
  expect_gt(res$z[1], 0)
  expect_lt(res$p[1], 1e-6)
})

test_that("enrich_corrected: liste e covariate non valide sono errori chiari", {
  toy <- make_toy(N = 300)
  expect_error(enrich_corrected(integer(0), toy$hits, toy$genes), "vuota")
  expect_error(enrich_corrected(c(1, 1, 2), toy$hits, toy$genes), "ripetuti")
  g2 <- toy$genes; g2$gc <- 0.5                  # gc costante: non si puo' standardizzare
  expect_error(enrich_corrected(1:60, toy$hits, g2), "non costanti")
})

test_that("enrich_corrected: modello nullo non stimabile -> tutti NA e un avviso", {
  toy <- make_toy(N = 300)
  toy$genes$gc <- sort(toy$genes$gc)             # gc cresce col numero di riga
  idx <- 296:300                                 # i 5 geni a gc piu' alto: gc li separa da soli
  expect_warning(res <- enrich_corrected(idx, toy$hits, toy$genes, min_expected = 0), "non e' stimabile")
  expect_true(all(is.na(res$z)))
})

# === min_expected: soglia sui geni attesi E ==================================
# Regola dichiarata prima del codice (esperimento H2 confermato: il test corretto e'
# tarato per E >= 10). E = n * K / N, calcolato moltiplicando prima di dividere (come
# `expected` in enrich_baseline), con n = geni della lista, K = geni con il sito,
# N = tutti i geni. Il motivo con E < min_expected (default 10) ha z, p, fdr = NA ma
# RESTA nelle righe; E == min_expected e' incluso (E >= min_expected). La colonna
# `expected` contiene E ed e' sempre riempita. Il BH e' calcolato solo sui motivi tenuti:
# gli NA non contano in m. min_expected deve essere un numero singolo, finito, >= 0
# (quindi NA e Inf sono errori).
# Fixture: N = 300 e lista 1:60, quindi E = K / 5.

test_that("enrich_corrected: min_expected vale 10 di default e `expected` = n*K/N", {
  expect_equal(formals(enrich_corrected)$min_expected, 10)
  toy <- make_toy(N = 300)
  res <- enrich_corrected(1:60, toy$hits, toy$genes)       # K = 40 per tutti: E = 8
  expect_equal(res$expected, rep(60 * 40 / 300, length(toy$hits)))
})

test_that("enrich_corrected: E < min_expected -> NA ma la riga resta; E = 10 esatto e' incluso", {
  toy <- make_toy(N = 300)
  toy$hits[[1]]$genes <- 101:149                 # K = 49 -> E = 9.8  (escluso)
  toy$hits[[2]]$genes <- 101:150                 # K = 50 -> E = 10   (incluso: E >= 10)
  toy$hits[[3]]$genes <- 101:151                 # K = 51 -> E = 10.2 (incluso)
  res <- enrich_corrected(1:60, toy$hits, toy$genes)
  expect_equal(nrow(res), length(toy$hits))      # nessuna riga persa
  expect_equal(res$expected[1:3], c(9.8, 10, 10.2))
  expect_true(is.na(res$z[1]) && is.na(res$p[1]) && is.na(res$fdr[1]))
  expect_true(all(is.finite(res$z[2:3])) && all(is.finite(res$p[2:3])))
  expect_true(all(is.na(res$z[4:8])))            # gli altri: K = 40, E = 8
})

test_that("enrich_corrected: min_expected e' un parametro (con 12 esce E = 11, resta E = 12)", {
  toy <- make_toy(N = 300)
  toy$hits[[1]]$genes <- 101:155                 # K = 55 -> E = 11
  toy$hits[[2]]$genes <- 101:160                 # K = 60 -> E = 12 (incluso)
  res <- enrich_corrected(1:60, toy$hits, toy$genes, min_expected = 12)
  expect_true(is.na(res$z[1]))
  expect_true(is.finite(res$z[2]))
})

test_that("enrich_corrected: il BH e' calcolato solo sui motivi tenuti (gli NA non contano in m)", {
  toy <- make_toy(N = 300)
  toy$hits[[1]]$genes <- 1:80                    # E = 16, segnale forte
  toy$hits[[2]]$genes <- 1:55                    # E = 11, segnale forte
  toy$hits[[3]]$genes <- 101:160                 # E = 12, nessun gene della lista
  res <- enrich_corrected(1:60, toy$hits, toy$genes)   # motivi 4..8: K = 40, E = 8 -> NA
  expect_equal(sum(!is.na(res$fdr)), 3)          # m = 3, non 8
  expect_equal(res$fdr, bh_hand(res$p))
  # con m = 8 (tutti i motivi) l'fdr dei motivi forti sarebbe piu' grande: il test lo distingue
  fdr_m8 <- p.adjust(res$p, "BH", n = length(toy$hits))
  expect_lt(res$fdr[1], fdr_m8[1])
})

# Regressione: con min_expected = 0 il comportamento e' quello precedente alla regola.
# Riferimenti calcolati con la versione SENZA min_expected (12/11): make_toy(N = 300),
# motivi 1-5 con K = 0, N, N - 1 (gene 1 mancante), N - 1 (gene 300 mancante), 1; lista 1:60.
test_that("enrich_corrected: min_expected = 0 riproduce il comportamento precedente", {
  toy <- make_toy(N = 300)
  toy$hits[[1]]$genes <- integer(0)
  toy$hits[[2]]$genes <- seq_len(300)
  toy$hits[[3]]$genes <- 2:300
  toy$hits[[4]]$genes <- 1:299
  toy$hits[[5]]$genes <- 7L
  res <- enrich_corrected(1:60, toy$hits, toy$genes, min_expected = 0)
  z_ref <- c(NA, NA, -1.9377799, 0.54261249, 1.7259265, -0.94329847, 1.2062664, -1.388163)
  p_ref <- c(NA, NA, 0.97367496, 0.29369832, 0.042180315, 0.82723587, 0.11385739, 0.91745629)
  fdr_ref <- c(NA, NA, 0.97367496, 0.58739664, 0.25308189, 0.97367496, 0.34157218, 0.97367496)
  expect_equal(res$z, z_ref, tolerance = 1e-7)
  expect_equal(res$p, p_ref, tolerance = 1e-7)
  expect_equal(res$fdr, fdr_ref, tolerance = 1e-7)
})

test_that("enrich_corrected: min_expected non valido e' un errore con messaggio chiaro", {
  toy <- make_toy(N = 300)
  # messaggio atteso: "`min_expected` deve essere un numero singolo finito >= 0"
  for (bad in list(-1, "10", NA_real_, c(5, 10), TRUE, NULL, Inf))
    expect_error(enrich_corrected(1:60, toy$hits, toy$genes, min_expected = bad),
                 "min_expected.*numero singolo finito", info = paste(deparse(bad), collapse = ""))
})

test_that("la soglia SCORE_VAR_MIN sta nel divario tra rumore e informazione vera", {
  toy <- make_toy(N = 300)
  toy$hits[[1]]$genes <- seq_len(300)            # K = N: solo rumore
  toy$hits[[2]]$genes <- 2:300                   # K = N - 1
  toy$hits[[3]]$genes <- 7L                      # K = 1
  sf <- score_fit(1:60, toy$hits, toy$genes)     # funzione interna: z e v / sum(w)
  expect_equal(SCORE_VAR_MIN, 1e-8)
  expect_lt(abs(sf$v_rel[1]), SCORE_VAR_MIN / 100)   # rumore: < 1e-10
  expect_gt(sf$v_rel[2], SCORE_VAR_MIN * 100)        # informazione vera: > 1e-6
  expect_gt(sf$v_rel[3], SCORE_VAR_MIN * 100)
  expect_gt(min(sf$v_rel[-(1:3)]), 1e-3)             # motivi normali (K = 40 su 300)
})

# Test esatto: z^2 dello score = statistica "Rao" di anova(modello nullo, modello pieno).
# Tolleranza relativa 1e-4. Include un motivo con x = 0 (anche Rao usa il modello nullo).
test_that("enrich_corrected: z^2 coincide con il test di Rao di anova()", {
  toy <- make_toy(N = 300)
  toy$hits[[1]]$genes <- 1:80                    # segnale piantato
  toy$hits[[2]]$genes <- 201:240                 # nessun gene della lista (x = 0)
  idx <- c(1:50, 101:110)
  res <- enrich_corrected(idx, toy$hits, toy$genes, min_expected = 0)
  d <- data.frame(y = as.numeric(seq_len(300) %in% idx), gc = toy$genes$gc,
                  ls = log1p(toy$genes$total_sites_noext))
  fit_null <- glm(y ~ gc + ls, binomial, d)
  for (j in seq_along(toy$hits)) {
    d$has <- as.numeric(seq_len(300) %in% toy$hits[[j]]$genes)
    fit_full <- glm(y ~ has + gc + ls, binomial, d)
    rao <- anova(fit_null, fit_full, test = "Rao")$Rao[2]
    expect_equal(res$z[j]^2, rao, tolerance = 1e-4, info = paste("motivo", j))
  }
  expect_gt(res$z[1], 0)                         # il segno si controlla a parte
  expect_lt(res$z[2], 0)
})

# TARATURA (dichiarata prima di eseguire). Fixture: N = 2000, M = 20 motivi indipendenti
# con K = 400 geni, L = 200 liste casuali di n = 100 geni: 4000 test per metodo, quota
# attesa di p < 0.05 pari a 0.05. I motivi nella stessa lista possono essere correlati:
# effetto di disegno 1 + (M - 1) * rho con rho = 0.1 (prudente) = 2.9, quindi
# SE = sqrt(0.05 * 0.95 / 4000 * 2.9) = 0.0059 e intervallo = 0.05 +- 4 SE.
#   corretto: [0.027, 0.073]     baseline (test discreto, esatto, conservativo): [0.02, 0.073]
# Se fallisce non cambio l'intervallo: indago.
test_that("taratura: su liste casuali p < 0.05 capita nella quota attesa", {
  toy <- make_toy(N = 2000, M = 20, K = 400, seed = 11)
  set.seed(100)
  pb <- pc <- numeric(0)
  for (l in 1:200) {
    idx <- sort(sample(2000, 100))
    pb <- c(pb, enrich_baseline(idx, toy$hits, 2000)$p)
    pc <- c(pc, enrich_corrected(idx, toy$hits, toy$genes)$p)
  }
  expect_length(pc, 4000)
  expect_false(anyNA(pc))
  rb <- mean(pb < 0.05); rc <- mean(pc < 0.05)
  expect_true(rc >= 0.027 && rc <= 0.073, info = paste("corretto:", rc))
  expect_true(rb >= 0.02 && rb <= 0.073, info = paste("baseline:", rb))
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

# --- Passo 5: regola del peso, direzione, errori -------------------------------
# Dichiarato prima del codice. Peso = (valore / media)^power, come nello script della
# griglia (data-raw/validation/step4_bias_calibration.R): power > 0 favorisce i valori
# alti, power < 0 i bassi. Default: biased = FALSE (uniforme); con biased = TRUE il
# default e' by = "gc" e power = 2. Il risultato e' ordinato, come nella griglia.
# Con lo stesso seed e gli stessi argomenti, random_list deve coincidere con
#   sort(sample.int(N, n, prob = (x / mean(x))^power))     (uniforme: sort(sample.int(N, n)))
# Un valore <= 0 della colonna vale come il minimo positivo, per qualunque segno di power.

test_that("random_list: coincide con sort(sample.int(N, n, prob = (x/mean(x))^power))", {
  toy <- make_toy(N = 500)
  N <- 500; g <- toy$genes; d <- g$total_sites_noext
  set.seed(11); a <- random_list(g, 50, biased = TRUE)           # default: by = "gc", power = 2
  set.seed(11); b <- sort(sample.int(N, 50, prob = (g$gc / mean(g$gc))^2))
  expect_identical(a, b)
  set.seed(12); a <- random_list(g, 50, biased = TRUE, by = "gc", power = -5)
  set.seed(12); b <- sort(sample.int(N, 50, prob = (g$gc / mean(g$gc))^-5))
  expect_identical(a, b)
  set.seed(13); a <- random_list(g, 50, biased = TRUE, by = "total_sites_noext", power = 10)
  set.seed(13); b <- sort(sample.int(N, 50, prob = (d / mean(d))^10))
  expect_identical(a, b)
  set.seed(14); a <- random_list(g, 50)                          # uniforme
  set.seed(14); b <- sort(sample.int(N, 50))
  expect_identical(a, b)
})

test_that("random_list: power vale 2 di default", {
  expect_equal(formals(random_list)$power, 2)
})

# Verso il basso. Atteso, dichiarato prima di eseguire (specchio dei test verso l'alto):
# con power = -2 il GC medio della lista scende di piu' di 0.03 rispetto all'uniforme;
# per total_sites_noext la densita' media scende.
test_that("random_list biased power < 0: GC medio piu' basso dell'uniforme", {
  toy <- make_toy(N = 2000)
  set.seed(15)
  m_u <- mean(replicate(20, mean(toy$genes$gc[random_list(toy$genes, 300)])))
  m_b <- mean(replicate(20, mean(toy$genes$gc[random_list(toy$genes, 300, biased = TRUE,
                                                          power = -2)])))
  expect_lt(m_b, m_u - 0.03)
})

test_that("random_list biased power < 0 by total_sites_noext: densita' media piu' bassa", {
  toy <- make_toy(N = 2000)
  d <- toy$genes$total_sites_noext
  set.seed(16)
  m_u <- mean(replicate(20, mean(d[random_list(toy$genes, 300)])))
  m_b <- mean(replicate(20, mean(d[random_list(toy$genes, 300, biased = TRUE,
                                              by = "total_sites_noext", power = -2)])))
  expect_lt(m_b, m_u)
})

test_that("random_list power > 0: un gene con valore 0 non ha peso 0 (n = N funziona)", {
  toy <- make_toy(N = 500)
  toy$genes$total_sites_noext[1] <- 0L           # con peso 0, sample.int(500, 500) fallirebbe
  x <- random_list(toy$genes, 500, biased = TRUE, by = "total_sites_noext", power = 10)
  expect_identical(x, 1:500)
})

test_that("random_list: n = 1 e n = N sono validi", {
  toy <- make_toy(N = 500)
  expect_length(random_list(toy$genes, 1), 1)
  expect_identical(random_list(toy$genes, 500), 1:500)
})

# Messaggi attesi: "`n` deve essere un intero tra 1 e 500"
#                  "`by` deve essere \"gc\" o \"total_sites_noext\""  (nome esatto, niente
#                  corrispondenza parziale: "total_sites" e' una colonna ma NON e' valida)
test_that("random_list: n fuori range o non valido e' un errore con messaggio chiaro", {
  toy <- make_toy(N = 500)
  for (bad in list(0, -3, 501, 2.5, NA_real_, c(5, 10)))
    expect_error(random_list(toy$genes, bad), "n. deve essere un intero tra 1 e 500",
                 info = paste(deparse(bad), collapse = ""))
})

test_that("random_list: by non valido e' un errore con messaggio chiaro", {
  toy <- make_toy(N = 500)
  for (bad in list("total_sites", "GC", "foo", NA_character_, 1))
    expect_error(random_list(toy$genes, 50, biased = TRUE, by = bad),
                 "by. deve essere.*gc.*total_sites_noext", info = paste(deparse(bad), collapse = ""))
})

# Messaggio atteso: "`power` deve essere un numero singolo finito"
test_that("random_list: power non valido e' un errore con messaggio chiaro", {
  toy <- make_toy(N = 500)
  for (bad in list(NA_real_, "2", Inf, c(1, 2)))
    expect_error(random_list(toy$genes, 50, biased = TRUE, power = bad),
                 "power. deve essere un numero singolo finito", info = paste(deparse(bad), collapse = ""))
})

# === odds_ratio_one =========================================================
# odds_ratio_one(idx, motif_genes, genes): UN solo motivo, con glm(). Stesso modello di
# enrich_corrected: in_list ~ has_site + gc + log1p(total_sites_noext). Restituisce una
# riga con beta (coefficiente di has_site), se, or = exp(beta), or_low e or_high =
# exp(beta -/+ qnorm(0.975) * se) (Wald 95%) e note. Dichiarato prima del codice:
#   note = "ok"                              caso normale
#   note = "separazione o non convergenza"   glm avverte (fitted probabilities 0 o 1) o non converge
#   note = "has_site costante"               K = 0 oppure K = N
# Negli ultimi due casi i cinque campi numerici sono NA; nessun errore e nessun avviso esce
# dalla funzione. Il riferimento e' un glm fatto a mano con la stessa formula (covariate non
# standardizzate); tolleranza relativa 1e-6.

test_that("odds_ratio_one: coincide con un glm fatto a mano (beta, se, or, intervallo Wald 95%)", {
  toy <- make_toy(N = 300)
  idx <- 1:60
  motif <- c(1:40, 101:140)                      # K = 80, 40 nella lista: arricchimento moderato
  d <- data.frame(y = as.numeric(seq_len(300) %in% idx), has = as.numeric(seq_len(300) %in% motif),
                  gc = toy$genes$gc, ls = log1p(toy$genes$total_sites_noext))
  fit <- glm(y ~ has + gc + ls, binomial, d)
  b <- unname(coef(fit)["has"]); se <- unname(summary(fit)$coefficients["has", "Std. Error"])
  res <- odds_ratio_one(idx, motif, toy$genes)
  expect_named(res, c("beta", "se", "or", "or_low", "or_high", "note"))
  expect_equal(nrow(res), 1)
  expect_equal(res$beta, b, tolerance = 1e-6)
  expect_equal(res$se, se, tolerance = 1e-6)
  expect_equal(res$or, exp(b), tolerance = 1e-6)
  expect_equal(res$or_low, exp(b - qnorm(0.975) * se), tolerance = 1e-6)
  expect_equal(res$or_high, exp(b + qnorm(0.975) * se), tolerance = 1e-6)
  expect_equal(res$note, "ok")
})

test_that("odds_ratio_one: motivo = esattamente la lista (separazione) -> NA e nota, nessun avviso", {
  toy <- make_toy(N = 300)
  res <- NULL
  expect_no_warning(res <- odds_ratio_one(1:60, 1:60, toy$genes))   # un errore farebbe fallire il test
  expect_true(all(is.na(res[, c("beta", "se", "or", "or_low", "or_high")])))
  expect_equal(res$note, "separazione o non convergenza")
})

test_that("odds_ratio_one: has_site costante (K = 0 e K = N) -> NA e nota, nessun avviso", {
  toy <- make_toy(N = 300)
  for (motif in list(integer(0), seq_len(300))) {
    res <- NULL
    expect_no_warning(res <- odds_ratio_one(1:60, motif, toy$genes))
    expect_true(all(is.na(res[, c("beta", "se", "or", "or_low", "or_high")])))
    expect_equal(res$note, "has_site costante")
  }
})

test_that("odds_ratio_one: il segno di log(or) coincide con quello di z di enrich_corrected", {
  toy <- make_toy(N = 300)
  # segnale piantato: z > 0, log(or) > 0
  toy$hits[[1]]$genes <- 1:80
  idx <- c(1:50, 101:110)
  z <- enrich_corrected(idx, toy$hits, toy$genes, min_expected = 0)$z[1]
  res <- odds_ratio_one(idx, toy$hits[[1]]$genes, toy$genes)
  expect_equal(res$note, "ok")
  expect_gt(z, 0); expect_gt(log(res$or), 0)
  expect_equal(sign(log(res$or)), sign(z))
  # impoverimento: solo 2 geni della lista su K = 42 (attesi 8.4): z < 0, log(or) < 0
  toy$hits[[1]]$genes <- c(1:2, 101:140)
  z <- enrich_corrected(1:60, toy$hits, toy$genes, min_expected = 0)$z[1]
  res <- odds_ratio_one(1:60, toy$hits[[1]]$genes, toy$genes)
  expect_equal(res$note, "ok")
  expect_lt(z, 0); expect_lt(log(res$or), 0)
  expect_equal(sign(log(res$or)), sign(z))
})

test_that("odds_ratio_one: lista non valida e' un errore chiaro", {
  toy <- make_toy(N = 300)
  expect_error(odds_ratio_one(integer(0), 1:10, toy$genes), "vuota")
  expect_error(odds_ratio_one(c(1, 1, 2), 1:10, toy$genes), "ripetuti")
})

# --- odds_ratio_one: celle vuote della tabella 2x2 has_site x in_list ------------
# Dichiarato prima del codice. Con k = geni della lista con il sito, n = geni della lista,
# K = geni con il sito, N = tutti i geni, le quattro celle sono
#     k,   n - k,   K - k,   N - K - (n - k).
# Se UNA e' zero la stima di massima verosimiglianza non esiste (separazione completa o
# parziale): i cinque campi numerici sono NA e note = "separazione o non convergenza",
# SENZA chiamare glm. has_site costante (K = 0 o K = N) resta "has_site costante" e ha la
# precedenza. Il caso 59 su 60 (nessuna cella vuota) resta "ok" con un OR finito, enorme
# ma esistente. I controlli su converged e fitted.values restano.
cells2x2 <- function(idx, motif, N) {
  n <- length(idx); K <- length(motif); k <- sum(motif %in% idx)
  c(k = k, n_k = n - k, K_k = K - k, resto = N - K - (n - k))
}

test_that("odds_ratio_one: una cella vuota della tabella 2x2 -> NA e nota, nessun avviso", {
  toy <- make_toy(N = 300)
  cases <- list(
    list(zero = "k",     idx = 1:60,  motif = 101:140),   # k = 0: nessun gene della lista ha il sito
    list(zero = "n_k",   idx = 1:60,  motif = 1:80),      # n - k = 0: lista tutta dentro il motivo
    list(zero = "K_k",   idx = 1:60,  motif = 1:40),      # K - k = 0: motivo contenuto nella lista
    list(zero = "resto", idx = 1:200, motif = 151:300))   # N - K - (n - k) = 0: tutti i geni senza sito sono in lista
  for (cs in cases) {
    expect_equal(unname(cells2x2(cs$idx, cs$motif, 300)[cs$zero]), 0, info = paste("fixture:", cs$zero))
    res <- NULL
    expect_no_warning(res <- odds_ratio_one(cs$idx, cs$motif, toy$genes))
    expect_true(all(is.na(res[, c("beta", "se", "or", "or_low", "or_high")])), info = cs$zero)
    expect_equal(res$note, "separazione o non convergenza", info = cs$zero)
  }
})

test_that("odds_ratio_one: 59 su 60 in lista (nessuna cella vuota) resta 'ok' con OR finito", {
  toy <- make_toy(N = 300)
  motif <- c(1:59, 101)                          # k = 59, n - k = 1, K - k = 1, resto = 239
  expect_true(all(cells2x2(1:60, motif, 300) > 0))
  res <- odds_ratio_one(1:60, motif, toy$genes)
  expect_equal(res$note, "ok")
  expect_true(is.finite(res$or))
  expect_gt(res$or, 1)
})

# glm non deve essere chiamato quando una cella e' vuota. Si osserva mettendo una copia di
# glm che conta le chiamate nell'ambiente globale (dove source() ha definito le funzioni) e
# togliendola alla fine. Controllo positivo: nel caso normale glm e' chiamato 1 volta, quindi
# il conteggio funziona davvero.
test_that("odds_ratio_one: con una cella vuota glm non viene chiamato (controllo positivo incluso)", {
  toy <- make_toy(N = 300)
  calls <- 0
  assign("glm", function(...) { calls <<- calls + 1; stats::glm(...) }, envir = globalenv())
  tryCatch({
    odds_ratio_one(1:60, c(1:40, 101:140), toy$genes)       # normale
    expect_equal(calls, 1)
    calls <- 0
    odds_ratio_one(1:60, 101:140, toy$genes)                # k = 0
    odds_ratio_one(1:60, 1:80, toy$genes)                   # n - k = 0
    odds_ratio_one(1:60, 1:40, toy$genes)                   # K - k = 0
    odds_ratio_one(1:200, 151:300, toy$genes)               # N - K - (n - k) = 0
    expect_equal(calls, 0)
  }, finally = rm("glm", envir = globalenv()))
})
