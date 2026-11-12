# tests/testthat/test-compare.R
# Test scritti PRIMA del codice (Fase B, passo B3a): dati del grafico di confronto, soglie FDR,
# conteggi e oggetto ggplot. Fixture finte, nessuna dipendenza da data/*.rds.
#
# Interfaccia (R/compare.R, nuovo; bh_subset sta in R/enrich.R):
#   FDR_LEVEL                      0.05
#   compare_data(table)            -> list(points = data.frame(id, tf, x, y, x_clipped, y_clipped), cap)
#   fdr_threshold(p, fdr, level = FDR_LEVEL) -> il p piu' grande con fdr < level (NA se nessuno)
#   fdr_counts(table)              -> list(base_all, base_sub, corrected, n_all, n_sub)
#   compare_plot(table)            -> ggplot
# `table` e' la tabella di analyze_list (colonne id, tf, k, K, expected, p_baseline, fdr_baseline,
# z, p_corrected, fdr_corrected).

source(testthat::test_path("..", "..", "R", "enrich.R"))
source(testthat::test_path("..", "..", "R", "compare.R"))
library(ggplot2)

# --- fixture --------------------------------------------------------------------
# mk(): tabella con le colonne di analyze_list. I motivi con p_corr NA hanno z NA (test corretto
# non calcolato). Gli FDR si calcolano con p.adjust (che esclude gli NA dal numero di test):
# fdr_baseline su TUTTI i motivi, fdr_corrected sui soli motivi con p_corr non NA, come in
# analyze_list.
mk <- function(p_base, p_corr) {
  n <- length(p_base)
  data.frame(id = sprintf("MA%04d.1", seq_len(n)), tf = paste0("TF", seq_len(n)),
             k = 5L, K = 100L, expected = 20,
             p_baseline = p_base, fdr_baseline = p.adjust(p_base, "BH"),
             z = ifelse(is.na(p_corr), NA_real_, qnorm(p_corr, lower.tail = FALSE)),
             p_corrected = p_corr, fdr_corrected = p.adjust(p_corr, "BH"),
             stringsAsFactors = FALSE)
}

# Fixture principale, calcolata a mano. 10 motivi, il corretto e' calcolato solo sui motivi
# 1, 2, 3 e 10 (gli altri hanno atteso troppo basso: z NA).
#   p baseline  .004 .012 .02 .2 .3 .4 .5 .6 .8 .9
#   BH sui 10 (m = 10): p*10/rango = .04 .06 .0667 .5 .6 .667 .714 .75 .889 .9   -> FDR < .05: 1 motivo (A = 1)
#   BH sui soli 4 (m = 4), p .004 .012 .02 .9: .016 .024 .0267 .9                -> FDR < .05: 3 motivi (B = 3)
#   corretto (m = 4), p .001 .03 .2 .5: p*4/rango = .004 .06 .267 .5              -> FDR < .05: 1 motivo (C = 1)
# Le soglie in scala p: baseline sul sottoinsieme p* = .02; corretto p* = .001;
# (baseline su tutti i motivi avrebbe p* = .004: e' la scelta sbagliata che i test devono escludere).
bh_p <- c(.004, .012, .02, .2, .3, .4, .5, .6, .8, .9)
fx <- mk(bh_p, c(.001, .03, .2, rep(NA, 6), .5))

# === compare_data ===============================================================

test_that("compare_data: restituisce points e cap; solo motivi con z non NA, ordine e id conservati", {
  d <- compare_data(fx)
  expect_named(d, c("points", "cap"))
  pts <- d$points
  expect_true(all(c("id", "tf", "x", "y", "x_clipped", "y_clipped") %in% names(pts)))
  expect_identical(pts$id, c("MA0001.1", "MA0002.1", "MA0003.1", "MA0010.1"))
  expect_identical(pts$tf, c("TF1", "TF2", "TF3", "TF10"))
})

test_that("compare_data: x = -log10(p baseline), y = -log10(p corretto), p = 1 vale 0", {
  pts <- compare_data(fx)$points
  expect_equal(pts$x, -log10(c(.004, .012, .02, .9)))
  expect_equal(pts$y, -log10(c(.001, .03, .2, .5)))
  expect_false(any(pts$x_clipped | pts$y_clipped))
  one <- compare_data(mk(c(1, .5), c(1, .5)))$points
  expect_equal(one$x, c(0, -log10(.5)))
  expect_false(any(is.infinite(one$x)))
})

test_that("compare_data: nessun motivo con z non NA -> 0 punti, colonne presenti, tetto finito", {
  d <- compare_data(mk(rep(.5, 3), rep(NA_real_, 3)))
  expect_equal(nrow(d$points), 0)
  expect_true(all(c("id", "x", "y", "x_clipped", "y_clipped") %in% names(d$points)))
  expect_equal(d$cap, 1)
})

# === underflow ==================================================================

test_that("underflow (p = 0 o p < 2.2e-308): disegnato al tetto comune, flag per asse", {
  # finiti: x = 2 (riga 3), .301 (riga 4); y = 1.699, 3, 5 (righe 1, 2, 4): massimo 5 -> tetto 1.1*5 + 1
  d <- compare_data(mk(c(0, 1e-320, .01, .5), c(.02, .001, 0, 1e-5)))
  pts <- d$points
  expect_equal(d$cap, 1.1 * 5 + 1)
  expect_identical(pts$x_clipped, c(TRUE, TRUE, FALSE, FALSE))
  expect_identical(pts$y_clipped, c(FALSE, FALSE, TRUE, FALSE))
  expect_equal(pts$x[1:2], c(d$cap, d$cap))          # anche il denormale 1e-320 va al tetto
  expect_equal(pts$y[3], d$cap)
  expect_equal(pts$x[3:4], -log10(c(.01, .5)))        # gli altri restano al valore vero
  expect_equal(pts$y[c(1, 2, 4)], -log10(c(.02, .001, 1e-5)))
  expect_true(all(is.finite(c(pts$x, pts$y))))
  expect_true(all(c(pts$x, pts$y) <= d$cap))
})

test_that("underflow: il confine e' p < .Machine$double.xmin (uguale NON e' troncato)", {
  pts <- compare_data(mk(c(.Machine$double.xmin, .Machine$double.xmin / 2), c(.5, .5)))$points
  expect_identical(pts$x_clipped, c(FALSE, TRUE))
  expect_equal(pts$x[1], -log10(.Machine$double.xmin))
})

test_that("underflow: tutti i p a 0 -> tetto 1 (nessun valore finito), niente Inf ne' NaN", {
  d <- compare_data(mk(c(0, 0), c(0, 0)))
  expect_equal(d$cap, 1)
  expect_equal(d$points$x, c(1, 1))
  expect_equal(d$points$y, c(1, 1))
  expect_true(all(d$points$x_clipped & d$points$y_clipped))
})

# === fdr_threshold ==============================================================

test_that("FDR_LEVEL e' 0.05", {
  expect_equal(FDR_LEVEL, 0.05)
})

test_that("fdr_threshold: il p piu' grande tra i motivi con FDR < livello (input non ordinato, ripetuti)", {
  expect_equal(fdr_threshold(c(.03, .01, .01, .2), c(.04, .04, .04, .5)), .03)
  expect_equal(fdr_threshold(c(.5, .004, .02, .012), c(.9, .04, .0667, .06)), .004)
  expect_equal(fdr_threshold(c(.5, .004, .02, .012), c(.9, .04, .0667, .06), level = .1), .02)
})

test_that("fdr_threshold: nessun significativo -> NA (e FDR uguale al livello non e' < livello)", {
  expect_identical(fdr_threshold(c(.01, .02), c(.05, .06)), NA_real_)
  expect_identical(fdr_threshold(c(.01, .02), c(.5, .6)), NA_real_)
  expect_identical(fdr_threshold(numeric(0), numeric(0)), NA_real_)
})

test_that("fdr_threshold: gli NA vengono ignorati", {
  expect_equal(fdr_threshold(c(.01, NA, .02), c(.04, NA, NA)), .01)
  expect_identical(fdr_threshold(c(NA_real_, NA_real_), c(NA_real_, NA_real_)), NA_real_)
})

test_that("fdr_threshold: {fdr < 0.05} e {p <= p*} sono lo stesso insieme (BH e' monotono in p)", {
  for (s in 1:20) {
    set.seed(100 + s)
    p <- c(runif(5, 0, 1e-3), runif(25)); fdr <- p.adjust(p, "BH")
    pstar <- fdr_threshold(p, fdr); sig <- fdr < .05
    if (any(sig)) expect_identical(p <= pstar, sig) else expect_true(is.na(pstar))
  }
})

# === fdr_counts =================================================================

test_that("fdr_counts: A (baseline su tutti) = 1, B (baseline sul sottoinsieme) = 3, C (corretto) = 1", {
  cnt <- fdr_counts(fx)
  expect_equal(cnt$base_all, 1)       # m = 10
  expect_equal(cnt$base_sub, 3)       # m = 4: stessi p, famiglia diversa, risultato diverso
  expect_equal(cnt$corrected, 1)      # m = 4
  expect_equal(cnt$n_all, 10)
  expect_equal(cnt$n_sub, 4)
  expect_false(cnt$base_all == cnt$base_sub)
})

test_that("fdr_counts: B e C stanno nella stessa famiglia (stesso insieme, stesso m)", {
  # se i p del baseline e del corretto sono gli stessi sugli stessi motivi, B e C devono coincidere
  for (s in 1:10) {
    set.seed(200 + s)
    p <- c(runif(6, 0, .002), runif(34)); keep <- runif(40) < .5
    tab <- mk(p, ifelse(keep, p, NA_real_))
    cnt <- fdr_counts(tab)
    expect_equal(cnt$base_sub, cnt$corrected)
    expect_equal(cnt$n_sub, sum(!is.na(tab$fdr_corrected)))
    expect_equal(cnt$n_all, 40)
  }
})

test_that("fdr_counts: nessun significativo -> 0 (non NA)", {
  cnt <- fdr_counts(mk(rep(.5, 5), rep(.5, 5)))
  expect_identical(c(cnt$base_all, cnt$base_sub, cnt$corrected), c(0L, 0L, 0L))
})

test_that("fdr_counts: nessun motivo con test corretto -> B e C NA (non 0), A calcolato", {
  cnt <- fdr_counts(mk(rep(.01, 4), rep(NA_real_, 4)))   # BH su 4 p uguali a .01 -> FDR .01
  expect_equal(cnt$base_all, 4)
  expect_true(is.na(cnt$base_sub))
  expect_true(is.na(cnt$corrected))
  expect_equal(cnt$n_sub, 0)
  expect_equal(cnt$n_all, 4)
})

# === compare_plot ===============================================================

geoms_of <- function(p) unname(vapply(p$layers, function(l) class(l$geom)[1], ""))

test_that("compare_plot: strati punti, diagonale, linea orizzontale (corretto) e verticale (baseline)", {
  p <- compare_plot(fx)
  expect_s3_class(p, "ggplot")
  expect_identical(geoms_of(p), c("GeomPoint", "GeomAbline", "GeomHline", "GeomVline"))
  expect_equal(layer_data(p, 2)$slope, 1)
  expect_equal(layer_data(p, 2)$intercept, 0)
  expect_equal(layer_data(p, 3)$yintercept, -log10(.001))     # corretto: p* = .001
})

test_that("compare_plot: la linea del baseline usa l'FDR sul sottoinsieme, NON fdr_baseline di tutti i motivi", {
  x <- layer_data(compare_plot(fx), 4)$xintercept
  expect_equal(x, -log10(.02))                                # sottoinsieme (m = 4): p* = .02
  expect_false(isTRUE(all.equal(x, -log10(.004))))           # tutti i motivi (m = 10): p* = .004
})

test_that("compare_plot: senza significativi nessuna linea (solo punti e diagonale)", {
  p <- compare_plot(mk(rep(.5, 5), rep(.5, 5)))
  expect_identical(geoms_of(p), c("GeomPoint", "GeomAbline"))
})

test_that("compare_plot: una linea per test, solo se quel test ha significativi", {
  # baseline: BH (m = 5) .005 .005 .005 .6 .6 -> 3 significativi, p* = .003; corretto: nessuno
  p <- compare_plot(mk(c(.001, .002, .003, .5, .6), c(.5, .6, .7, .8, .9)))
  expect_identical(geoms_of(p), c("GeomPoint", "GeomAbline", "GeomVline"))
  expect_equal(layer_data(p, 3)$xintercept, -log10(.003))
})

test_that("compare_plot: se la soglia cade su p = 0 la linea sta al tetto (finita)", {
  tab <- mk(c(0, .5, .6), c(0, .5, .6))
  p <- compare_plot(tab)
  cap <- compare_data(tab)$cap
  expect_identical(geoms_of(p), c("GeomPoint", "GeomAbline", "GeomHline", "GeomVline"))
  expect_equal(layer_data(p, 3)$yintercept, cap)
  expect_equal(layer_data(p, 4)$xintercept, cap)
})

test_that("compare_plot: etichette italiane, assi con gli stessi limiti (0, tetto), rapporto 1:1", {
  p <- compare_plot(fx)
  expect_equal(get_labs(p)$x, "-log10 p baseline")
  expect_equal(get_labs(p)$y, "-log10 p corretto")
  cap <- compare_data(fx)$cap
  expect_equal(p$scales$get_scales("x")$limits, c(0, cap))
  expect_equal(p$scales$get_scales("y")$limits, c(0, cap))
  expect_equal(p$coordinates$ratio, 1)
})

test_that("compare_plot: forma diversa solo per i punti troncati; nessun Inf in nessuno strato", {
  tab <- mk(c(0, .01, .02), c(.5, .001, 0))         # punti 1 e 3 troncati, punto 2 no
  p <- compare_plot(tab)
  sh <- layer_data(p, 1)$shape
  expect_equal(length(unique(sh)), 2)
  expect_equal(sh[1], sh[3])
  expect_false(sh[2] == sh[1])
  expect_equal(length(unique(layer_data(compare_plot(fx), 1)$shape)), 1)   # nessun troncato: una forma
  for (i in seq_along(p$layers)) {
    d <- layer_data(p, i)
    cols <- intersect(c("x", "y", "xintercept", "yintercept"), names(d))
    expect_true(all(is.finite(unlist(d[cols]))), info = paste("strato", i))
  }
})
