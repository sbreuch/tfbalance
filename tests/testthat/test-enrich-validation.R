# tests/testthat/test-enrich-validation.R
# Validazione su DATI REALI (data/genes.rds e data/motif_hits.rds). Salta se mancano.
# Test scritti PRIMA di odds_ratio_one. Le soglie sono quelle gia' dichiarate nel piano
# (soglie "score vs Wald, 20 motivi" e Revisione 2, punto 5); se un test fallisce si
# riporta quale motivo, i due z, K e x, e NON si cambia la soglia.
#
# FIXTURE (seed fissi): 20 motivi a caso tra gli 852 del filtro di default (seed 9100);
# lista CASUALE di 300 geni (seed 9101); lista ARRICCHITA = 150 geni del primo dei 20 motivi
# + 150 geni a caso, senza duplicati (seed 9102).
#
# DEFINIZIONI (mia lettura delle soglie del piano, da confermare):
#   x = geni della lista con il sito; K = geni con il sito.
#   A = motivi con x > 0 e 0 < K < N (il Wald e' inaffidabile per separazione fuori da A).
#   B = A, piu' |z - qnorm(0.95)| > 0.25 in ENTRAMBI i metodi (lontani dalla soglia p = 0.05).
#   Si richiedono almeno 15 motivi su 20 in B (lettura piu' severa di "confrontabili").
#   Wald di riferimento: summary(glm(in_list ~ has_site + gc + log1p(noext)))$z value.
#   Score: enrich_corrected(..., min_expected = 0) (identita' statistica, senza filtro su E).
# SOGLIE
#   lista casuale:   max |z_score - z_wald| <= 0.1 su A; stessa decisione a p < 0.05 su B.
#   lista arricchita: Spearman(z_score, z_wald) >= 0.99 su A; stesso segno su A; stessa
#                     decisione su B; il motivo arricchito ha p < 1e-6 in entrambi.
#   baseline = fisher.test(greater) su 20 motivi reali (lista casuale), tolleranza 1e-10.
#   odds_ratio_one: lista arricchita, motivo con z piu' alto -> note == "ok", OR finito e > 1,
#     nessun avviso che esca. SE QUI glm AVVERTE ("fitted probabilities numerically 0 or 1")
#     su un segnale forte: ci si FERMA e lo si segnala prima di restituire NA.

source(testthat::test_path("..", "..", "R", "enrich.R"))

.real <- new.env()
get_real <- function() {
  gp <- testthat::test_path("..", "..", "data", "genes.rds")
  hp <- testthat::test_path("..", "..", "data", "motif_hits.rds")
  testthat::skip_if_not(file.exists(gp) && file.exists(hp), "data/genes.rds o data/motif_hits.rds mancano")
  if (is.null(.real$g)) {                        # si caricano una volta sola (31 MB)
    g <- readRDS(gp); h <- readRDS(hp); N <- nrow(g)
    h <- h[select_motifs(h, N)]                  # 852 motivi
    set.seed(9100); mot <- h[sort(sample.int(length(h), 20))]
    set.seed(9101); idx_rand <- sort(sample.int(N, 300))
    set.seed(9102); idx_enr <- sort(unique(c(sample(mot[[1]]$genes, 150), sample.int(N, 150))))
    .real$g <- g; .real$mot <- mot; .real$idx_rand <- idx_rand; .real$idx_enr <- idx_enr
  }
  .real
}

# z di Wald di has_site da un glm vero (riferimento; gli avvisi di separazione dei motivi
# esclusi da A non interessano qui)
wald_z <- function(idx, motif_genes, genes) {
  N <- nrow(genes)
  d <- data.frame(y = as.numeric(seq_len(N) %in% idx), has = as.numeric(seq_len(N) %in% motif_genes),
                  gc = genes$gc, ls = log1p(genes$total_sites_noext))
  suppressWarnings(summary(glm(y ~ has + gc + ls, family = binomial(), data = d))$coefficients["has", "z value"])
}

compare_scores <- function(idx, mot, g) {
  N <- nrow(g); thr <- qnorm(0.95)
  zs <- enrich_corrected(idx, mot, g, min_expected = 0)$z
  zw <- vapply(mot, function(m) wald_z(idx, m$genes, g), 0)
  K <- lengths(lapply(mot, `[[`, "genes")); x <- vapply(mot, function(m) sum(m$genes %in% idx), 0L)
  A <- x > 0 & K > 0 & K < N & !is.na(zs) & !is.na(zw)
  B <- A & abs(zs - thr) > 0.25 & abs(zw - thr) > 0.25
  list(id = vapply(mot, function(m) m$id, ""), zs = zs, zw = zw, K = K, x = x, A = A, B = B, thr = thr)
}

describe_worst <- function(cmp) {                # per il messaggio se un test fallisce
  dz <- ifelse(cmp$A, abs(cmp$zs - cmp$zw), -Inf); j <- which.max(dz)
  sprintf("peggiore: %s z_score=%.3f z_wald=%.3f K=%d x=%d | |A|=%d |B|=%d", cmp$id[j], cmp$zs[j], cmp$zw[j],
          cmp$K[j], cmp$x[j], sum(cmp$A), sum(cmp$B))
}

test_that("validazione: baseline = fisher.test(greater) su 20 motivi reali", {
  r <- get_real(); N <- nrow(r$g); n <- length(r$idx_rand)
  res <- enrich_baseline(r$idx_rand, r$mot, N)
  expect_equal(nrow(res), 20)
  for (j in seq_along(r$mot)) {
    K <- length(r$mot[[j]]$genes); k <- sum(r$mot[[j]]$genes %in% r$idx_rand)
    tab <- matrix(c(k, K - k, n - k, N - K - (n - k)), nrow = 2)
    p_ref <- fisher.test(tab, alternative = "greater")$p.value
    expect_equal(res$p[j], p_ref, tolerance = 1e-10, info = r$mot[[j]]$id)
    expect_equal(res$k[j], k)
  }
})

test_that("validazione: score vs Wald su 20 motivi reali, lista casuale", {
  r <- get_real(); cmp <- compare_scores(r$idx_rand, r$mot, r$g)
  expect_gte(sum(cmp$B), 15)
  expect_lte(max(abs(cmp$zs - cmp$zw)[cmp$A]), 0.1, label = describe_worst(cmp))
  same <- (cmp$zs > cmp$thr) == (cmp$zw > cmp$thr)
  expect_true(all(same[cmp$B]), info = describe_worst(cmp))
})

test_that("validazione: score vs Wald su 20 motivi reali, lista arricchita", {
  r <- get_real(); cmp <- compare_scores(r$idx_enr, r$mot, r$g)
  expect_gte(sum(cmp$B), 15)
  expect_gte(cor(cmp$zs[cmp$A], cmp$zw[cmp$A], method = "spearman"), 0.99)
  expect_true(all(sign(cmp$zs[cmp$A]) == sign(cmp$zw[cmp$A])), info = describe_worst(cmp))
  same <- (cmp$zs > cmp$thr) == (cmp$zw > cmp$thr)
  expect_true(all(same[cmp$B]), info = describe_worst(cmp))
  expect_lt(pnorm(cmp$zs[1], lower.tail = FALSE), 1e-6)   # il motivo arricchito e' il primo
  expect_lt(pnorm(cmp$zw[1], lower.tail = FALSE), 1e-6)
})

test_that("validazione: lista arricchita, motivo con z piu' alto -> odds_ratio_one 'ok', OR finito", {
  r <- get_real()
  zs <- enrich_corrected(r$idx_enr, r$mot, r$g, min_expected = 0)$z
  top <- which.max(zs)
  expect_equal(top, 1L)                          # il motivo piantato deve essere quello col z piu' alto
  res <- NULL
  expect_no_warning(res <- odds_ratio_one(r$idx_enr, r$mot[[top]]$genes, r$g))
  expect_equal(res$note, "ok")
  expect_true(is.finite(res$or))
  expect_gt(res$or, 1)
})
