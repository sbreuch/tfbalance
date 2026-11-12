# Passo 4b - script usa-e-getta (fuori dal progetto), SEPARATO dal resto del Passo 4:
# liste di 3, 5 e 10 geni reali. Il modello nullo di enrich_corrected converge?
# Quante volte escono tutti NA?  NON ancora eseguito.
#
# DICHIARATO PRIMA DI ESEGUIRE
#  (a) Liste casuali di geni reali (righe di genes.rds), n in {3, 5, 10}, R = 500 liste
#      per n, SEED set.seed(4200 + n). Per ogni lista si registra se il modello nullo
#      (intercetta + gc + log1p(noext), senza motivo) da' un avviso di glm.fit, non converge,
#      o e' ok; "tutti NA" = enrich_corrected restituisce NA su tutti i motivi (score_fit NULL).
#      Atteso: ~0 liste "tutti NA" (3 punti devono essere tutti estremi nello spazio gc-densita'
#      per essere separabili dai 35.147 restanti). SOGLIA DI PROBLEMA: piu' dell'1% (> 5 liste su
#      500) "tutti NA" per un qualunque n.
#  (b) Sulle prime 100 liste di ogni n, enrich_corrected completo: quanti motivi NA (oltre
#      K=0 / K=N, che con filtro di default non ci sono), minimo di v_rel, e quota di
#      p < 0.05 per motivo (SECONDARIA, per liste cosi' piccole l'approssimazione normale e'
#      scarsa; nessuna soglia, si riporta e basta).
#  (c) Liste AVVERSARIE deterministiche (non casuali), per mostrare il limite: i n geni con
#      gc piu' alto, piu' basso, e con total_sites_noext piu' alto, piu' basso. Qui ci si
#      aspetta tutti NA (separazione): e' un'informazione, non un pass/fail.
#  Il classificatore "avviso / non converge" replica le righe del modello nullo di score_fit;
#  "tutti NA" invece usa la funzione vera.

setwd("/home/giuliop/projects/tfbalance")
source("R/enrich.R")
g <- readRDS("data/genes.rds"); h <- readRDS("data/motif_hits.rds"); N <- nrow(g)
h <- h[select_motifs(h, N)]                          # 852 motivi
Z <- cbind(1, as.numeric(scale(g$gc)), as.numeric(scale(log1p(g$total_sites_noext))))

null_status <- function(idx) {                       # replica del modello nullo di score_fit
  y <- numeric(N); y[idx] <- 1
  warned <- FALSE
  f <- withCallingHandlers(glm.fit(Z, y, family = binomial()),
                           warning = function(w) { warned <<- TRUE; invokeRestart("muffleWarning") })
  c(warned = warned, converged = f$converged)
}
all_na <- function(idx) is.null(score_fit(idx, h[1], g))   # funzione vera (un solo motivo basta: NULL dipende dal nullo)

cat("== (a) liste casuali di geni reali: modello nullo, R = 500 per n ==\n")
ra <- do.call(rbind, lapply(c(3, 5, 10), function(n) {
  set.seed(4200 + n)
  lists <- lapply(1:500, function(i) sort(sample.int(N, n)))
  st <- t(vapply(lists, null_status, c(warned = 0, converged = 0)))
  na_all <- vapply(lists, all_na, TRUE)
  data.frame(n = n, liste = 500, avviso = sum(st[, "warned"]), non_converge = sum(st[, "converged"] == 0),
             tutti_NA = sum(na_all), `% tutti_NA` = round(100 * mean(na_all), 2), check.names = FALSE)
}))
print(ra, row.names = FALSE)
cat("Problema (> 1% tutti NA) per qualche n:", any(ra$tutti_NA > 5), "\n")

cat("\n== (b) enrich_corrected completo sulle prime 100 liste di ogni n ==\n")
rb <- do.call(rbind, lapply(c(3, 5, 10), function(n) {
  set.seed(4200 + n)
  lists <- lapply(1:500, function(i) sort(sample.int(N, n)))[1:100]
  r <- lapply(lists, function(idx) suppressWarnings({
    sf <- score_fit(idx, h, g)
    if (is.null(sf)) return(c(ok = 0, na_motivi = NA, min_v_rel = NA, p05 = NA))
    cc <- enrich_corrected(idx, h, g)
    c(ok = 1, na_motivi = sum(is.na(cc$p)), min_v_rel = min(sf$v_rel), p05 = mean(cc$p < 0.05, na.rm = TRUE))
  }))
  M <- do.call(rbind, r)
  data.frame(n = n, liste_ok = sum(M[, "ok"]), motivi_NA_max = max(M[, "na_motivi"], na.rm = TRUE),
             min_v_rel = signif(min(M[, "min_v_rel"], na.rm = TRUE), 2),
             quota_p05_media = round(mean(M[, "p05"], na.rm = TRUE), 3))
}))
print(rb, row.names = FALSE)

cat("\n== (c) liste avversarie (deterministiche) ==\n")
adv <- do.call(rbind, lapply(c(3, 5, 10), function(n) {
  mk <- list(`gc piu' alto` = order(-g$gc)[1:n], `gc piu' basso` = order(g$gc)[1:n],
             `noext piu' alto` = order(-g$total_sites_noext)[1:n], `noext piu' basso` = order(g$total_sites_noext)[1:n])
  do.call(rbind, lapply(names(mk), function(nm) {
    st <- null_status(sort(mk[[nm]]))
    data.frame(n = n, lista = nm, avviso = as.logical(st["warned"]), converge = as.logical(st["converged"]),
               tutti_NA = all_na(sort(mk[[nm]])))
  }))
}))
print(adv, row.names = FALSE)
