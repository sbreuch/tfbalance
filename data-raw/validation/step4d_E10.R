# Passo 4d - esperimento H2 (script usa-e-getta, NON ancora eseguito).
# Si lancia dalla radice del progetto, un solo processo R:
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript data-raw/validation/step4d_E10.R | tee data-raw/validation/step4d_output.txt
#
# IPOTESI H2: sui motivi con E >= 10 (E = n * K / N geni attesi con il sito) il test corretto
# (score, enrich_corrected) e' tarato. A n = 300, E >= 10 equivale a K >= 1172.
#
# DATI: 300 liste uniformi NUOVE di n = 300 geni: set.seed(4300), poi 300 chiamate consecutive
# a sort(sample.int(N, 300)). Non sono le liste della griglia ne' del 4c (seed 4101).
# Solo enrich_corrected, sui 852 motivi del filtro di default (enrich_corrected non usa numeri
# casuali). Si salvano p, z, K, E per motivo e per lista in step4d_uniform_scores.rds.
#
# CRITERI (fissati ora, non si cambiano dopo aver visto i risultati)
#   (i)  q_c10 <= 0.09: frazione di liste con almeno un motivo a FDR < 0.05, considerando SOLO
#        i motivi con E >= 10 e con il BH RICALCOLATO su quel sottoinsieme.
#        (0.09 = 0.05 + 3 SE, SE = sqrt(0.05 * 0.95 / 300) = 0.0126)
#   (ii) frazione di p < 0.001 sulle coppie (lista, motivo) con E >= 10 <= 0.002 (2 volte l'attesa 0.001).
#   H2 CONFERMATA se valgono (i) E (ii). Se (i) o (ii) falliscono: H2 NON CONFERMATA; si riporta
#   cosi' e basta, senza provare altre spiegazioni.
#
# SI RIPORTA SENZA CRITERIO
#   * fasce E in [10, 20), [20, 50), [50, inf): motivi, test, eventi con p < 0.001, eventi attesi
#     (0.001 * test), frazione, SE tra liste (solo informativo), rapporto con l'attesa;
#   * q_c su tutti i 852 motivi (FDR originale di enrich_corrected, BH su 852), con SE binomiale.
#
# NOTE DICHIARATE PRIMA DI ESEGUIRE
#   * E >= 10 a n = 300 e' il 85% circa dei motivi (728 su 852 secondo la tabulazione del 4c):
#     218.400 test, 218.4 eventi attesi con p < 0.001; (ii) richiede <= 437 eventi.
#   * Sulle liste del 4c la frazione di p < 0.001 in E >= 10 era 0.00119; qui e' una replica su
#     liste nuove. q per E >= 10 non e' stato calcolato prima su nessun dato.
#   * Le coppie (lista, motivo) non sono indipendenti (stessi motivi in tutte le liste, motivi
#     che condividono geni): le frazioni sono stime senza IC binomiale; q ha SE binomiale sulle liste.

source_dir <- "data-raw/validation"
source("R/enrich.R")
g <- readRDS("data/genes.rds"); h <- readRDS("data/motif_hits.rds"); N <- nrow(g)
h <- h[select_motifs(h, N)]                                   # 852 motivi
n_list <- 300; L <- 300
K <- lengths(lapply(h, `[[`, "genes")); E <- n_list * K / N

set.seed(4300)
lists <- lapply(seq_len(L), function(l) sort(sample.int(N, n_list)))
t0 <- Sys.time(); q_all <- numeric(L)
scores <- lapply(seq_len(L), function(l) {
  cc <- enrich_corrected(lists[[l]], h, g)
  q_all[l] <<- any(cc$fdr < 0.05, na.rm = TRUE)               # q su tutti i 852 (solo riportato)
  data.frame(lista = l, id = cc$id, K = K, E = E, z = cc$z, p = cc$p)
})
sc <- do.call(rbind, scores)
saveRDS(sc, file.path(source_dir, "step4d_uniform_scores.rds"))
cat(sprintf("%d liste x %d motivi calcolati in %.0f s | p mancanti (NA): %d\n", L, length(h),
            as.numeric(difftime(Sys.time(), t0, units = "secs")), sum(is.na(sc$p))))

s10 <- sc[sc$E >= 10, ]
cat(sprintf("motivi con E >= 10: %d (test: %d)\n", length(unique(s10$id)), nrow(s10)))

# criterio (i): BH ricalcolato sui soli motivi E >= 10
q10 <- vapply(split(s10, s10$lista), function(s) any(p.adjust(s$p, "BH") < 0.05, na.rm = TRUE), TRUE)
q_c10 <- mean(q10)
cat(sprintf("\n(i)  q_c10 (BH su %d motivi): %.3f (SE %.3f)  soglia <= 0.09\n",
            length(unique(s10$id)), q_c10, sqrt(q_c10 * (1 - q_c10) / L)))

# criterio (ii): frazione di p < 0.001 su E >= 10
ev <- s10$p < 0.001
frac10 <- mean(ev, na.rm = TRUE)
cat(sprintf("(ii) frazione di p < 0.001 su E >= 10: %.5f (eventi %d, attesi %.1f, rapporto %.2f)  soglia <= 0.002\n",
            frac10, sum(ev, na.rm = TRUE), 0.001 * nrow(s10), frac10 / 0.001))

# senza criterio: fasce di E
s10$band <- cut(s10$E, c(10, 20, 50, Inf), right = FALSE, labels = c("10<=E<20", "20<=E<50", "E>=50"))
tab <- do.call(rbind, lapply(levels(s10$band), function(b) {
  s <- s10[s10$band == b, ]; e <- s$p < 0.001
  data.frame(fascia = b, motivi = length(unique(s$id)), test = nrow(s), eventi = sum(e, na.rm = TRUE),
             attesi = round(0.001 * nrow(s), 1), frazione = round(mean(e, na.rm = TRUE), 5),
             se_liste = round(sd(tapply(e, s$lista, mean)) / sqrt(L), 5),
             rapporto = round(mean(e, na.rm = TRUE) / 0.001, 2))
}))
cat("\n== (senza criterio) frazione di p < 0.001 per fascia di E >= 10 ==\n"); print(tab, row.names = FALSE)

q_c_all <- mean(q_all)
cat(sprintf("\n(senza criterio) q_c su tutti i %d motivi (BH su %d): %.3f (SE %.3f)\n", length(h), length(h),
            q_c_all, sqrt(q_c_all * (1 - q_c_all) / L)))

ok <- c(i = q_c10 <= 0.09, ii = frac10 <= 0.002)
cat("\nCRITERI: (i) q_c10 <= 0.09:", ok["i"], "| (ii) frazione p<0.001 (E>=10) <= 0.002:", ok["ii"], "\n")
cat(if (all(ok)) "H2 CONFERMATA\n" else "H2 NON CONFERMATA (si riporta e ci si ferma)\n")
