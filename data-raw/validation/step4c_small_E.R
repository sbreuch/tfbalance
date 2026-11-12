# Passo 4c - Opzione A (script usa-e-getta, NON ancora eseguito).
# Si lancia dalla radice del progetto, un solo processo R: Rscript data-raw/validation/step4c_small_E.R
#
# IPOTESI H1: i motivi con pochi geni attesi E = n * K / N hanno la coda superiore dello
# score troppo ottimista (p troppo piccoli), e questo spiega l'eccesso di q_c = 0.183 nel
# controllo uniforme del Passo 4. (K = geni con il sito, n = geni nella lista, N = 35.150.)
#
# DATI: le stesse 300 liste uniformi della griglia (set.seed(4101), 300 chiamate consecutive
# a sort(sample.int(N, 300))), SOLO enrich_corrected sui 852 motivi del filtro di default.
# enrich_corrected non usa numeri casuali, quindi le liste sono identiche a quelle della
# griglia. VERIFICA DI IDENTITA' (prima di tutto il resto): r_c e q_c ricalcolati qui devono
# coincidere con la riga "uniform" di step4_results.rds (tolleranza 1e-10); se no, lo script
# si ferma e non si interpreta nulla.
# Si salvano p, z, K, E per motivo e per lista (255.600 righe) in step4c_uniform_scores.rds.
#
# METRICHE (fasce di E: <2, 2-5, 5-10, >=10; bordi: E = 2 -> K = 234; 5 -> 586; 10 -> 1172 a n = 300)
#   - per fascia: frazione di p < 0.001 su tutte le coppie (lista, motivo) della fascia;
#     attesa 0.001 se calibrato; si riportano motivi, test, eventi, eventi attesi, rapporto
#     = frazione / 0.001, e un SE tra liste (solo informativo: i criteri sotto non lo usano).
#   - q_c5: frazione di liste con almeno un motivo a FDR < 0.05 considerando SOLO i motivi con
#     E >= 5, con il BH RICALCOLATO su quel sottoinsieme (m = numero di motivi con E >= 5).
#   - conteggio dei 852 motivi per fascia, per n = 300, 50 e 10 (solo tabulazione di K).
#
# CRITERI (fissati ora, non si cambiano dopo aver visto i risultati)
#   (a)  fascia E < 2:   frazione di p < 0.001 >= 3 * 0.001 = 0.003
#   (a') fascia accorpata E < 5 (E<2 piu' 2<=E<5): frazione di p < 0.001 >= 0.003
#        [aggiunto prima del run, su richiesta: la fascia E<2 a n = 300 ha solo 9 motivi]
#   (b)  fascia E >= 10: frazione di p < 0.001 <= 2 * 0.001 = 0.002
#   (c)  q_c5 <= 0.09    (0.05 + 3 SE, SE = sqrt(0.05 * 0.95 / 300) = 0.0126)
#   H1 CONFERMATA se valgono (b) e (c) E ALMENO UNO tra (a) e (a').
#   (a) e (a') si riportano comunque separatamente.
#   Altrimenti H1 NON CONFERMATA: si riporta cosi' e basta, senza provare altre spiegazioni.
#   Le fasce 2-5 e 5-10 si riportano senza criterio.
#
# LIMITI NOTI PRIMA DI ESEGUIRE (da dichiarare con i risultati)
#   * A n = 300 la fascia E < 2 contiene solo 9 motivi (K tra 191 e 234, quindi E tra 1.6 e 2.0):
#     2700 test, 2.7 eventi attesi con p < 0.001 se calibrato; il criterio (i) richiede >= 9 eventi
#     (la soglia 0.003 su 2700 test = 8.1 -> servono 9). Test valido ma poco potente.
#   * Con n = 300 non esistono motivi con E < 1.6. I valori di E davvero piccoli (n = 50: 137 motivi
#     con E < 2; n = 10: 364 con E < 2 e nessuno con E >= 10) NON sono coperti da questo esperimento,
#     quindi l'esito non dice nulla sulle liste piccole.
#   * Le coppie (lista, motivo) non sono indipendenti (stessi motivi in tutte le liste, motivi
#     che condividono geni): le frazioni sono stime, non hanno un IC binomiale.

source_dir <- "data-raw/validation"
source("R/enrich.R")
g <- readRDS("data/genes.rds"); h <- readRDS("data/motif_hits.rds"); N <- nrow(g)
h <- h[select_motifs(h, N)]                                   # 852 motivi
n_list <- 300; L <- 300
K <- lengths(lapply(h, `[[`, "genes")); ids <- vapply(h, function(m) m$id, "")
band <- function(E) cut(E, c(-Inf, 2, 5, 10, Inf), right = FALSE,
                        labels = c("E<2", "2<=E<5", "5<=E<10", "E>=10"))

cat("== motivi per fascia di E (852 motivi) ==\n")
print(sapply(c(`n=300` = 300, `n=50` = 50, `n=10` = 10), function(n) table(band(n * K / N))))

set.seed(4101)
lists <- lapply(seq_len(L), function(l) sort(sample.int(N, n_list)))
t0 <- Sys.time(); rc <- qc <- numeric(L)
scores <- lapply(seq_len(L), function(l) {
  cc <- enrich_corrected(lists[[l]], h, g)
  rc[l] <<- mean(cc$p < 0.05, na.rm = TRUE); qc[l] <<- any(cc$fdr < 0.05, na.rm = TRUE)
  data.frame(lista = l, id = cc$id, K = K, E = n_list * K / N, z = cc$z, p = cc$p)
})
sc <- do.call(rbind, scores)
saveRDS(sc, file.path(source_dir, "step4c_uniform_scores.rds"))
cat(sprintf("\n%d liste x %d motivi calcolati in %.0f s\n", L, length(h), as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# verifica di identita' con la griglia
grid <- do.call(rbind, readRDS(file.path(source_dir, "step4_results.rds")))
u <- grid[grid$var == "uniform", ]
same <- isTRUE(all.equal(mean(rc), u$r_c, tolerance = 1e-10)) && isTRUE(all.equal(mean(qc), u$q_c, tolerance = 1e-10))
cat(sprintf("identita' con la griglia: r_c %.6f vs %.6f | q_c %.4f vs %.4f -> %s\n",
            mean(rc), u$r_c, mean(qc), u$q_c, if (same) "OK" else "DIVERSE"))
if (!same) stop("le liste non coincidono con quelle della griglia: nessuna interpretazione")

sc$band <- band(sc$E)
tab <- do.call(rbind, lapply(levels(sc$band), function(b) {
  s <- sc[sc$band == b, ]; ev <- s$p < 0.001
  per_list <- tapply(ev, s$lista, mean)
  data.frame(fascia = b, motivi = length(unique(s$id)), test = nrow(s), eventi = sum(ev),
             attesi = round(0.001 * nrow(s), 1), frazione = round(mean(ev), 5),
             se_liste = round(sd(per_list) / sqrt(L), 5), rapporto = round(mean(ev) / 0.001, 2))
}))
cat("\n== frazione di p < 0.001 per fascia di E (attesa 0.001) ==\n"); print(tab, row.names = FALSE)

s5 <- sc[sc$E >= 5, ]
q5 <- vapply(split(s5, s5$lista), function(s) any(p.adjust(s$p, "BH") < 0.05, na.rm = TRUE), TRUE)
q_c5 <- mean(q5)
cat(sprintf("\nq_c5 (solo E >= 5, BH ricalcolato su %d motivi): %.3f (SE %.3f) | q_c originale (852 motivi): %.3f\n",
            length(unique(s5$id)), q_c5, sqrt(q_c5 * (1 - q_c5) / L), mean(qc)))

# fascia accorpata E < 5 (criterio (a'))
ev5 <- sc$p[sc$E < 5] < 0.001
cat(sprintf("\nfascia accorpata E<5: motivi %d, test %d, eventi %d (attesi %.1f), frazione %.5f, rapporto %.2f\n",
            length(unique(sc$id[sc$E < 5])), length(ev5), sum(ev5), 0.001 * length(ev5), mean(ev5), mean(ev5) / 0.001))

cr <- c(a = tab$frazione[tab$fascia == "E<2"] >= 0.003, a1 = mean(ev5) >= 0.003,
        b = tab$frazione[tab$fascia == "E>=10"] <= 0.002, c = q_c5 <= 0.09)
cat("\nCRITERI: (a) E<2 frazione >= 0.003:", cr["a"], "| (a') E<5 frazione >= 0.003:", cr["a1"],
    "| (b) E>=10 frazione <= 0.002:", cr["b"], "| (c) q_c5 <= 0.09:", cr["c"], "\n")
cat("H1 = (b) e (c) e almeno uno tra (a) e (a'):", cr["b"] && cr["c"] && (cr["a"] || cr["a1"]), "\n")
cat(if (cr["b"] && cr["c"] && (cr["a"] || cr["a1"])) "H1 CONFERMATA\n" else "H1 NON CONFERMATA\n")
