# Passo 4 - script usa-e-getta (fuori dal progetto): liste distorte per gc e per
# total_sites_noext, verso alto e verso basso, piu' liste uniformi di controllo.
# NON ancora eseguito. Criteri dichiarati PRIMA di eseguire (anche nel piano, Revisione 4).
#
# DATI E FUNZIONI: dati reali e funzioni vere di R/enrich.R. Si testano i 852 motivi
# del filtro di default; il BH e' su quelli. Lista di n = 300 geni.
#
# UNITA' DI REPLICA = LA LISTA, non il motivo: i motivi condividono geni (e molti sono
# varianti dello stesso TF), quindi i test sulla stessa lista sono correlati e un errore
# standard binomiale sui 852 motivi sarebbe troppo ottimista. Gli SE sono calcolati dalla
# variazione tra liste.
#
# METRICHE (per impostazione, su L liste)
#   r_b, r_c : frazione di motivi con p < 0.05, media sulle liste (baseline, corretto)
#              SE = sd tra liste / sqrt(L). Atteso 0.05 se non c'e' confondimento.
#   q_b, q_c : frazione di liste con almeno un motivo a FDR < 0.05 (BH sui 852).
#              Atteso <= ~0.05 se non c'e' nulla da trovare (guida, non test esatto:
#              con test dipendenti e BH il valore esatto non e' garantito).
#   shift_gc, shift_nx : GC medio e total_sites_noext medio della lista rispetto alla
#              popolazione, in deviazioni standard (misura della distorsione ottenuta).
#
# CRITERI (soglie fissate ora, NON si cambiano dopo aver visto i risultati)
#   baseline GONFIATO : r_b >= 0.10 (almeno il doppio del nominale) E q_b >= 0.50
#                       (meta' delle liste mostra almeno una "scoperta" falsa).
#   corretto TARATO   : 0.03 <= r_c <= 0.075 E q_c <= 0.15.
#                       r_c > 0.075 = anti-conservativo (gonfiato); r_c < 0.03 = conservativo.
#                       q_c <= 0.15 = 0.05 piu' ~4.6 SE binomiali a L = 100.
#   CONTROLLO         : la riga uniforme deve avere "corretto TARATO"; se no, il metodo o lo
#                       script sono sbagliati e le altre righe non si interpretano.
#   DEFAULT di power  : la forza |power| piu' grande della griglia per cui, in tutti e 4 i
#                       casi (gc/noext x alto/basso), baseline GONFIATO e corretto TARATO.
#                       Dove non succede, si riporta caso per caso.
#   Nota: il modello corretto assume logit lineare in gc e in log1p(noext). Un bias
#   "peso = gc^power" e' log-lineare in gc, quindi per gc c'e' una lieve mancata
#   corrispondenza: se il corretto si gonfia per gc ad alta forza, e' un risultato, non un errore.
#
# DISEGNO: 1 controllo uniforme (L = 300) + 12 impostazioni distorte (2 variabili x
#   power in {2, 5, 10, -2, -5, -10}, L = 100) = 1500 liste, ~0.6 s l'una -> ~15 minuti.
#   SEED: set.seed(4100 + k) all'inizio dell'impostazione k = 1..13 (1 = controllo).
#   PESO di campionamento: (valore / media)^power (stesso rapporto relativo di valore^power);
#   power > 0 favorisce valori alti, power < 0 bassi. random_list useranno la stessa regola.

setwd("/home/giuliop/projects/tfbalance")
source("R/enrich.R")
out_dir <- "/tmp/claude-1000/-home-giuliop-projects-tfbalance/647d1c88-2cf2-4dd9-8ed7-dbf562cc4f43/scratchpad"

g <- readRDS("data/genes.rds"); h <- readRDS("data/motif_hits.rds"); N <- nrow(g)
h <- h[select_motifs(h, N)]                       # 852 motivi
n_list <- 300; L_unif <- 300; L_bias <- 100

grid <- rbind(data.frame(var = "uniform", power = 0, stringsAsFactors = FALSE),
              expand.grid(var = c("gc", "total_sites_noext"),
                          power = c(2, 5, 10, -2, -5, -10), stringsAsFactors = FALSE))

draw <- function(var, power) {
  if (power == 0) return(sort(sample.int(N, n_list)))
  x <- g[[var]]
  sort(sample.int(N, n_list, prob = (x / mean(x))^power))
}

shift <- function(idx, v) (mean(g[[v]][idx]) - mean(g[[v]])) / sd(g[[v]])

one_list <- function(var, power) {
  idx <- draw(var, power)
  b <- enrich_baseline(idx, h, N)
  cc <- enrich_corrected(idx, h, g)               # se il modello nullo fallisce: tutti NA + avviso
  c(rb = mean(b$p < 0.05), rc = mean(cc$p < 0.05, na.rm = TRUE),
    qb = any(b$fdr < 0.05), qc = any(cc$fdr < 0.05, na.rm = TRUE),
    nullfail = all(is.na(cc$p)),
    shift_gc = shift(idx, "gc"), shift_nx = shift(idx, "total_sites_noext"))
}
proto <- c(rb = 0, rc = 0, qb = 0, qc = 0, nullfail = 0, shift_gc = 0, shift_nx = 0)

summ <- function(M) data.frame(
  L = nrow(M), shift_gc = mean(M[, "shift_gc"]), shift_nx = mean(M[, "shift_nx"]),
  r_b = mean(M[, "rb"]), se_rb = sd(M[, "rb"]) / sqrt(nrow(M)), q_b = mean(M[, "qb"]),
  r_c = mean(M[, "rc"]), se_rc = sd(M[, "rc"]) / sqrt(nrow(M)), q_c = mean(M[, "qc"]),
  nullfail = sum(M[, "nullfail"]))

results <- vector("list", nrow(grid))
for (k in seq_len(nrow(grid))) {
  set.seed(4100 + k)
  L <- if (grid$var[k] == "uniform") L_unif else L_bias
  t0 <- Sys.time()
  suppressWarnings(M <- t(vapply(seq_len(L), function(l) one_list(grid$var[k], grid$power[k]), proto)))
  results[[k]] <- cbind(grid[k, ], summ(M))
  saveRDS(results, file.path(out_dir, "step4_results.rds"))   # salvataggio dopo ogni impostazione
  cat(sprintf("[%2d/%d] %-18s power=%4g  L=%3d  %4.0f s\n", k, nrow(grid), grid$var[k], grid$power[k],
              L, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

res <- do.call(rbind, results)
res$baseline_gonfiato <- res$r_b >= 0.10 & res$q_b >= 0.50
res$corretto_tarato <- res$r_c >= 0.03 & res$r_c <= 0.075 & res$q_c <= 0.15
options(width = 200)
print(cbind(res[, 1:2], round(res[, c("L", "shift_gc", "shift_nx", "r_b", "se_rb", "q_b", "r_c", "se_rc", "q_c")], 3),
            res[, c("nullfail", "baseline_gonfiato", "corretto_tarato")]), row.names = FALSE)
cat("\nCONTROLLO uniforme (corretto TARATO):", res$corretto_tarato[res$var == "uniform"], "\n")
