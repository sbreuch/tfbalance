# Passo 7 - tempi (script usa-e-getta). Dalla radice del progetto, un solo processo R:
#   Rscript data-raw/validation/step7_timing.R | tee data-raw/validation/step7_timing_output.txt
#
# Con le funzioni vere di R/enrich.R misura enrich_baseline + enrich_corrected (default:
# min_expected = 10) su liste di geni reali (righe di genes.rds) di 10, 50, 300 e 5000 geni
# (seed 7000 + n, stesse liste a ogni esecuzione), con i motivi del filtro di default (852) e
# con il filtro disattivato (select_motifs(..., c(0, 1)) = 1019). 3 ripetizioni per caso:
# si riportano la mediana e il minimo-massimo del tempo totale (secondi, tempo reale).
# Per ogni caso: motivi con z non NA nel corretto e avvisi emessi.

source("R/enrich.R")
g <- readRDS("data/genes.rds"); h <- readRDS("data/motif_hits.rds"); N <- nrow(g)
cat("R", as.character(getRversion()), "| BLAS:", extSoftVersion()[["BLAS"]], "| OPENBLAS_NUM_THREADS =",
    Sys.getenv("OPENBLAS_NUM_THREADS", "(non impostata)"), "\n\n")

sizes <- c(10, 50, 300, 5000)
lists <- lapply(sizes, function(n) { set.seed(7000 + n); sort(sample.int(N, n)) }); names(lists) <- sizes
filters <- list(`default (852)` = select_motifs(h, N), `disattivato (1019)` = select_motifs(h, N, c(0, 1)))

reps <- 3
rows <- list()
for (fn in names(filters)) for (n in sizes) {
  hits <- h[filters[[fn]]]; idx <- lists[[as.character(n)]]
  tb <- tc <- numeric(reps); nw <- 0L
  for (r in seq_len(reps)) {
    tb[r] <- system.time(b <- enrich_baseline(idx, hits, N))[["elapsed"]]
    tc[r] <- system.time(cc <- withCallingHandlers(enrich_corrected(idx, hits, g),
               warning = function(w) { nw <<- nw + 1L; invokeRestart("muffleWarning") }))[["elapsed"]]
  }
  tot <- tb + tc
  rows[[length(rows) + 1]] <- data.frame(filtro = fn, n = n, motivi = length(hits),
    baseline_s = round(median(tb), 3), corretto_s = round(median(tc), 3), totale_s = round(median(tot), 3),
    totale_min_max = sprintf("%.2f-%.2f", min(tot), max(tot)),
    z_non_NA = sum(!is.na(cc$z)), avvisi = nw)
}
res <- do.call(rbind, rows)
options(width = 200)
print(res, row.names = FALSE)
cat("\nz_non_NA = motivi con z non NA nel corretto con min_expected = 10 (default); avvisi = totale su", reps, "ripetizioni\n")
