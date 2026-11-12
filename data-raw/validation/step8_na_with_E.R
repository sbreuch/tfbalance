# Passo 8 - motivi con z NA nonostante E >= min_expected (script usa-e-getta, SOLA LETTURA).
# Dalla radice del progetto:
#   Rscript data-raw/validation/step8_na_with_E.R | tee data-raw/validation/step8_na_with_E_output.txt
#
# Domanda: la didascalia dell'app dice "NA nel test corretto = atteso < N geni". E' vera?
# Cioe': esistono motivi con z NA anche se E = n * K / N >= min_expected (10)? E se si', perche'?
#
# Con le funzioni vere di R/enrich.R (enrich_corrected per i NA veri; score_fit, interna,
# solo per capire la CAUSA). Liste di geni reali (righe di genes.rds, uniformi), n = 12, 20,
# 50, 300, 5000, 20 liste per n, seed = 800000 + 100 * n + k (k = 1..20). Due filtri: default
# (select_motifs, 852 motivi) e c(0, 1) (1019 motivi, filtro disattivato).
# In piu': la lista n = 5000 del Passo 7 (seed 7000 + n), filtro disattivato, per chiudere i NA
# di quel passo.
#
# Cause (nell'ordine in cui le assegna enrich_corrected):
#   modello    : modello senza motivo non stimabile (score_fit = NULL, con avviso) -> tutti NA
#   K0_KN      : K = 0 oppure K = N (has_site costante)
#   guardia    : v_rel < SCORE_VAR_MIN (o v_rel NA / <= 0): il motivo e' quasi spiegato dalle covariate
#   altro      : z NA senza nessuna delle cause sopra (non dovrebbe esistere)

source("R/enrich.R")
g <- readRDS("data/genes.rds"); h <- readRDS("data/motif_hits.rds"); N <- nrow(g)
ME <- formals(enrich_corrected)$min_expected          # 10: letto dalla firma, non scritto a mano
K_all <- unname(lengths(lapply(h, `[[`, "genes")))   # unname: lengths() porta i nomi dei motivi
cat("N =", N, "| motivi =", length(h), "| min_expected =", ME, "| SCORE_VAR_MIN =", SCORE_VAR_MIN, "\n")
cat("motivi con K = 0:", sum(K_all == 0), "| K = N:", sum(K_all == N), "\n\n")

filters <- list(default = select_motifs(h, N), `c(0,1)` = select_motifs(h, N, c(0, 1)))
sizes <- c(12, 20, 50, 300, 5000)

# Una lista: restituisce i conteggi e, per i NA con E >= ME, id/K/E/v_rel/causa.
one_list <- function(idx, keep) {
  hits <- h[keep]; K <- K_all[keep]; E <- length(idx) * K / N
  nw <- 0L
  cc <- withCallingHandlers(enrich_corrected(idx, hits, g),
          warning = function(w) { nw <<- nw + 1L; invokeRestart("muffleWarning") })
  sf <- score_fit(idx, hits, g)                       # NULL se il modello non e' stimabile
  na <- is.na(cc$z)
  cause <- rep(NA_character_, length(hits)); v_rel <- rep(NA_real_, length(hits))
  if (is.null(sf)) {
    cause[] <- "modello"
  } else {
    v_rel <- sf$v_rel
    bad_K <- sf$K == 0 | sf$K == N
    guard <- !bad_K & !(sf$v_rel >= SCORE_VAR_MIN)
    cause[bad_K] <- "K0_KN"; cause[guard] <- "guardia"
    # z NA per la sola regola E < min_expected: nessuna causa "statistica"
  }
  low <- E < ME
  # controllo di coerenza: i NA veri = (E bassa) oppure (causa statistica)
  stat_na <- !is.na(cause)
  stopifnot(identical(unname(na), unname(low | stat_na)), (is.null(sf)) == (nw > 0))
  cause[na & !low & is.na(cause)] <- "altro"
  list(n_motifs = length(hits), n_E_ok = sum(!low), n_na_E_ok = sum(na & !low),
       n_na_E_low = sum(na & low), nw = nw,
       detail = if (any(na & !low)) data.frame(id = cc$id[na & !low], K = K[na & !low],
                  E = E[na & !low], v_rel = v_rel[na & !low], causa = cause[na & !low],
                  stringsAsFactors = FALSE) else NULL)
}

rows <- list(); details <- list()
run_case <- function(fn, n, label, idx_list) {
  res <- lapply(idx_list, one_list, keep = filters[[fn]])
  d <- do.call(rbind, lapply(seq_along(res), function(k) if (!is.null(res[[k]]$detail)) cbind(lista = k, res[[k]]$detail)))
  ca <- if (is.null(d)) character(0) else d$causa
  rows[[length(rows) + 1]] <<- data.frame(
    filtro = fn, caso = label, liste = length(res),
    coppie_E_ok = sum(vapply(res, `[[`, 0L, "n_E_ok")),       # (lista, motivo) con E >= min_expected
    NA_con_E_ok = sum(vapply(res, `[[`, 0L, "n_na_E_ok")),
    liste_con_NA = sum(vapply(res, function(r) r$n_na_E_ok > 0, TRUE)),
    modello = sum(ca == "modello"), K0_KN = sum(ca == "K0_KN"), guardia = sum(ca == "guardia"),
    altro = sum(ca == "altro"), avvisi = sum(vapply(res, `[[`, 0L, "nw")))
  if (!is.null(d)) details[[length(details) + 1]] <<- cbind(filtro = fn, caso = label, d)
}

t0 <- Sys.time()
for (fn in names(filters)) for (n in sizes) {
  idx_list <- lapply(1:20, function(k) { set.seed(800000 + 100 * n + k); sort(sample.int(N, n)) })
  run_case(fn, n, paste0("n=", n), idx_list)
}
set.seed(7000 + 5000)                                   # la lista del Passo 7 (n = 5000)
run_case("c(0,1)", 5000, "Passo7 n=5000", list(sort(sample.int(N, 5000))))

options(width = 220)
cat("== riepilogo (conteggi su tutte le liste del caso) ==\n")
print(do.call(rbind, rows), row.names = FALSE)
cat("\ncoppie_E_ok = coppie (lista, motivo) con E >= min_expected; NA_con_E_ok = di queste, quelle con z NA;\n",
    "liste_con_NA = liste (su 20) con almeno un NA con E >= min_expected; ultime colonne = causa.\n", sep = "")

if (length(details)) {
  dd <- do.call(rbind, details)
  cat("\n== NA con E >= min_expected: dettaglio per (filtro, caso, motivo), motivi ripetuti raggruppati ==\n")
  key <- paste(dd$filtro, dd$caso, dd$id, dd$K, dd$causa, sep = " | ")
  out <- data.frame(chiave = names(table(key)), volte = as.integer(table(key)),
                    E_min = as.numeric(tapply(dd$E, key, min)[names(table(key))]),
                    E_max = as.numeric(tapply(dd$E, key, max)[names(table(key))]),
                    v_rel_min = as.numeric(tapply(dd$v_rel, key, min)[names(table(key))]),
                    v_rel_max = as.numeric(tapply(dd$v_rel, key, max)[names(table(key))]))
  print(head(out[order(-out$volte), ], 40), row.names = FALSE)
  cat("\nmotivi diversi con NA e E >= min_expected:", length(unique(paste(dd$filtro, dd$id))), "su", length(h), "\n")
} else cat("\nNessun NA con E >= min_expected in nessun caso.\n")
cat("\ntempo:", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "minuti\n")
