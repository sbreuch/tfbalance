# R/enrich.R
# Tutta la statistica dell'app sta qui, in funzioni semplici.
# Le funzioni NON leggono file: ricevono `genes` (data.frame di genes.rds) e
# `hits` (lista di motif_hits.rds) come argomenti, cosi' si possono testare
# con piccoli dati finti.

# map_ids: trasforma i nomi scritti dall'utente in numeri di riga di `genes`.
# Riconosce simbolo, Entrez ed Ensembl (anche con ".versione"), maiuscole o
# minuscole. Restituisce list(idx = righe trovate, not_found = nomi non trovati).
map_ids <- function(x, genes) {
  if (!all(c("entrez", "symbol", "ensembl") %in% names(genes)))
    stop("`genes` deve avere le colonne entrez, symbol ed ensembl", call. = FALSE)

  x <- x[!is.na(x)]
  # Un Entrez numerico come 100000 diventerebbe "1e+05" con as.character():
  # format(scientific = FALSE) lo scrive per esteso.
  if (is.numeric(x)) x <- format(x, scientific = FALSE, trim = TRUE)
  x <- trimws(as.character(x))
  x <- x[x != ""]

  # "chiave" di confronto: maiuscola e, per Ensembl, senza ".versione"
  key <- toupper(sub("^(ENSG[0-9]+)\\.[0-9]+$", "\\1", x, ignore.case = TRUE))
  keep <- !duplicated(key)            # stesso gene scritto in piu' modi: una volta
  x <- x[keep]
  key <- key[keep]

  # match() da' la posizione della prima corrispondenza (NA se non c'e').
  # Provo prima il simbolo, poi Ensembl, poi Entrez.
  i <- match(key, toupper(genes$symbol))
  miss <- is.na(i)
  i[miss] <- match(key[miss], toupper(genes$ensembl))
  miss <- is.na(i)
  i[miss] <- match(key[miss], as.character(genes$entrez))

  list(idx = sort(unique(i[!is.na(i)])), not_found = x[is.na(i)])
}

# check_list (uso interno): controlla che la lista di geni abbia senso.
# Un errore chiaro qui e' meglio di un risultato strano piu' avanti.
check_list <- function(idx, n_genes) {
  if (length(idx) == 0)
    stop("la lista e' vuota: nessun gene riconosciuto", call. = FALSE)
  if (anyNA(idx) || any(idx < 1 | idx > n_genes | idx != round(idx)))
    stop("gli indici dei geni devono essere interi tra 1 e ", n_genes, call. = FALSE)
  if (anyDuplicated(idx))
    stop("la lista contiene indici ripetuti", call. = FALSE)
  if (length(idx) >= n_genes)
    stop("la lista contiene tutti i geni: il test non ha senso", call. = FALSE)
  invisible(idx)
}

# select_motifs: sceglie i motivi da testare in base a `frac`, la frazione di
# geni con almeno un sito. Fuori da [0.005, 0.9] il test e' poco informativo:
# sotto 0.5% ci sono pochi geni (poca potenza), sopra 90% il motivo non
# distingue i geni. Restituisce un vettore TRUE/FALSE, uno per motivo.
# Il BH si calcola poi solo sui motivi tenuti (quelli che passi alle funzioni).
select_motifs <- function(hits, n_genes, frac_range = c(0.005, 0.9)) {
  if (!is.numeric(frac_range) || length(frac_range) != 2 || frac_range[1] > frac_range[2])
    stop("`frac_range` deve essere c(minimo, massimo) con minimo <= massimo", call. = FALSE)
  # lengths() da' la lunghezza di ogni elemento: qui, quanti geni hanno il sito
  frac <- lengths(lapply(hits, `[[`, "genes")) / n_genes
  frac >= frac_range[1] & frac <= frac_range[2]
}

# enrich_baseline: test ipergeometrico, una coda (arricchimento), sui "geni con
# almeno un sito". Per ogni motivo, tra i n_genes geni:
#   K = geni con sito, n = geni nella lista, k = geni della lista con sito.
# p = P(X >= k) se si pescano n geni a caso senza reinserimento.
enrich_baseline <- function(idx, hits, n_genes) {
  check_list(idx, n_genes)
  n <- length(idx)
  in_list <- logical(n_genes)            # TRUE per i geni della lista
  in_list[idx] <- TRUE

  K <- lengths(lapply(hits, `[[`, "genes"))
  # vapply = come sapply ma dichiara il tipo del risultato (integer): piu' sicuro
  k <- vapply(hits, function(h) sum(in_list[h$genes]), integer(1))

  # phyper(q, m, n, k, lower.tail = FALSE) da' P(X > q): per P(X >= k) si usa q = k - 1
  p <- phyper(k - 1, K, n_genes - K, n, lower.tail = FALSE)
  data.frame(id = vapply(hits, function(h) h$id, ""),
             tf = vapply(hits, function(h) h$tf, ""),
             k = k, K = K, expected = n * K / n_genes,
             p = p, fdr = p.adjust(p, "BH"),
             stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Test corretto: test di punteggio (score) per il modello logistico
#   in_list ~ has_site + gc + log1p(total_sites_noext)
# Perche' lo score: glm() per ognuno dei ~1000 motivi sui 35.150 geni richiede
# ~90 s; lo score adatta UNA volta il modello senza has_site e poi, per ogni
# motivo, somma qualche numero sui soli geni che hanno il sito (niente matrici
# 35.150 x 1019).

# Soglia della guardia: un motivo vale NA se v / sum(w) < SCORE_VAR_MIN.
# v = informazione di has_site dopo aver tolto quella spiegata dalle covariate;
# sum(w) = informazione totale del modello nullo. E' relativa perche' v scala con
# la lista; 1e-8 sta ~8 ordini sopra il rumore di arrotondamento (~1e-16) e ben
# sotto l'informazione di un singolo gene. Vedi i test per le misure.
SCORE_VAR_MIN <- 1e-8

# score_fit (uso interno): per ogni motivo restituisce z (senza guardie), v_rel =
# v / sum(w) e K. Restituisce NULL se il modello senza motivo non e' stimabile
# (non converge o glm.fit avverte, ad es. separazione perfetta).
score_fit <- function(idx, hits, genes) {
  N <- nrow(genes)
  y <- numeric(N)                        # 1 = gene nella lista
  y[idx] <- 1
  # covariate standardizzate (media 0, sd 1): stesso test, ma calcoli piu' stabili
  Z <- cbind(1, as.numeric(scale(genes$gc)),
             as.numeric(scale(log1p(genes$total_sites_noext))))
  if (anyNA(Z))
    stop("gc e total_sites_noext devono essere numeri non costanti e senza NA", call. = FALSE)

  # modello nullo: solo intercetta, gc e densita' (senza has_site)
  warned <- FALSE
  fit0 <- withCallingHandlers(
    glm.fit(Z, y, family = binomial()),
    warning = function(w) { warned <<- TRUE; invokeRestart("muffleWarning") })
  if (warned || !fit0$converged) return(NULL)

  mu <- fit0$fitted.values               # probabilita' stimate di stare in lista
  w <- mu * (1 - mu)                     # pesi del modello logistico
  r <- y - mu                            # residui
  wZ <- w * Z                            # N x 3: l'unica matrice, non N x motivi
  Iinv <- solve(crossprod(Z, wZ))        # inversa dell'informazione delle covariate
  I_tot <- sum(w)

  one <- function(h) {
    i <- h$genes                         # solo i geni con il sito: si lavora sugli indici
    b <- colSums(wZ[i, , drop = FALSE])
    v <- sum(w[i]) - drop(b %*% Iinv %*% b)
    c(z = if (v > 0) sum(r[i]) / sqrt(v) else NA_real_, v_rel = v / I_tot)
  }
  out <- as.data.frame(t(vapply(hits, one, c(z = 0, v_rel = 0))))
  out$K <- lengths(lapply(hits, `[[`, "genes"))
  out
}

# enrich_corrected: z e p (una coda: arricchimento) per ogni motivo, e fdr = BH.
#
# min_expected (default 10): i motivi con E = n * K / N < min_expected (n = geni della
# lista, K = geni con il sito, N = tutti i geni) hanno z, p, fdr = NA ma restano nelle
# righe; E == min_expected e' incluso. La colonna `expected` e' sempre E. Il BH e'
# calcolato solo sui motivi tenuti: gli NA (E bassa, K = 0, K = N, guardia sulla
# varianza, modello nullo non stimabile) non contano in m. min_expected = 0 non
# esclude nessun motivo.
#
# Perche' 10: esperimento H2 (data-raw/validation/step4d_E10.R, output in
# step4d_output.txt). 300 liste uniformi nuove di n = 300 geni reali (seed 4300), 728
# motivi con E >= 10: q_c10 = 0.053 (soglia dichiarata <= 0.09) e frazione di p < 0.001
# = 0.00132 (soglia <= 0.002), cioe' test tarato. Sotto 10 non e' validato (step4c: con
# E < 5 la frazione di p < 0.001 era 4.04 volte l'attesa).
#
# Residuo noto (step4d, senza criterio): nelle fasce 10 <= E < 20 e 20 <= E < 50 la
# frazione di p < 0.001 e' 2.3 e 1.8 volte l'attesa (E >= 50: 1.06). Validato solo con
# n = 300 e liste uniformi. Con liste piccole pochi motivi raggiungono E >= 10: a n = 10
# nessuno dei 852 motivi del filtro di default (step4c), quindi tutto NA.
enrich_corrected <- function(idx, hits, genes, min_expected = 10) {
  if (!is.numeric(min_expected) || length(min_expected) != 1 ||
      !is.finite(min_expected) || min_expected < 0)
    stop("`min_expected` deve essere un numero singolo finito >= 0", call. = FALSE)
  N <- nrow(genes)
  check_list(idx, N)
  K <- lengths(lapply(hits, `[[`, "genes"))
  E <- length(idx) * K / N               # geni attesi con il sito (n * K / N, come in enrich_baseline)
  sf <- score_fit(idx, hits, genes)
  if (is.null(sf)) {
    warning("il modello senza motivo non e' stimabile (lista troppo piccola o ",
            "separabile dalle covariate): risultati NA", call. = FALSE)
    z <- rep(NA_real_, length(hits))
  } else {
    z <- sf$z
    # K = 0 o K = N: has_site e' costante. Poi la guardia relativa sulla varianza.
    z[sf$K == 0 | sf$K == N | !(sf$v_rel >= SCORE_VAR_MIN)] <- NA
  }
  z[E < min_expected] <- NA              # pochi geni attesi: la riga resta, il test no
  p <- pnorm(z, lower.tail = FALSE)      # P(Z > z): grande z = arricchimento
  data.frame(id = vapply(hits, function(h) h$id, ""),
             tf = vapply(hits, function(h) h$tf, ""),
             expected = E, z = z, p = p, fdr = p.adjust(p, "BH"),
             stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# random_list: lista dimostrativa di n geni (indici di riga di `genes`, ordinati).
#
# biased = FALSE: n geni a caso con la stessa probabilita'.
# biased = TRUE: ogni gene ha un peso (valore / media)^power, dove "valore" e' la
# colonna `by` ("gc" o "total_sites_noext"); si pesca senza reinserimento con
# sample.int(prob = peso). power > 0 favorisce i valori alti, power < 0 i bassi,
# power = 0 equivale all'uniforme. E' la stessa regola della griglia del Passo 4
# (data-raw/validation/step4_bias_calibration.R): a parita' di seed e argomenti
# il risultato coincide con sort(sample.int(N, n, prob = (x / mean(x))^power)).
# Un valore <= 0 nella colonna vale come il minimo valore positivo, per evitare
# pesi infiniti (power < 0) o nulli (power > 0). Nei dati reali non capita
# (minimo di total_sites_noext = 239, di gc = 0.09).
#
# Il default power = 2 e' una scelta dimostrativa, NON l'esito della regola di
# calibrazione del Passo 4, che non ha dato un default (il controllo uniforme non
# risultava tarato). Con power = 2 su gc la lista si sposta di circa +0.38 DS in
# GC: il baseline si gonfia (frazione di motivi con p < 0.05 = 0.151) e il test
# corretto resta vicino al nominale (0.057).
random_list <- function(genes, n, biased = FALSE, by = "gc", power = 2) {
  N <- nrow(genes)
  if (!is.numeric(n) || length(n) != 1 || is.na(n) || n != round(n) || n < 1 || n > N)
    stop("`n` deve essere un intero tra 1 e ", N, call. = FALSE)
  if (!is.character(by) || length(by) != 1 || !(by %in% c("gc", "total_sites_noext")))
    stop("`by` deve essere \"gc\" o \"total_sites_noext\"", call. = FALSE)
  if (!is.numeric(power) || length(power) != 1 || !is.finite(power))
    stop("`power` deve essere un numero singolo finito", call. = FALSE)

  if (!biased) return(sort(sample.int(N, n)))

  x <- genes[[by]]
  x[x <= 0] <- min(x[x > 0])             # valore <= 0 vale come il minimo positivo
  sort(sample.int(N, n, prob = (x / mean(x))^power))
}

# ---------------------------------------------------------------------------
# odds_ratio_one: odds ratio di UN motivo, con glm() vero (circa 0.1 s: per tutti i
# motivi sarebbe troppo lento, vedi lo score test sopra). Stesso modello di
# enrich_corrected: in_list ~ has_site + gc + log1p(total_sites_noext).
# Restituisce una riga: beta (coefficiente di has_site), se, or = exp(beta),
# or_low e or_high = exp(beta -/+ qnorm(0.975) * se) (Wald 95%), note.
#
# note = "ok"                             stima valida.
# note = "has_site costante"              il motivo ha K = 0 o K = N geni: niente da stimare.
# note = "separazione o non convergenza"  una cella della tabella 2x2 has_site x in_list e'
#                                         zero (stima inesistente, glm NON e' chiamato), oppure
#                                         glm non converge, o le probabilita' stimate sono
#                                         numericamente 0 o 1, o il coefficiente non e' finito.
# Negli ultimi due casi beta, se, or, or_low, or_high sono NA.
#
# Come si riconosce la separazione. (1) Prima di glm, dalle quattro celle k, n - k, K - k,
# N - K - (n - k): una cella vuota basta. Cattura anche la separazione parziale, che glm.fit
# non segnala (misurato: x = 0 dava OR = 3e-8 con se = 1029 e nota "ok"). Nessuna soglia su
# beta o se. (2) Dopo glm, con gli stessi due criteri di stats::glm.fit (letto in R 4.6.2):
# fit$converged, e per la binomiale any(mu < eps | mu > 1 - eps) con
# eps = 10 * .Machine$double.eps. fit$boundary (avviso "stopped at boundary value") NON e'
# usato. Mai dal testo degli avvisi di glm (i messaggi di R sono tradotti nella lingua del
# sistema). Gli avvisi di glm sono silenziati (muffleWarning) e non escono dalla funzione.
# "has_site costante" ha la precedenza sulle celle vuote.
# Gli errori di glm NON sono catturati: in casi estremi glm.fit puo' dare un errore
# ("no valid set of coefficients"), che qui uscirebbe.
# Indici di motif_genes fuori da 1..N sono ignorati (si usa %in%).
odds_ratio_one <- function(idx, motif_genes, genes) {
  N <- nrow(genes)
  check_list(idx, N)
  has <- seq_len(N) %in% motif_genes
  z975 <- qnorm(0.975)
  out <- function(note, beta = NA_real_, se = NA_real_)
    data.frame(beta = beta, se = se, or = exp(beta),
               or_low = exp(beta - z975 * se), or_high = exp(beta + z975 * se),
               note = note, stringsAsFactors = FALSE)

  if (sum(has) == 0 || sum(has) == N) return(out("has_site costante"))

  # Tabella 2x2 has_site x in_list: k = geni della lista con il sito, n = geni della
  # lista, K = geni con il sito. Se una delle quattro celle e' zero la stima di massima
  # verosimiglianza non esiste (separazione completa o parziale): glm darebbe un numero
  # enorme senza avvisare. Allora NA, senza chiamare glm.
  n <- length(idx); K <- sum(has); k <- sum(has[idx])
  if (k == 0 || n - k == 0 || K - k == 0 || N - K - (n - k) == 0)
    return(out("separazione o non convergenza"))

  muffle <- function(w) invokeRestart("muffleWarning")
  d <- data.frame(y = as.numeric(seq_len(N) %in% idx), has = as.numeric(has),
                  gc = genes$gc, ls = log1p(genes$total_sites_noext))
  fit <- withCallingHandlers(glm(y ~ has + gc + ls, family = binomial(), data = d),
                             warning = muffle)
  eps <- 10 * .Machine$double.eps        # lo stesso di glm.fit
  mu <- fit$fitted.values
  beta <- unname(coef(fit)["has"])
  if (!fit$converged || any(mu < eps | mu > 1 - eps) || is.na(beta))
    return(out("separazione o non convergenza"))
  se <- withCallingHandlers(summary(fit)$coefficients["has", "Std. Error"], warning = muffle)
  if (!is.finite(beta) || !is.finite(se)) return(out("separazione o non convergenza"))
  out("ok", beta = beta, se = unname(se))
}
