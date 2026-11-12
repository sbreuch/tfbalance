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

# enrich_corrected: z e p (una coda: arricchimento) per ogni motivo, e fdr = BH
# calcolato sui motivi passati; i motivi non testabili (NA) non contano in m.
enrich_corrected <- function(idx, hits, genes) {
  N <- nrow(genes)
  check_list(idx, N)
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
  p <- pnorm(z, lower.tail = FALSE)      # P(Z > z): grande z = arricchimento
  data.frame(id = vapply(hits, function(h) h$id, ""),
             tf = vapply(hits, function(h) h$tf, ""),
             z = z, p = p, fdr = p.adjust(p, "BH"),
             stringsAsFactors = FALSE)
}
