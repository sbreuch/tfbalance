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
