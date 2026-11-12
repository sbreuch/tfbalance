# R/analyze.R
# Logica dell'analisi per l'app: dal testo scritto dall'utente alla tabella dei
# risultati. NON contiene statistica (sta tutta in R/enrich.R): qui si leggono i
# nomi, si controllano i limiti, si applica il filtro dei motivi e si mettono in
# ordine i risultati. Nessuna UI, nessuna lettura di file: `genes` e `hits` sono
# argomenti. Richiede che R/enrich.R sia gia' stato caricato (source).

# Massimo di geni riconosciuti: e' il massimo misurato nel Passo 7 (tempi).
MAX_LIST <- 5000L
# Massimo di token grezzi nel testo (ripetizioni e ignoti inclusi): protegge da
# incolla enormi prima ancora di cercare i geni. Controllo fatto PRIMA di map_ids.
MAX_TOKENS <- 50000L

# analyze_list: `text` e' un unico testo con i nomi dei geni (simboli, Entrez o
# Ensembl) separati da spazi, tab, virgole, punto e virgola o a capo.
# Restituisce una lista:
#   table           un riga per motivo tenuto dal filtro, ordinata per p_corrected
#                   crescente (NA in fondo, a parita' decide p_baseline)
#   not_found       token non riconosciuti, senza ripetizioni
#   n_list          geni riconosciuti (unici)
#   n_duplicates    token "in piu'": ripetizioni, o lo stesso gene scritto in modi
#                   diversi (anche ignoti ripetuti)
#   n_motifs, n_excluded_low, n_excluded_high   motivi tenuti / esclusi dal filtro
#   n_corrected     motivi con il test corretto calcolato (z non NA)
#   corrected_note  NULL se n_corrected > 0, altrimenti una spiegazione
analyze_list <- function(text, genes, hits, frac_range = c(0.005, 0.9), min_expected = 10) {
  # --- 1. parsing -----------------------------------------------------------
  # strsplit() spezza il testo; "[[:space:],;]+" = uno o piu' tra spazio, tab, a capo
  # (compreso \r di \r\n), virgola e punto e virgola. I pezzi vuoti si scartano.
  text <- as.character(text)
  text <- text[!is.na(text)]
  tokens <- unlist(strsplit(paste(text, collapse = "\n"), "[[:space:],;]+"))
  tokens <- tokens[nzchar(tokens)]

  if (length(tokens) > MAX_TOKENS)
    stop("Too many names in the text (", length(tokens), "): the maximum is ",
         MAX_TOKENS, " (repeats included)", call. = FALSE)

  # --- 2. nomi -> righe di `genes` ------------------------------------------
  m <- map_ids(tokens, genes)          # una sola chiamata: costa ~0.1 s anche con 5000 token
  idx <- unique(m$idx)                 # (map_ids li da' gia' unici: qui per sicurezza)
  if (length(idx) == 0)
    stop("No gene recognised: check that the names are symbols, Entrez IDs ",
         "or Ensembl IDs, separated by spaces, commas, semicolons or new lines",
         call. = FALSE)
  if (length(idx) > MAX_LIST)
    stop("The list has ", length(idx), " recognised genes: the maximum is ", MAX_LIST,
         call. = FALSE)
  N <- nrow(genes)
  check_list(idx, N)                   # errori (es. "all genes") passano com'e'

  # token in piu' = token puliti - geni unici - ignoti unici (map_ids toglie da solo
  # le ripetizioni, quindi la differenza e' proprio il numero di ripetizioni)
  n_duplicates <- max(0L, length(tokens) - length(idx) - length(m$not_found))

  # --- 3. filtro dei motivi -------------------------------------------------
  keep <- select_motifs(hits, N, frac_range)   # controlla anche frac_range
  # esclusi sotto il minimo: motivi che NON passano il filtro c(minimo, 1);
  # esclusi sopra il massimo: motivi che NON passano il filtro c(0, massimo).
  # Cosi' la regola della frazione resta solo in select_motifs.
  n_excluded_low <- sum(!select_motifs(hits, N, c(frac_range[1], 1)))
  n_excluded_high <- sum(!select_motifs(hits, N, c(0, frac_range[2])))
  if (!any(keep))
    stop("The filter excludes all motifs: no motif has a fraction of genes with a ",
         "site between ", frac_range[1], " and ", frac_range[2], "; widen the range",
         call. = FALSE)
  hits <- hits[keep]

  # --- 4. statistica (tutta in enrich.R) ------------------------------------
  b <- enrich_baseline(idx, hits, N)
  cc <- enrich_corrected(idx, hits, genes, min_expected)   # stesso ordine di b

  tab <- data.frame(id = b$id, tf = b$tf, k = b$k, K = b$K, expected = b$expected,
                    p_baseline = b$p, fdr_baseline = b$fdr,
                    z = cc$z, p_corrected = cc$p, fdr_corrected = cc$fdr,
                    stringsAsFactors = FALSE)
  # order(): prima p_corrected, poi p_baseline per gli spareggi; na.last = TRUE
  # mette gli NA in fondo (e dentro gli NA decide p_baseline)
  tab <- tab[order(tab$p_corrected, tab$p_baseline, na.last = TRUE), ]
  rownames(tab) <- NULL

  n_corrected <- sum(!is.na(tab$z))
  corrected_note <- NULL
  if (n_corrected == 0) {
    corrected_note <- if (all(tab$expected < min_expected))
      paste0("Corrected test not computed: no motif has E >= ", format(min_expected),
             " (E = genes of the list expected to have the site). This happens with ",
             "small lists; use the baseline test or a larger list.")
    else
      "Corrected test not computed: the model cannot be fitted for this list."
  }

  list(table = tab, not_found = m$not_found, n_list = length(idx),
       n_duplicates = n_duplicates, n_motifs = nrow(tab),
       n_excluded_low = n_excluded_low, n_excluded_high = n_excluded_high,
       n_corrected = n_corrected, corrected_note = corrected_note)
}
