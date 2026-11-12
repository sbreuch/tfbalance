# prepare_data.R -- da eseguire UNA volta sola, fuori dall'app.
# Produce data/genes.rds e data/motif_hits.rds (vedi i passi sotto).
# Genome build: hg38 = GRCh38 (TxDb.Hsapiens.UCSC.hg38.knownGene).

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(GenomeInfoDb)
})

# ---- Passo 1: geni umani sui cromosomi standard ---------------------------
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# genes() ritorna un intervallo per gene; il nome e' l'Entrez ID.
# Scarta (con un messaggio) i geni con esoni su piu' filamenti/cromosomi.
g <- genes(txdb)

# Solo chr1-22, X, Y (esclude contig alternativi e chrM).
standard_chr <- paste0("chr", c(1:22, "X", "Y"))
g <- keepSeqlevels(g, standard_chr, pruning.mode = "coarse")

stopifnot(identical(unique(genome(g)[standard_chr]), "hg38"))
cat("Geni su cromosomi standard:", length(g), "\n")

# ---- Passo 2: un promotore per gene (1000 bp a monte, 100 a valle del TSS) --
# promoters() rispetta il filamento: il TSS e' start() per "+" e end() per "-".
# L'intervallo e' start-1000 ... start+99 (1100 bp; la base del TSS e' la
# prima delle 100 a valle). Il TSS e' quello piu' a monte tra i trascritti.
p <- promoters(g, upstream = 1000, downstream = 100)

# Rimuovo (non accorcio con trim()) i promotori che escono dal cromosoma.
chr_len <- seqlengths(p)[as.character(seqnames(p))]
inside <- start(p) >= 1 & end(p) <= chr_len
cat("Promotori fuori dai bordi rimossi:", sum(!inside), "\n")
p <- p[inside]

stopifnot(all(width(p) == 1100), !anyNA(names(p)))
cat("Promotori tenuti:", length(p), "\n")

# ---- Passo 3: sequenze dei promotori e contenuto in GC ---------------------
suppressPackageStartupMessages({
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(Biostrings)
})

# getSeq() con un GRanges rispetta il filamento: le sequenze sono gia'
# orientate come il gene (5' -> 3').
seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, p)

# Frazione di G+C su tutte le 1100 basi (as.prob = TRUE divide per la lunghezza).
gc <- as.numeric(letterFrequency(seqs, "GC", as.prob = TRUE))

stopifnot(length(seqs) == length(p), all(width(seqs) == 1100),
          !anyNA(gc), all(gc >= 0 & gc <= 1))
# TP53 (Entrez 7157): deve stare su chr17, filamento "-".
stopifnot(as.character(seqnames(p["7157"])) == "chr17",
          as.character(strand(p["7157"])) == "-")
cat("Sequenze:", length(seqs), "| GC medio:", round(mean(gc), 3), "\n")

# Escludo i promotori con almeno una N (basi non sequenziate): falserebbero
# GC e punteggi PWM. Ne vengono tolti 11 (35161 -> 35150).
has_n <- alphabetFrequency(seqs)[, "N"] > 0
cat("Promotori con N tolti:", sum(has_n), "(Entrez:", names(p)[has_n], ")\n")
p <- p[!has_n]; seqs <- seqs[!has_n]; gc <- gc[!has_n]
stopifnot(sum(has_n) == 11, length(p) == 35150, length(seqs) == 35150,
          length(gc) == 35150, "7157" %in% names(p), !anyNA(gc),
          all(width(seqs) == 1100), all(alphabetFrequency(seqs)[, "N"] == 0))

# ---- Passo 4: simboli ed Ensembl ID da org.Hs.eg.db ------------------------
suppressPackageStartupMessages(library(org.Hs.eg.db))
# mapIds(): un valore per Entrez ID; se ce n'e' piu' d'uno prende il primo.
# Gli Ensembl ID mancanti restano NA (non scartiamo geni per questo).
symbol  <- mapIds(org.Hs.eg.db, names(p), "SYMBOL",  "ENTREZID", multiVals = "first")
ensembl <- mapIds(org.Hs.eg.db, names(p), "ENSEMBL", "ENTREZID", multiVals = "first")
stopifnot(identical(names(symbol), names(p)), identical(names(ensembl), names(p)))
cat("Senza simbolo:", sum(is.na(symbol)), "| senza Ensembl:", sum(is.na(ensembl)),
    "| simboli duplicati:", sum(duplicated(symbol[!is.na(symbol)])), "\n")

# ---- Passo 5a: lettura delle matrici JASPAR (formato pfm) -------------------
# Ogni motivo occupa 5 righe: ">ID<TAB>nome" e poi le righe A, C, G, T con i
# conteggi per posizione tra [ ]. Solo base R: nessun pacchetto serve qui.
jaspar_file <- "data-raw/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
lines  <- readLines(jaspar_file)
head_i <- grep("^>", lines)
stopifnot(length(head_i) * 5 == length(lines))
header   <- strsplit(sub("^>", "", lines[head_i]), "\t")
motif_id <- vapply(header, `[`, "", 1)
motif_tf <- vapply(header, `[`, "", 2)
stopifnot(!anyDuplicated(motif_id), !anyNA(motif_tf))

read_pfm <- function(i) {
  rows <- lines[head_i[i] + 1:4]
  stopifnot(substr(rows, 1, 1) == c("A", "C", "G", "T"))
  nums <- strsplit(trimws(gsub("^[ACGT]|[][]", "", rows)), " +")
  stopifnot(length(unique(lengths(nums))) == 1)       # 4 righe, stesse colonne
  m <- matrix(as.numeric(unlist(nums)), nrow = 4, byrow = TRUE,
              dimnames = list(c("A", "C", "G", "T"), NULL))
  stopifnot(!anyNA(m), all(m >= 0), all(m == round(m)))
  m
}
pfms <- setNames(lapply(seq_along(head_i), read_pfm), motif_id)
stopifnot(identical(pfms[["MA0069.1"]]["A", 1:3], c(2, 2, 4)))   # PAX6, dal file
w <- vapply(pfms, ncol, 1L)
cat("Motivi letti:", length(pfms), "| larghezza min/mediana/max:",
    min(w), median(w), max(w), "\n")

# ---- Passo 5b: colonne con la stessa somma (riscalamento) -------------------
# PWM() vuole una matrice di interi con tutte le colonne della stessa somma, ma
# in 528 motivi JASPAR non e' cosi'. rel_diff = (max-min)/max delle somme di
# colonna PRIMA del riscalamento: 0 vuol dire che il motivo e' gia' a posto.
rel_diff <- vapply(pfms, function(m) { s <- colSums(m); (max(s) - min(s)) / max(s) }, 0)
rescaled <- rel_diff > 0
stopifnot(all(vapply(pfms, function(m) all(colSums(m) > 0), TRUE)))

# Ogni colonna viene portata alla somma media del motivo (arrotondata): si
# moltiplica per target/somma e si arrotonda; il resto dell'arrotondamento va
# sul valore piu' grande della colonna, cosi' la somma e' esattamente target.
rescale_pfm <- function(m) {
  s <- colSums(m); target <- round(mean(s))
  out <- round(sweep(m, 2, target / s, `*`))
  top <- cbind(max.col(t(out), ties.method = "first"), seq_len(ncol(out)))
  out[top] <- out[top] + (target - colSums(out))
  out
}
pfms_eq <- pfms
pfms_eq[rescaled] <- lapply(pfms[rescaled], rescale_pfm)
pfms_eq <- lapply(pfms_eq, function(m) { storage.mode(m) <- "integer"; m })

stopifnot(all(vapply(pfms_eq, function(m) all(m >= 0) && length(unique(colSums(m))) == 1, TRUE)),
          identical(pfms_eq[!rescaled], lapply(pfms[!rescaled], function(m) { storage.mode(m) <- "integer"; m })))
# Controllo che il riscalamento non sposti le proporzioni per colonna:
shift <- max(mapply(function(a, b) max(abs(sweep(a, 2, colSums(a), "/") -
                                           sweep(b, 2, colSums(b), "/"))),
                    pfms_eq[rescaled], pfms[rescaled]))
cat("Motivi riscalati:", sum(rescaled), "su", length(pfms),
    "| spostamento max di una proporzione:", signif(shift, 2), "\n")
band <- cut(rel_diff, c(-Inf, 0, 0.01, 0.1, Inf), labels = c("0", "0-1%", "1-10%", ">10%"))
cat("rel_diff per fascia (estremo destro incluso):\n"); print(table(band))

# ---- Passo 6a: PWM e verifica della soglia min.score = "90%" ---------------
pwms <- lapply(pfms_eq, PWM)     # log2probratio, riscalata tra 0 e 1 (unitScale)

# Pezzo di sequenza reale per i test: i primi 2000 promotori messi di seguito.
test_subj <- unlist(seqs[1:2000])

# matchPWM calcola la soglia come maxScore * 90/100, non come
# minScore + 0.9*(maxScore-minScore): coincidono solo se minScore = 0.
check_threshold <- function(id) {
  pwm <- pwms[[id]]; mn <- minScore(pwm); mx <- maxScore(pwm)
  thr_matchpwm <- mx * 90 / 100
  thr_formula  <- mn + 0.9 * (mx - mn)
  n_warn <- 0
  hits <- withCallingHandlers(
    matchPWM(pwm, test_subj, min.score = "90%", with.score = TRUE),
    warning = function(w) { n_warn <<- n_warn + 1; invokeRestart("muffleWarning") })
  all_sc <- PWMscoreStartingAt(pwm, test_subj,
                               starting.at = seq_len(length(test_subj) - ncol(pwm) + 1))
  ref <- which(all_sc >= thr_formula)         # tutte le posizioni sopra soglia
  stopifnot(abs(thr_matchpwm - thr_formula) < 1e-9, n_warn == 0,
            identical(start(hits), ref),
            isTRUE(all.equal(mcols(hits)$score, all_sc[ref])))
  data.frame(id = id, tf = motif_tf[match(id, motif_id)], L = ncol(pwm),
             rel_diff = round(rel_diff[id], 3), minScore = signif(mn, 3),
             maxScore = mx, thr_matchPWM = thr_matchpwm, thr_formula = thr_formula,
             hits = length(hits), hits_ref = length(ref),
             min_hit_score = round(if (length(hits)) min(mcols(hits)$score) else NA, 4),
             warnings = n_warn, row.names = NULL)
}
cat("Verifica soglia 90% su", length(test_subj), "bp:\n")
print(do.call(rbind, lapply(c("MA0069.1", "MA1930.2", "MA1483.3"), check_threshold)))

# ---- Passo 6b: test sintetico dei due filamenti -----------------------------
# Un sito "in avanti" lo trova la PWM originale; uno sull'altro filamento (il
# reverse complement del consenso) lo trova reverseComplement(pwm).
consensus_of <- function(pwm)       # base col peso massimo in ogni posizione
  paste(rownames(pwm)[max.col(t(pwm), ties.method = "first")], collapse = "")

# Sequenza casuale (seed fisso) con stringhe note scritte in posizioni note.
plant <- function(strings, positions, n = 300, seed = 42) {
  set.seed(seed)
  s <- paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
  for (i in seq_along(strings))
    substr(s, positions[i], positions[i] + nchar(strings[i]) - 1) <- strings[i]
  DNAString(s)
}

test_two_strands <- function(id, p_fwd, p_rev) {
  pwm <- pwms[[id]]; tf <- motif_tf[match(id, motif_id)]
  L <- ncol(pwm); cons <- consensus_of(pwm)
  rc_cons <- as.character(reverseComplement(DNAString(cons)))
  subj <- plant(c(cons, rc_cons), c(p_fwd, p_rev))
  hit_fwd <- start(matchPWM(pwm, subj, min.score = "90%"))
  hit_rev <- start(matchPWM(reverseComplement(pwm), subj, min.score = "90%"))
  # Alternativa: PWM originale sul reverse complement dell'intera sequenza.
  hit_rc <- start(matchPWM(pwm, reverseComplement(subj), min.score = "90%"))
  back   <- sort(length(subj) - (hit_rc + L - 1L) + 1L)   # coordinate originali
  stopifnot(abs(PWMscoreStartingAt(pwm, DNAString(cons)) - 1) < 1e-9,
            cons != rc_cons,                              # non palindromo
            identical(hit_fwd, as.integer(p_fwd)), identical(hit_rev, as.integer(p_rev)),
            length(hit_rc) == length(hit_rev), identical(back, hit_rev))
  cat(sprintf("%s (L=%d) consenso %s\n", tf, L, cons))
  data.frame(motif = tf, planted_fwd = p_fwd, planted_rev = p_rev,
             found_fwd = hit_fwd, found_rev = hit_rev, hits_on_rc_seq = length(hit_rc),
             rc_hits_mapped_back = back)
}
cat("Test sintetico (n = 300 bp, soglia 90%):\n")
res <- rbind(test_two_strands("MA0069.1", 41, 181),
             test_two_strands("MA1930.2", 61, 201))
print(res, row.names = FALSE)

# Motivo palindromo (GAATTC e' uguale al suo reverse complement): un solo sito
# piantato viene trovato da ENTRAMBI i filamenti, quindi contato due volte.
pal <- PWM(DNAStringSet(rep("GAATTC", 20)))
pal_subj <- plant("GAATTC", 100)
pal_fwd <- start(matchPWM(pal, pal_subj, min.score = "90%"))
pal_rev <- start(matchPWM(reverseComplement(pal), pal_subj, min.score = "90%"))
stopifnot(isTRUE(all.equal(pal, reverseComplement(pal))),
          identical(pal_fwd, 100L), identical(pal_rev, 100L))
cat("Palindromo GAATTC: 1 sito piantato, trovato in avanti:", length(pal_fwd),
    "| sul filamento opposto:", length(pal_rev), "| contati in totale:",
    length(pal_fwd) + length(pal_rev), "\n")

# ---- Passo 6c: promotori di seguito, assegnazione ai geni, scarto dei confini --
# Tutti i promotori sono larghi W = 1100 e uniti senza spazi: il promotore k
# occupa le posizioni (k-1)*W+1 ... k*W. Quindi basta il resto della divisione:
#   promotore = (start-1) %/% W + 1;  posizione nel promotore = (start-1) %% W;
# il sito (largo L) sta tutto dentro se (start-1) %% W + L <= W.
# (findInterval servirebbe solo con larghezze diverse.)
W <- 1100L
stopifnot(all(width(seqs) == W))
all_seq <- unlist(seqs)                       # 35150 * 1100 bp, un'unica sequenza

# Scansiona i due filamenti e ritorna il promotore (indice di riga) di ogni
# sito tenuto; l'attributo n_discarded conta i siti a cavallo, scartati.
scan_motif <- function(pwm, subject = all_seq) {
  L <- ncol(pwm)
  st <- c(start(matchPWM(pwm, subject, min.score = "90%")),
          start(matchPWM(reverseComplement(pwm), subject, min.score = "90%")))
  keep <- (st - 1L) %% W + L <= W
  structure((st[keep] - 1L) %/% W + 1L, n_discarded = sum(!keep))
}

# Test sintetico: 3 promotori da 1100 bp, PAX6 (L=14), tre siti piantati.
#   2194: attraversa il confine tra promotore 2 e 3 -> scartato
#   1101: inizia alla prima base del promotore 2     -> tenuto, promotore 2
#   3287: finisce all'ultima base del promotore 3 (sul filamento opposto) -> promotore 3
pax6 <- pwms[["MA0069.1"]]; cons <- consensus_of(pax6)
syn <- plant(c(cons, cons, as.character(reverseComplement(DNAString(cons)))),
             c(2194, 1101, 3287), n = 3 * W)
syn_hits <- scan_motif(pax6, syn)
stopifnot(identical(sort(as.integer(syn_hits)), c(2L, 3L)), attr(syn_hits, "n_discarded") == 1)
cat("Test sintetico: siti tenuti nei promotori", sort(syn_hits), "| scartati a cavallo:",
    attr(syn_hits, "n_discarded"), "\n")

# Confronto su ~200 promotori reali: concatenazione contro un promotore alla volta.
test_ids <- c("MA0069.1", "MA1930.2", "MA1483.3")
full <- lapply(test_ids, function(id) scan_motif(pwms[[id]]))   # su tutti i promotori
# Campione: 150 promotori a caso + fino a 17 con almeno un sito PER OGNI motivo
# (cosi' anche i motivi rari, come CTCF, hanno casi non banali).
set.seed(7)
pick <- function(g) { u <- unique(g); u[sample.int(length(u), min(17L, length(u)))] }
idx <- sort(unique(c(sample(length(seqs), 150), unlist(lapply(full, pick)))))
sub <- seqs[idx]; sub_cat <- unlist(sub)
cmp <- do.call(rbind, lapply(test_ids, function(id) {
  pwm <- pwms[[id]]
  by_cat <- tabulate(scan_motif(pwm, sub_cat), nbins = length(sub))
  by_one <- vapply(seq_along(sub), function(i)    # un promotore alla volta
    length(matchPWM(pwm, sub[[i]], min.score = "90%")) +
    length(matchPWM(reverseComplement(pwm), sub[[i]], min.score = "90%")), 1L)
  stopifnot(identical(by_cat, by_one))
  data.frame(motif = motif_tf[match(id, motif_id)], promoters = length(sub),
             genes_with_site = sum(by_cat > 0), sites_concat = sum(by_cat),
             sites_one_by_one = sum(by_one), identical = identical(by_cat, by_one),
             discarded_crossing = attr(scan_motif(pwm, sub_cat), "n_discarded"))
}))
print(cmp, row.names = FALSE)

# ---- Passo 7: prova su 10 motivi (nessun file viene salvato) ----------------
# Il piu' corto, il piu' lungo, ELF2 (rel_diff alto) e 7 a caso, uno per ogni
# fascia di lunghezza (le fasce hanno lo stesso numero di motivi).
set.seed(2026)
n_genes <- length(seqs)
fixed  <- c(which.min(w), which.max(w), match("MA1483.3", motif_id))
rest   <- setdiff(order(w), fixed)                     # per larghezza crescente
strata <- split(rest, cut(seq_along(rest), 7, labels = FALSE))
trial  <- c(fixed, vapply(strata, function(s) s[sample.int(length(s), 1)], 1L))
stopifnot(length(trial) == 10, !anyDuplicated(trial))

invisible(scan_motif(pwms[[trial[1]]]))                # riscaldamento, non cronometrato
trial_res <- do.call(rbind, lapply(trial, function(k) {
  secs <- system.time(g <- scan_motif(pwms[[k]]))[["elapsed"]]   # SOLO scan_motif()
  data.frame(id = motif_id[k], tf = motif_tf[k], L = w[[k]], rescaled = rescaled[[k]],
             rel_diff = round(rel_diff[[k]], 3), sites = length(g),
             frac_genes = length(unique(g)) / n_genes, seconds = round(secs, 2))
}))
trial_res$flag <- ifelse(trial_res$frac_genes > 0.90, "ALTA (>90%)",
                  ifelse(trial_res$frac_genes < 0.005, "bassa (<0.5%)", ""))
trial_res$frac_genes <- round(trial_res$frac_genes, 4)
print(trial_res[order(trial_res$L), ], row.names = FALSE)

# Stima del tempo totale: il costo cresce con L, quindi retta seconds ~ L,
# applicata alla larghezza di ognuno dei 1019 motivi (intervallo di confidenza 95%
# sulla retta, sommato: i parametri sono gli stessi per tutti i motivi).
fit <- lm(seconds ~ L, data = trial_res)
tot <- colSums(predict(fit, data.frame(L = w), interval = "confidence"))
cat(sprintf("Retta: %.2f s + %.3f s per posizione (R2 = %.2f)\n", coef(fit)[1], coef(fit)[2],
            summary(fit)$r.squared))
cat(sprintf("Stima per %d motivi: %.1f min (intervallo %.1f - %.1f min) | media x N: %.1f min\n",
            length(w), tot[["fit"]] / 60, tot[["lwr"]] / 60, tot[["upr"]] / 60,
            mean(trial_res$seconds) * length(w) / 60))

# ---- Passo 8a: scansione di tutti i motivi ----------------------------------
# Per ogni motivo: geni con almeno un sito (unique) e conteggio per gene
# (tabulate, sommato in total_sites). I match si scartano subito dopo.
# total_sites_noext e' uguale ma senza i motivi con frazione di geni > 90%
# ("estesi": troppo poco selettivi). Il risultato grezzo va in /tmp (vedi sotto).
scan_all <- function(ids, every = 100L) {
  hits <- vector("list", length(ids)); sites <- integer(length(ids))
  total <- total_noext <- numeric(n_genes); t0 <- Sys.time()
  for (i in seq_along(ids)) {
    g <- scan_motif(pwms[[ids[i]]])
    hits[[i]] <- sort(unique(g)); sites[i] <- length(g)
    counts <- tabulate(g, nbins = n_genes)
    total <- total + counts
    if (length(hits[[i]]) / n_genes <= 0.90) total_noext <- total_noext + counts
    if (i %% every == 0)
      cat(sprintf("  %d/%d motivi, %.1f min\n", i, length(ids),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  ext <- lengths(hits) / n_genes > 0.90
  stopifnot(sum(total) == sum(sites),               # ogni sito conta una volta
            sum(total_noext) == sum(sites[!ext]))
  list(hits = hits, sites = sites, total_sites = total, total_sites_noext = total_noext,
       motif_id = motif_id[ids], n_genes = n_genes,
       secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

# Con TFBALANCE_REUSE=1 e il file presente si rilegge il risultato senza scansionare.
scan_file <- Sys.getenv("TFBALANCE_SCAN_FILE", "/tmp/tfbalance_scan.rds")
reuse <- Sys.getenv("TFBALANCE_REUSE") == "1" && file.exists(scan_file)

if (!reuse) {
  # Prova su 20 motivi: i 10 del passo 7 + 10 sparsi lungo le lunghezze.
  ids20 <- c(trial, head(setdiff(order(w)[round(seq(1, length(w), length.out = 30))], trial), 10))
  stopifnot(length(ids20) == 20, !anyDuplicated(ids20))
  r20 <- scan_all(ids20, every = 1e9L)
  cat(sprintf("Prova su 20 motivi: %.1f s misurati | %.1f s previsti dalla retta del passo 7\n",
              r20$secs, sum(predict(fit, data.frame(L = w[ids20])))))
  # Coerenza con il passo 7: stessi 10 motivi, stessi siti e stessa frazione di geni.
  stopifnot(identical(r20$sites[1:10], trial_res$sites),
            all(round(lengths(r20$hits[1:10]) / n_genes, 4) == trial_res$frac_genes))
  if (Sys.getenv("TFBALANCE_TRIAL_ONLY") == "1") {
    cat("Solo prova (TFBALANCE_TRIAL_ONLY=1): mi fermo qui.\n"); quit(save = "no", status = 0)
  }
  cat("Ciclo completo su", length(pwms), "motivi:\n")
  all_res <- scan_all(seq_along(pwms))
  saveRDS(all_res, scan_file)                      # risultato grezzo, fuori dal progetto
  cat(sprintf("Fatto in %.1f min. Risultato grezzo salvato in %s\n", all_res$secs / 60, scan_file))
} else {
  all_res <- readRDS(scan_file)
  stopifnot(identical(all_res$motif_id, motif_id), all_res$n_genes == n_genes,
            length(all_res$hits) == length(pwms), length(all_res$sites) == length(pwms),
            length(all_res$total_sites) == n_genes, length(all_res$total_sites_noext) == n_genes,
            sum(all_res$total_sites) == sum(all_res$sites),
            sum(all_res$total_sites_noext) ==
              sum(all_res$sites[lengths(all_res$hits) / n_genes <= 0.90]))
  # Stesso numero di motivi e di geni non basta: il file potrebbe venire da promotori
  # o PWM diversi. Riscansiono 3 motivi (corto, lungo, ELF2: pochi secondi) e confronto.
  for (k in trial[1:3]) {
    g <- scan_motif(pwms[[k]])
    stopifnot(identical(sort(unique(g)), all_res$hits[[k]]), length(g) == all_res$sites[k])
  }
  cat("Risultato grezzo riletto da", scan_file, "(controllato su 3 motivi, nessuna scansione completa)\n")
}

# ---- Riassunti ---------------------------------------------------------------
frac <- lengths(all_res$hits) / n_genes           # frazione di geni con >=1 sito
fb <- cut(frac, c(-Inf, 0.005, 0.10, 0.50, 0.90, Inf), right = FALSE,
          labels = c("<0.5%", "0.5-10%", "10-50%", "50-90%", ">90%"))
cat("Motivi per frazione di geni con >=1 sito (estremo sinistro incluso):\n"); print(table(fb))
cat("Mediana della frazione per fascia di rel_diff:\n")
print(data.frame(rel_diff = levels(band), motivi = as.vector(table(band)),
                 mediana_frazione = round(as.vector(tapply(frac, band, median)), 4)),
      row.names = FALSE)
top10 <- order(all_res$sites, decreasing = TRUE)[1:10]
print(data.frame(id = motif_id[top10], tf = motif_tf[top10], L = w[top10],
                 siti = all_res$sites[top10],
                 pct_total_sites = round(100 * all_res$sites[top10] / sum(all_res$sites), 1)),
      row.names = FALSE)
cat(sprintf("I 10 motivi con piu' siti: %.1f%% di tutti i siti (%d su %d)\n",
            100 * sum(all_res$sites[top10]) / sum(all_res$sites), sum(all_res$sites[top10]),
            sum(all_res$sites)))
cat(sprintf("Correlazione gc ~ log1p(total_sites): Pearson %.3f, Spearman %.3f\n",
            cor(gc, log1p(all_res$total_sites)),
            cor(gc, log1p(all_res$total_sites), method = "spearman")))
cat(sprintf("Motivi con frazione > 90%% esclusi da total_sites_noext: %d su %d\n",
            sum(frac > 0.90), length(frac)))
cat(sprintf("Correlazione gc ~ log1p(total_sites_noext): Pearson %.3f, Spearman %.3f\n",
            cor(gc, log1p(all_res$total_sites_noext)),
            cor(gc, log1p(all_res$total_sites_noext), method = "spearman")))

# ---- Passo 8b: costruzione e salvataggio dei due file in data/ ---------------
out_dir <- Sys.getenv("TFBALANCE_OUT_DIR", "data")        # di norma "data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
as_count <- function(x) { stopifnot(!anyNA(x), all(x == round(x))); as.integer(x) }

# Una riga per promotore, nello stesso ordine di p (quindi di seqs, gc e indici).
# symbol ed ensembl possono essere NA (geni senza annotazione): e' voluto.
genes_df <- data.frame(entrez = names(p), symbol = unname(symbol), ensembl = unname(ensembl),
                       gc = gc, total_sites = as_count(all_res$total_sites),
                       total_sites_noext = as_count(all_res$total_sites_noext),
                       stringsAsFactors = FALSE)
# Un elemento per motivo, nell'ordine del file JASPAR; genes = indici di riga di genes_df.
motif_hits <- lapply(seq_along(pwms), function(k)
  list(id = motif_id[k], tf = motif_tf[k], L = w[[k]], genes = all_res$hits[[k]],
       rescaled = rescaled[[k]], rel_diff = rel_diff[[k]],
       frac = length(all_res$hits[[k]]) / n_genes))
names(motif_hits) <- motif_id

# Ordine dei motivi riletto dal file JASPAR, indipendente da motif_id.
jaspar_ids <- sub("^>", "", vapply(strsplit(grep("^>", readLines(jaspar_file), value = TRUE),
                                            "\t"), `[`, "", 1))
all_idx <- unlist(lapply(motif_hits, `[[`, "genes"), use.names = FALSE)
stopifnot(nrow(genes_df) == length(p), nrow(genes_df) == 35150, length(motif_hits) == 1019,
          identical(names(motif_hits), jaspar_ids),
          identical(unname(vapply(motif_hits, `[[`, "", "id")), jaspar_ids),
          all(vapply(motif_hits, function(m) !is.unsorted(m$genes, strictly = TRUE), TRUE)),
          all(all_idx >= 1L & all_idx <= nrow(genes_df)),
          !anyNA(genes_df$gc), !anyNA(genes_df$total_sites), !anyNA(genes_df$total_sites_noext),
          !anyDuplicated(genes_df$entrez),
          all(genes_df$total_sites_noext <= genes_df$total_sites),
          all(tabulate(all_idx, nrow(genes_df)) <= genes_df$total_sites),  # >=1 sito per motivo
          sum(genes_df$total_sites) == sum(all_res$sites),
          genes_df$symbol[genes_df$entrez == "7157"] == "TP53")

f_genes <- file.path(out_dir, "genes.rds"); f_hits <- file.path(out_dir, "motif_hits.rds")
saveRDS(genes_df, f_genes); saveRDS(motif_hits, f_hits)
stopifnot(identical(readRDS(f_genes), genes_df), identical(readRDS(f_hits), motif_hits))
cat(sprintf("%s: %.2f MB\n", c(f_genes, f_hits), file.size(c(f_genes, f_hits)) / 1e6), sep = "")
