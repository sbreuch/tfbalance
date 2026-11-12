# tfbalance

Shiny app for transcription-factor (TF) binding-site enrichment on human
promoters, with a test that corrects for GC content and site density.

**Status: work in progress.**

## Idea

A standard enrichment test asks whether the genes in a list have more predicted
binding sites for a TF than expected. Promoters differ in GC content and in the
total number of predicted sites, which can inflate enrichment for some TFs.
tfbalance (planned) shows two tests side by side:

- baseline: hypergeometric test on "genes with at least one site"
- corrected: logistic regression `in_list ~ has_site + gc + log1p(total_sites)`

P-values are adjusted with Benjamini-Hochberg.

## Data preparation

`data-raw/prepare_data.R` runs once, outside the app, and writes
`data/genes.rds` and `data/motif_hits.rds`. The numbers below are from its first
complete run. Two of them, 35,356 (genes returned by `genes()`) and 115 (motifs
with a palindromic consensus), were measured separately and are not printed by
the script.

**Genes and promoters** (GRCh38, UCSC hg38 knownGene)

- `genes()` dropped 2,169 genes whose exons lie on both strands or on several
  chromosomes, and returned 35,356; 35,161 of them are on chr1-22, X, Y.
- One promoter per gene: 1000 bp upstream to 100 bp downstream of the TSS
  (1100 bp, strand respected), using the most upstream TSS of the gene. None
  fell outside its chromosome.
- 11 promoters containing at least one `N` were removed, leaving **35,150
  promoters**. Of these, 21 have no gene symbol and 7,752 have no Ensembl ID
  (kept, with `NA`).

**Motifs and scan**

- 1019 JASPAR 2026 CORE vertebrates non-redundant motifs, 4 to 33 positions
  (median 9).
- `Biostrings::PWM` needs equal column sums, but in 528 motifs the JASPAR column
  sums differ. In these, each column was rescaled to the mean column sum of its
  motif (largest change of a base proportion: 0.017). In 409 of the 528 the
  columns differ by more than 10% (`rel_diff` = (max - min) / max of the column
  sums); the other 119 differ by 10% or less.
- Scan: `Biostrings::matchPWM`, minimum score 90% of the maximum, both strands
  (`reverseComplement(pwm)` for the opposite one). Promoters are scanned
  concatenated and sites straddling two promoters are discarded. 57,009,554
  sites in total.
- A palindromic site is found on both strands and counted twice in
  `total_sites` (115 motifs have a palindromic consensus). It does not change
  "genes with at least one site".

**Fraction of genes with at least one site, per motif**

| < 0.5% | 0.5-10% | 10-50% | 50-90% | > 90% |
|-------:|--------:|-------:|-------:|------:|
|     58 |     243 |    383 |    226 |   109 |

The 109 motifs above 90% are not selective at this threshold (for example, the
4-position ZGLP1 matches 99.98% of promoters), and the 58 below 0.5% (for
example CTCF, 16 genes) leave little statistical power.

**Output**

- `data/genes.rds` (0.54 MB): 35,150 rows, columns `entrez`, `symbol`,
  `ensembl`, `gc`, `total_sites`, `total_sites_noext`. `total_sites` sums the
  sites of all motifs on both strands; `total_sites_noext` leaves out the 109
  motifs above 90%.
- `data/motif_hits.rds` (31 MB): list of 1019 elements in the order of the
  JASPAR file, each with `id`, `tf`, `L`, `genes`, `rescaled`, `rel_diff`,
  `frac`. `genes` are row indices of `genes.rds`, so its rows must not be
  reordered.
- `gc` and `total_sites` are strongly related: correlation of `gc` with
  `log1p(total_sites)` is -0.749 (Spearman -0.772), and -0.185 (Spearman -0.212)
  with `log1p(total_sites_noext)`.

**How to run**

The JASPAR file is not stored in this repository. Download it with:

    curl -L -o data-raw/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt \
      "https://jaspar.elixir.no/download/data/2026/CORE/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt"

Then, from the project root, run the script once (about 20 minutes, of which
16 are the scan; run a single R process at a time):

    Rscript data-raw/prepare_data.R

The raw scan result is cached in `/tmp/tfbalance_scan.rds`. With
`TFBALANCE_REUSE=1` the script reloads it instead of scanning again. The header
of the script lists the other `TFBALANCE_*` variables.

## Requirements

- App: R with `shiny`, `bslib`, `DT`, `ggplot2`
- Data preparation only: Bioconductor packages `GenomicFeatures`, `Biostrings`,
  `BSgenome`, `AnnotationDbi`, `org.Hs.eg.db`, `TxDb.Hsapiens.UCSC.hg38.knownGene`,
  `BSgenome.Hsapiens.UCSC.hg38`

## Assumptions and limitations

- The background is the set of promoters kept after filtering; it includes
  non-coding genes. Genes with exons on several strands or chromosomes are
  excluded. In the first run the background is 35,150 promoters, and `genes()`
  excluded 2,169 genes.
- Palindromic sites are counted once per strand, so twice in total.
- The scanning threshold follows the 0.90 relative score used in the literature,
  but the pipeline is not identical to the original published tools.

## References

To be added.
