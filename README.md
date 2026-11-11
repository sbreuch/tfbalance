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

## Data preparation (planned, not run yet)

- Genome and genes: GRCh38 (UCSC hg38 knownGene), one promoter per gene
- Promoter: 1000 bp upstream and 100 bp downstream of the TSS, strand respected
- Motifs: JASPAR 2026 CORE vertebrates, non-redundant, PFM format
- Scan: `Biostrings::matchPWM`, minimum score 90%, both strands
- Output: `data/genes.rds` and `data/motif_hits.rds`

The JASPAR file is not stored in this repository. Download it with:

    curl -L -o data-raw/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt \
      "https://jaspar.elixir.no/download/data/2026/CORE/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt"

Then run the preparation script once (it will be added to `data-raw/`):

    Rscript data-raw/prepare_data.R

## Requirements

- App: R with `shiny`, `bslib`, `DT`, `ggplot2`
- Data preparation only: Bioconductor packages `GenomicFeatures`, `Biostrings`,
  `BSgenome`, `AnnotationDbi`, `org.Hs.eg.db`, `TxDb.Hsapiens.UCSC.hg38.knownGene`,
  `BSgenome.Hsapiens.UCSC.hg38`

## Assumptions and limitations

- The background is the set of promoters kept after filtering; it includes
  non-coding genes. Genes with exons on several strands or chromosomes are
  excluded. Exact numbers will be added after the run.
- Palindromic sites are counted once per strand.
- The scanning threshold follows the 0.90 relative score used in the literature,
  but the pipeline is not identical to the original published tools.

## References

To be added.
