# tfbalance

Shiny app for transcription-factor (TF) binding-site enrichment on human
promoters, with a test that corrects for GC content and site density.

**Status: work in progress.** The data preparation, the statistics (`R/enrich.R`),
the analysis logic (`R/analyze.R`) and a first version of the app (`app.R`) are
written and tested; the statistics are also checked on real data (see "Statistics"
and "Validation and known limits"). The app runs locally and is not published. See
"The app" for what it does today and what is planned.

## Idea

A standard enrichment test asks whether the genes in a list have more predicted
binding sites for a TF than expected. Promoters differ in GC content and in the
total number of predicted sites, which can inflate enrichment for some TFs.
tfbalance shows two tests side by side:

- baseline: one-sided hypergeometric test on "genes with at least one site"
- corrected: one-sided test of `has_site` in the logistic model
  `in_list ~ has_site + gc + log1p(total_sites_noext)`, computed as a score
  test (see "Statistics"). The model uses `total_sites_noext`, not
  `total_sites`.

Both only look for enrichment, not for depletion. P-values are adjusted with
Benjamini-Hochberg.

## The app

Paste a list of human genes and press **Analizza** (the interface is in Italian).
Genes can be symbols, Entrez IDs or Ensembl IDs (also `ENSG....7`), separated by
spaces, tabs, commas, semicolons or new lines. Both tests are computed on the 852
motifs of the default filter and shown side by side, one row per motif.

**Current state**

- Tab *Analisi*: text box and button; a summary box (recognised genes, repeated or
  ignored tokens, unrecognised names, motifs tested and excluded by the filter,
  motifs with the corrected test, and a note when the corrected test is not
  available); a table (DT) with `k`, `K`, expected genes, `z` and p and FDR of both
  tests, sorted by the corrected p-value (`NA` last), with a caption on the `NA`
  values. Errors, for example no gene recognised, are shown in red.
- Tab *Geni non riconosciuti*: number and list of the names not recognised.
- Tab *Come funziona*: placeholder, the text is still to be written.
- Planned, not in the app yet: a comparison plot (its data and the FDR counts are
  in `R/compare.R`, written and tested but not used by `app.R`) and demo buttons
  (`random_list` exists in `R/enrich.R`).
- The app always uses the default filter and `min_expected = 10`: there is no
  control for them.

**Run it**, from the project root (`data/` is in the repository):

    Rscript -e 'shiny::runApp()'

Shiny loads the files in `R/` by itself. It needs `shiny`, `bslib` and `DT`, plus
`ggplot2` for `R/compare.R`.

**Code:** `app.R` (interface and server); `R/analyze.R` (`analyze_list`: parsing,
limits, motif filter, both tests, sorting); `R/compare.R` (plot data, thresholds and
FDR counts); `R/enrich.R` (all the statistics).

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

## Statistics

All the statistics are in `R/enrich.R`, as plain functions. They read no files:
they receive `genes` (the data frame of `genes.rds`) and `hits` (the list of
`motif_hits.rds`). Lists of genes are passed as row indices of `genes`.

| Function | What it does |
|---|---|
| `map_ids(x, genes)` | symbols, Entrez or Ensembl IDs (also `ENSG....7`), any case, to row indices; also returns the names not found |
| `select_motifs(hits, n_genes, frac_range = c(0.005, 0.9))` | TRUE/FALSE per motif: drops the 58 with fewer than 0.5% and the 109 with more than 90% of genes with a site (852 stay) |
| `enrich_baseline(idx, hits, n_genes)` | hypergeometric test: `id, tf, k, K, expected, p, fdr` |
| `enrich_corrected(idx, hits, genes, min_expected = 10)` | score test: `id, tf, expected, z, p, fdr` |
| `bh_subset(p, keep)` | BH computed only on the motifs with `keep = TRUE`, `NA` for the others (used to put the baseline in the same family of motifs as the corrected test) |
| `odds_ratio_one(idx, motif_genes, genes)` | `glm()` for ONE motif: `beta, se, or, or_low, or_high` (Wald 95%), `note` |
| `random_list(genes, n, biased = FALSE, by = "gc", power = 2)` | demo list of `n` genes, uniform or biased (see below) |

Notation: `n` = genes in the list, `K` = genes with a site for the motif,
`k` = genes of the list with a site, `N` = 35,150, `E = n * K / N` = genes of
the list expected to have a site.

- **BH** is computed on the motifs passed to the function (852 with the default
  filter). Motifs that cannot be tested (`NA`) do not count in the number of tests.
- **Why a score test.** `glm()` for each of ~1000 motifs on 35,150 genes takes
  roughly 1.5-2 minutes, too slow for an app. The score test fits the model
  without `has_site` once, then for each motif adds up a few numbers over its
  genes: about 0.5-0.8 s for all motifs. It gives `z` and `p`, not an odds ratio;
  `odds_ratio_one` gives the odds ratio for a single motif (about 0.1 s). On the
  test data `z^2` equals the Rao statistic of `anova(null, full, test = "Rao")`
  (largest relative difference 9e-6).
- **`min_expected`.** A motif with `E < min_expected` gets `z`, `p`, `fdr` = `NA`
  but keeps its row; `E == min_expected` is included; `expected` is always filled.
  `min_expected = 0` turns the rule off. Other `NA` cases: `K = 0` or `K = N`
  (nothing to compare), a relative variance guard (`v / sum(w) < 1e-8`), and a
  null model that cannot be fitted (a warning and all `NA`).
- **`odds_ratio_one` notes** (the strings are in Italian): `"ok"`;
  `"has_site costante"` (`K = 0` or `K = N`); `"separazione o non convergenza"`.
  The last one is returned, with `NA` numbers and without calling `glm`, when
  any cell of the 2x2 table `has_site` x `in_list` is zero (`k`, `n - k`,
  `K - k`, `N - K - (n - k)`): the maximum-likelihood estimate does not exist.
  It is also returned if `glm` does not converge or the fitted probabilities
  are numerically 0 or 1 (same criterion as `stats::glm.fit`). Warning texts
  are never read, because R messages are translated. With 59 of 60 genes of
  the list having the site, no cell is empty and the answer is `"ok"` with a
  huge but existing odds ratio.
- **`random_list`.** With `biased = TRUE` each gene gets the weight
  `(value / mean(value))^power`, where `value` is `gc` or `total_sites_noext`;
  `power > 0` favours high values, `power < 0` low ones. A value `<= 0` counts as
  the smallest positive one. The result is sorted. The default `power = 2` is a
  choice for demos, not the outcome of the calibration (see below).

## Tests

    Rscript -e 'testthat::test_dir("tests/testthat")'

- `tests/testthat/test-enrich.R`: unit tests on small synthetic data (including
  `bh_subset`).
- `tests/testthat/test-enrich-validation.R`: on real data, skipped if
  `data/*.rds` are missing. Baseline against `fisher.test(alternative = "greater")`
  on 20 real motifs; score test against the Wald test of `glm()` on 20 real motifs,
  for a random and an enriched list; `odds_ratio_one` on the motif with the
  highest `z`.
- `tests/testthat/test-analyze.R`: `analyze_list` on synthetic data.
- `tests/testthat/test-compare.R`: `R/compare.R` on synthetic tables, with
  hand-computed values.
- `tests/testthat/test-app.R`: the app server run with `shiny::testServer()`, no
  browser, so the layout is not checked. Needs `data/*.rds`, skipped if missing.

## Validation and known limits

The experiments are in `data-raw/validation/` (scripts `step4*.R`, `step7_timing.R`, `step8_na_with_E.R`;
outputs `*_output.txt`; results `*.rds`). The criteria were written before each run
and not changed afterwards. Real data, the 852 motifs of the default filter, lists
of 300 genes unless stated.

- **Score test against the Wald test** (20 real motifs). Random list: largest
  |z_score - z_wald| = 0.015. Enriched list: Spearman 1.000 and the same decision at
  p < 0.05 for every motif away from the threshold. The planted motif has z = 11.29
  (score) and 10.62 (Wald): the two differ in the far tail, so p-values there should
  not be compared one to one.
- **Lists biased by GC or by `total_sites_noext`** (100 lists for each of: GC or
  density, high or low, `power` 2, 5, 10; plus 300 uniform lists). Fraction of motifs
  with p < 0.05 per list: baseline 0.041 on uniform lists and up to 0.375 with a GC
  bias (`power = -10`); corrected 0.053 on uniform lists and 0.045-0.066 on all biased
  ones. A low-density bias does not inflate the baseline (0.034-0.043): the test only
  sees enrichment. The corrected test keeps a small GC effect (0.066 at `power = 5`,
  0.045 at `power = -5`), because the model assumes a straight-line effect of GC on
  the log-odds.
- **The declared calibration rule failed.** On uniform lists, the fraction of lists with
  at least one motif at FDR < 0.05 was 0.183 for the corrected test (declared bound
  0.15; baseline 0.04). So no default for `power` came out of the rule, and
  `power = 2` is only a demo choice.
- **Few expected genes** (`E`). The hypothesis "low `E` explains that excess" was
  **not confirmed** by its declared criteria (lists with a discovery among motifs with
  `E >= 5`: 0.123, bound 0.09). With `E < 5` the fraction of p < 0.001 was 4.04 times
  the expected one. On 300 new uniform lists and the 728 motifs with `E >= 10`, the
  declared criteria were met: lists with a discovery 0.053 (SE 0.013, bound 0.09);
  fraction of p < 0.001 equal to 0.00132 (expected 0.001, bound 0.002). That is why
  `min_expected = 10`. Known residue: the fraction of p < 0.001 is 2.3 times the
  expected one for `10 <= E < 20` and 1.8 times for `20 <= E < 50` (1.06 for
  `E >= 50`). Over all 852 motifs, lists with a discovery: 0.130 (SE 0.019).
  Below `E = 10` the test is not validated.
- **Small lists.** On 500 random lists of real genes, the null model could not be
  fitted in 1 list of 3 genes, and in none of the lists of 5 and 10 genes. Lists made
  of extreme genes (highest or lowest GC or density) always fail, with a warning and
  all `NA`. With `min_expected = 10`, no motif of the 852 reaches `E >= 10` unless the
  list has at least 12 genes (the largest motif has 31,596 genes).
- **`NA` with `E >= 10`** (`step8_na_with_E.R`: 20 uniform random lists of real genes
  for each of 12, 20, 50, 300 and 5000 genes). With the default filter, none of the
  46,740 (list, motif) pairs with `E >= 10` had `z = NA`. With the filter off (1019
  motifs) the only ones were the three motifs present in all 35,150 genes (`K = N`:
  MA0498.3, MA0625.3, MA1536.2). Lists whose null model cannot be fitted (see "Small
  lists") are not part of this experiment: they give all `NA`, also with `E >= 10`;
  the app then shows a note instead of the warning.
- **Speed** (`enrich_baseline` + `enrich_corrected`, lists of 10 to 5000 genes, 852 or
  1019 motifs): median 0.55-0.90 s in total (slowest single run 0.96 s), of which the
  baseline is about 0.07 s. Motifs with
  a `z` in the default filter: 0 for 10 genes, 488 for 50, 728 for 300, 852 for 5000
  (with the filter off: 0, 594, 834, 989 of 1019).
- **Checked only for** lists of 300 genes drawn at random (and the biased lists above),
  plus the small lists above. Other sizes and other kinds of lists were not checked.

## Requirements

- App: R with `shiny`, `bslib`, `DT`, `ggplot2` (`ggplot2` is used by
  `R/compare.R`, not yet by `app.R`)
- Tests: `testthat` (the app tests also use `jsonlite` and `htmltools`, installed
  with `shiny`)
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
- `total_sites_noext` includes the sites of the motif being tested (for the 910 motifs
  with at most 90% of genes with a site). Their share of `total_sites_noext` is 0.05%
  at the median, 0.11% on average and 1.32% at most (Arid3b). The correlation of
  `has_site` with `log1p(total_sites_noext)` over the 852 motifs has median 0.02, 5th
  to 95th percentile -0.08 to 0.20 and maximum 0.40. The effect is expected to be
  small and, if anything, to make the corrected test slightly conservative. It cannot
  be removed without per-gene counts for each motif, which are not stored.
- Motifs share genes, and some TFs have several matrices (55 of the 1019 names repeat
  a TF already listed). Tests are therefore not independent and the assumptions of
  Benjamini-Hochberg may not hold exactly (see the 0.183 above).
- Only enrichment is tested (one-sided), not depletion.
- The app accepts at most 5000 recognised genes (`MAX_LIST` in `R/analyze.R`, the
  largest size timed in step 7) and 50,000 raw tokens (`MAX_TOKENS`); beyond that it
  shows an error.

## Use of AI

Developed with the help of Claude Code; method choices, test design and code review by the author.

## References

- JASPAR 2026 CORE, vertebrates, non-redundant matrices, downloaded from
  <https://jaspar.elixir.no/download/data/2026/CORE/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt>.
- Everything else: to be added.
