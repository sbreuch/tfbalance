# R/compare.R
# Dati e grafico per confrontare il test baseline e il test corretto, e i conteggi dei motivi
# significativi. NON contiene statistica nuova (il BH sul sottoinsieme e' bh_subset in
# R/enrich.R): qui si trasformano e si selezionano i numeri della tabella di analyze_list.
# Richiede che R/enrich.R sia gia' caricato e, per compare_plot, ggplot2 (library(ggplot2)).

# Livello di FDR: un motivo e' "significativo" se il suo FDR e' < 0.05.
FDR_LEVEL <- 0.05

# baseline_fdr_sub (uso interno): FDR del baseline calcolato solo sui motivi che hanno anche il
# test corretto (z non NA), cioe' nella stessa famiglia di motivi del corretto.
baseline_fdr_sub <- function(table) bh_subset(table$p_baseline, !is.na(table$z))

# compare_data: i punti del grafico, uno per motivo con z non NA.
#   x = -log10(p baseline), y = -log10(p corretto).
# Un p troppo piccolo (p = 0 o sotto il piu' piccolo numero che R sa scrivere, circa 2.2e-308) ha
# -log10 infinito o inaffidabile: quella coordinata si disegna a un tetto `cap` uguale per i due
# assi, 1.1 volte il massimo dei valori veri piu' 1 (il +1 serve se non ci sono valori veri).
# x_clipped / y_clipped dicono quale coordinata e' stata messa al tetto.
compare_data <- function(table) {
  t <- table[!is.na(table$z), , drop = FALSE]       # solo i motivi con il test corretto
  ux <- t$p_baseline < .Machine$double.xmin         # p troppo piccolo: underflow
  uy <- t$p_corrected < .Machine$double.xmin
  x <- -log10(t$p_baseline)
  y <- -log10(t$p_corrected)
  cap <- 1.1 * max(c(x[!ux], y[!uy], 0)) + 1
  x[ux] <- cap
  y[uy] <- cap
  list(points = data.frame(id = t$id, tf = t$tf, x = x, y = y,
                           x_clipped = ux, y_clipped = uy, stringsAsFactors = FALSE),
       cap = cap)
}

# fdr_threshold: la soglia in scala p. E' il p piu' grande tra i motivi con FDR < level; NA se
# nessun motivo e' significativo. Siccome l'FDR di BH non diminuisce al crescere di p, i motivi
# con FDR < level sono esattamente quelli con p <= soglia.
fdr_threshold <- function(p, fdr, level = FDR_LEVEL) {
  sig <- !is.na(p) & !is.na(fdr) & fdr < level
  if (any(sig)) max(p[sig]) else NA_real_
}

# fdr_counts: quanti motivi hanno FDR < 0.05.
#   base_all  baseline, su tutti i motivi della tabella (e' l'fdr_baseline che si vede in tabella)
#   base_sub  baseline, ricalcolato solo sui motivi con test corretto (stessa famiglia del corretto)
#   corrected corretto (il suo FDR e' gia' calcolato su quegli stessi motivi)
#   n_all, n_sub  quanti motivi ci sono nelle due famiglie
# base_sub e corrected sono NA (non 0) se nessun motivo ha il test corretto.
fdr_counts <- function(table) {
  has_z <- !is.na(table$z)
  list(
    base_all = sum(table$fdr_baseline < FDR_LEVEL, na.rm = TRUE),
    base_sub = if (any(has_z)) sum(baseline_fdr_sub(table) < FDR_LEVEL, na.rm = TRUE)
               else NA_integer_,
    corrected = if (any(has_z)) sum(table$fdr_corrected < FDR_LEVEL, na.rm = TRUE)
                else NA_integer_,
    n_all = nrow(table),
    n_sub = sum(has_z)
  )
}

# compare_plot: il grafico di confronto (un oggetto ggplot).
#   - un punto per motivo con z non NA; i punti troncati hanno forma di triangolo;
#   - diagonale y = x tratteggiata: sopra la diagonale il test corretto e' piu' significativo;
#   - linee a FDR 0.05 (scala p): orizzontale per il corretto, verticale per il baseline,
#     solo se quel test ha motivi significativi. Il baseline usa l'FDR sul sottoinsieme, cosi'
#     la soglia appartiene alla stessa famiglia dei punti disegnati;
#   - assi con gli stessi limiti e coord_fixed(1): la diagonale e' a 45 gradi.
compare_plot <- function(table) {
  d <- compare_data(table)
  pts <- d$points
  pts$clip <- pts$x_clipped | pts$y_clipped
  has_z <- !is.na(table$z)
  p_x <- fdr_threshold(table$p_baseline[has_z], baseline_fdr_sub(table)[has_z])
  p_y <- fdr_threshold(table$p_corrected[has_z], table$fdr_corrected[has_z])
  # posizione di una soglia in scala -log10; se la soglia e' un p in underflow sta al tetto
  pos <- function(p) if (p < .Machine$double.xmin) d$cap else -log10(p)

  g <- ggplot(pts, aes(x, y, shape = clip)) +
    geom_point(colour = "#2c7fb8", alpha = 0.7, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50")
  if (!is.na(p_y)) g <- g + geom_hline(yintercept = pos(p_y), colour = "firebrick")
  if (!is.na(p_x)) g <- g + geom_vline(xintercept = pos(p_x), colour = "firebrick")

  g <- g +
    scale_shape_manual(name = NULL, values = c(`FALSE` = 16, `TRUE` = 17),
                       labels = c(`FALSE` = "p calcolato",
                                  `TRUE` = "p = 0 (underflow): disegnato al limite del grafico")) +
    scale_x_continuous(limits = c(0, d$cap), expand = expansion(mult = 0.02)) +
    scale_y_continuous(limits = c(0, d$cap), expand = expansion(mult = 0.02)) +
    coord_fixed(1) +
    labs(x = "-log10 p baseline", y = "-log10 p corretto") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom")
  if (!any(pts$clip)) g <- g + guides(shape = "none")    # niente legenda se non c'e' nulla da spiegare
  g
}
