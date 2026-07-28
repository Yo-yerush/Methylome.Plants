DMRs_into_groups <- function(treatment, ann, context = "all",
                             DMRs_ann_dir = "./genome_annotation",
                             datasets_dir = NULL,
                             gene_sets_path = NULL) {
  if (is.null(gene_sets_path)) {
    stop("Functional groups require functional.gene_sets or a bundle-specific adapter.")
  }
  DMRs_into_generic_groups(
    treatment = treatment,
    ann = ann,
    context = context,
    gene_sets_path = gene_sets_path,
    DMRs_ann_dir = DMRs_ann_dir
  )
}

groups_barPlots <- function(x) {
  x %>%
    ggplot(aes(x = reorder(group, total), y = n, fill = regionType)) +
    geom_col() +
    geom_text(aes(y = total, label = total), hjust = -0.3, size = 3, check_overlap = TRUE) +
    facet_wrap(~context, ncol = 4, scales = "free") +
    coord_flip() +
    theme_bw() +
    scale_fill_manual(values = c(gain = "#d96c6c", loss = "#6c96d9")) +
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.15))) +
    labs(x = NULL, y = "Number of DMR-annotated genes", fill = "DMR direction")
}
