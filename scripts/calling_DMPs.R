calling_DMPs <- function(gr, min_cov = 6, fdr = 0.05) {
    gr <- gr[which(gr$readsN1 >= min_cov & gr$readsN2 >= min_cov)]
    gr$readsU1 <- gr$readsN1 - gr$readsM1
    gr$readsU2 <- gr$readsN2 - gr$readsM2
    gr$proportion1 <- gr$readsM1 / gr$readsN1
    gr$proportion2 <- gr$readsM2 / gr$readsN2
    gr$delta <- gr$proportion2 - gr$proportion1

    gr <- gr[which(
        (gr$context == "CG" & abs(gr$delta) >= 0.4) |
            (gr$context == "CHG" & abs(gr$delta) >= 0.2) |
            (gr$context == "CHH" & abs(gr$delta) >= 0.1)
    )]

    gr$pvalue <- apply(as.data.frame(mcols(gr)[, c("readsM2", "readsU2", "readsM1", "readsU1")]), 1, function(v) {
        fisher.test(matrix(c(v[1], v[2], v[3], v[4]), nrow = 2, byrow = TRUE))$p.value
    })

    gr$direction <- ifelse(gr$delta > 0, "gain", "loss")

    if (is.null(fdr)) {
        gr <- gr[which(gr$pvalue <= 0.05)]
    } else {
        gr$padj <- p.adjust(gr$pvalue, method = "BH")
        gr <- gr[which(gr$padj <= fdr)]
    }

    return(gr)
}
