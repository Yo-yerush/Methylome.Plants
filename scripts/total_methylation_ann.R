total_methylation_ann <- function(annotation_vec, var1, var2, meth_var1, meth_var2, context) {
  
  ######################## for CX annotation
  meth_edit <- function(x) {
    x = x[which(x$context == context)]
    x = x[which(x$readsN >= 6)]
    x$proportion = x$readsM / x$readsN
    return(x)
  }
  
  meth_var1_cx = meth_edit(meth_var1)
  meth_var2_cx = meth_edit(meth_var2)
  
  ann_CX_df = data.frame(type = c(names(annotation_vec)),
                         x = rep(NA, length(annotation_vec)),
                         y = rep(NA, length(annotation_vec)))
  
  names(ann_CX_df)[2:3] = c(var2, var1)
  
  ######################## main annotation loop
  for (i in 1:length(annotation_vec)) {
    
    # find overlaps for all CX with annotation file (gff3)
    ## var1 (control)
    m1 <- findOverlaps(meth_var1_cx, annotation_vec[[i]])
    CX_annotation1 <- meth_var1_cx[queryHits(m1)]
    ann_CX_df[i,var1] = length(CX_annotation1)
    
    ## var2 (treatment)
    m2 <- findOverlaps(meth_var2_cx, annotation_vec[[i]])
    CX_annotation2 <- meth_var2_cx[queryHits(m2)]
    ann_CX_df[i,var2] = length(CX_annotation2)
  }
  
  ann_CX_df$foldChange = round(ann_CX_df[,var2] / ann_CX_df[,var1], 4)
  write.csv(ann_CX_df, paste0(context,"_total_methylation_annotations.csv"), row.names = F)
}