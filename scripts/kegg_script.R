run_KEGG <- function(comparison_name, genome_ann_path, KEGG_path, n.cores,
                     reference_bundle = NULL) {
  kegg_organism <- bundle_get(reference_bundle, 'functional.kegg_organism')
  if (is.null(kegg_organism)) stop('KEGG analysis requires functional.kegg_organism in the reference bundle.')
  orgdb_package <- bundle_get(reference_bundle, 'functional.orgdb_package')
  orgdb_keytype <- bundle_get(reference_bundle, 'functional.orgdb_keytype')
  kegg_id_map_path <- bundle_get(reference_bundle, 'functional.kegg_id_map')
  kegg_id_map <- NULL
  if (!is.null(kegg_id_map_path)) {
    kegg_id_map <- data.table::fread(kegg_id_map_path, data.table = FALSE)
    if (!all(c('gene_id', 'kegg_id') %in% names(kegg_id_map))) {
      stop('kegg_id_map must contain gene_id and kegg_id columns.')
    }
  }
  universe_ids <- NULL
  if (!is.null(orgdb_package) && requireNamespace(orgdb_package, quietly = TRUE)) {
    orgdb <- getExportedValue(orgdb_package, orgdb_package)
    if (is.null(orgdb_keytype)) orgdb_keytype <- AnnotationDbi::keytypes(orgdb)[1]
    universe_ids <- AnnotationDbi::keys(orgdb, keytype = orgdb_keytype)
  }
  if (is.null(universe_ids) && !is.null(kegg_id_map)) universe_ids <- unique(kegg_id_map$gene_id)
  tryCatch(
    {
      pathways.list <- keggList('pathway', kegg_organism)
      # Pull all genes for each pathway
      pathway.codes <- sub("path:", "", names(pathways.list))
      genes.by.pathway.loop <- sapply(
        pathway.codes,
        function(pwid) {
          pw <- keggGet(pwid)
          if (is.null(pw[[1]]$GENE)) {
            return(NA)
          }
          pw2 <- pw[[1]]$GENE[c(TRUE, FALSE)] # may need to modify this to c(FALSE, TRUE) for other organisms
          pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x) x[1]))
          return(pw2)
        }
      )
      if (!is.null(kegg_id_map)) {
        genes.by.pathway.loop <- lapply(genes.by.pathway.loop, function(ids) {
          unique(kegg_id_map$gene_id[match(ids, kegg_id_map$kegg_id, nomatch = 0L)])
        })
      }
      message("create KEGG 'gene to pathway' dataset: successfully")
    },
    error = function(cond) {
      stop("create KEGG 'gene to pathway' dataset: fail")
    }
  )

  # run kegg pathway plot fun
  message("KEGG pathways for annotated DMRs...")
  n.cores.ann <- ifelse(n.cores >= 6, 6, n.cores)
  n.cores.cntx <- ifelse(n.cores >= 18, 3, n.cores)

  mclapply(c("Genes", "Promoters", "CDS", "Introns", "fiveUTRs", "threeUTRs"), function(ann_loop) {
    mclapply(c("CG", "CHG", "CHH"), function(contx_loop) {
      for (gain_loss_loop in c("gain", "loss")) {
        tryCatch(
          {
            kegg.pathway.fun(
              treatment = comparison_name,
              gain_OR_loss = gain_loss_loop,
              context = contx_loop,
              annotation = ann_loop,
              gene2pathway = genes.by.pathway.loop,
              pathways.list = pathways.list,
              genome_ann_path = genome_ann_path,
              path_for_results = KEGG_path,
              kegg_organism = kegg_organism,
              universe_ids = universe_ids
            )
          },
          error = function(cond) {
            message(paste0("Error in 'kegg.pathway.fun': ", paste(ann_loop, contx_loop, gain_loss_loop, sep = "-")))
          }
        )
      }
    }, mc.cores = n.cores.cntx)
    tryCatch(
      {
        KEGG_one_plot(comparison_name, ann_loop, KEGG_path)
        message(paste0("\tdone: ", ann_loop))
      },
      error = function(cond) {
        message(paste0("Error in kegg pathway plot: ", ann_loop))
      }
    )
  }, mc.cores = n.cores.ann)
}

########################################################

kegg.pathway.fun = function(treatment,
                            gain_OR_loss, # "gain" OR "loss"
                            context, # "CG", "CHG", "CHH" or "all"
                            annotation, # "Genes" for example
                            gene2pathway = NULL,
                            pathways.list,
                            pValue.kegg = 0.01,
                            genome_ann_path,
                            path_for_results,
                            save.files = T,
                            kegg_organism,
                            universe_ids = NULL) {
  
  ##########################################
  start_path = getwd()
  
  #####
  DMR_file = read.csv(paste0(genome_ann_path,"/",context,"/DMRs_",annotation,"_",context,"_genom_annotations.csv"))
  DMR_file = DMR_file[DMR_file$regionType == gain_OR_loss,
                      c("gene_id","pValue")]
  #####
  if (is.null(universe_ids)) universe_ids <- unique(DMR_file$gene_id)
  organism_ids <- data.frame(gene_id = universe_ids)
  all_genes = merge.data.frame(DMR_file, organism_ids, by = "gene_id", all.y = T)
  all_genes$pValue[is.na(all_genes$pValue)] <- 0.999
  all_genes$pValue[all_genes$pValue == 0] <- 1e-300
  
  geneList = all_genes$pValue
  names(geneList) = all_genes$gene_id
  
  # Pull all genes for each pathway
  if (is.null(gene2pathway)) {
    pathways.list <- keggList('pathway', kegg_organism)
    pathway.codes <- sub("path:", "", names(pathways.list)) 
    genes.by.pathway <- sapply(pathway.codes,
                               function(pwid){
                                 pw <- keggGet(pwid)
                                 if (is.null(pw[[1]]$GENE)) return(NA)
                                 pw2 <- pw[[1]]$GENE[c(TRUE,FALSE)] # may need to modify this to c(FALSE, TRUE) for other organisms
                                 pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
                                 return(pw2)
                               }
    )
  } else {genes.by.pathway = gene2pathway}
  
  # Wilcoxon test for each pathway
  pVals.by.pathway <- t(sapply(names(genes.by.pathway),
                               function(pathway) {
                                 pathway.genes <- genes.by.pathway[[pathway]]
                                 list.genes.in.pathway <- intersect(names(geneList), pathway.genes)
                                 list.genes.not.in.pathway <- setdiff(names(geneList), list.genes.in.pathway)
                                 scores.in.pathway <- geneList[list.genes.in.pathway]
                                 scores.not.in.pathway <- geneList[list.genes.not.in.pathway]
                                 
                                 if (length(scores.in.pathway) > 0){
                                   p.value <- wilcox.test(scores.in.pathway, scores.not.in.pathway, alternative = "less")$p.value
                                   significants.in.pathway <- scores.in.pathway[scores.in.pathway < 0.05]
                                 } else{
                                   p.value <- NA
                                   significants.in.pathway <- 0
                                 }
                                 return(c(p.value = p.value,
                                          Significant = length(significants.in.pathway),
                                          Annotated = length(list.genes.in.pathway)))
                               }
  ))
  
  # Assemble output table
  outdat <- data.frame(pathway.code = rownames(pVals.by.pathway))
  outdat$pathway.name <- pathways.list[outdat$pathway.code]
  outdat$p.value <- pVals.by.pathway[,"p.value"]
  outdat$Significant <- pVals.by.pathway[,"Significant"]
  outdat$Annotated <- pVals.by.pathway[,"Annotated"]
  outdat <- outdat[order(outdat$p.value),]
  outdat$pathway.name = sub(' - [^-]+$', '', outdat$pathway.name)
  outdat$pathway.name = gsub(",",";", outdat$pathway.name)
  outdat = outdat[outdat$p.value <= pValue.kegg,] %>% na.omit()
  
  
  ### save
  new_path_1 = paste(path_for_results,context, sep = "/")
  new_path = paste(path_for_results,context,annotation, sep = "/")
  if (save.files == T) {
    if (dir.exists(new_path_1) == F) {
      dir.create(new_path_1)
    }
    if (dir.exists(new_path) == F) {
      dir.create(new_path)
    }
    write.csv(outdat, paste0(new_path,"/",annotation,".",context,".",
                             gain_OR_loss,".",treatment,".KEGG.csv"),
              quote = F, row.names = F)
  }
  
  #message(paste(treatment,annotation,context,gain_OR_loss,">> done", sep = " - "))
  setwd(start_path)
}

########################################################

KEGG_one_plot <- function(treatment,
                          annotation,
                          path_for_results,
                          breaks_gain = NULL, # c(2,4,6),
                          limits_gain = NULL, # c(0.5,6.5)
                          breaks_loss = NULL,
                          limits_loss = NULL
) {
  
  cg_gain = read.csv(paste0(path_for_results,"/CG/",annotation,"/",annotation,".CG.gain.",treatment,".KEGG.csv"))
  cg_loss = read.csv(paste0(path_for_results,"/CG/",annotation,"/",annotation,".CG.loss.",treatment,".KEGG.csv"))
  chg_gain = read.csv(paste0(path_for_results,"/CHG/",annotation,"/",annotation,".CHG.gain.",treatment,".KEGG.csv"))
  chg_loss = read.csv(paste0(path_for_results,"/CHG/",annotation,"/",annotation,".CHG.loss.",treatment,".KEGG.csv"))
  chh_gain = read.csv(paste0(path_for_results,"/CHH/",annotation,"/",annotation,".CHH.gain.",treatment,".KEGG.csv"))
  chh_loss = read.csv(paste0(path_for_results,"/CHH/",annotation,"/",annotation,".CHH.loss.",treatment,".KEGG.csv"))
  
  if (length(cg_gain$pathway.code) != 0) {cg_gain$type = "CG"}
  if (length(cg_loss$pathway.code) != 0) {cg_loss$type = "CG"}
  if (length(chg_gain$pathway.code) != 0) {chg_gain$type = "CHG"}
  if (length(chg_loss$pathway.code) != 0) {chg_loss$type = "CHG"}
  if (length(chh_gain$pathway.code) != 0) {chh_gain$type = "CHH"}
  if (length(chh_loss$pathway.code) != 0) {chh_loss$type = "CHH"}
  
  gain_bind = rbind(cg_gain, chg_gain, chh_gain)
  loss_bind = rbind(cg_loss, chg_loss, chh_loss)
  
  gain_col = "#cf534c"
  loss_col = "#6397eb"
  
  bubble_gain = gain_bind %>% 
    ggplot(aes(Significant, reorder(pathway.name,Significant), color = p.value, size = Annotated)) + 
    scale_color_gradient("p.value", low = gain_col, high = "black") + #theme_classic() +
    labs(#title = paste0(type_name," - ",gain_loss,"regulated transcripts"),
      x = "Significant", y = "") + theme_bw() + 
    theme(
      #plot.title=element_text(hjust=0.5),
      legend.key.size = unit(0.25, 'cm'),
      legend.title = element_text(size=9.5),
      legend.position = "right",
      text = element_text(family = "serif")) + 
    geom_point() +
    facet_grid(rows = vars(type), scales = "free_y", space = "free_y") + 
    guides(color = guide_colorbar(order = 1, barheight = 4))
  
  if (is.null(breaks_gain) == F & is.null(limits_gain) == F) {
    bubble_gain = bubble_gain + scale_x_continuous(breaks=breaks_gain,
                                                   limits=limits_gain)
  }
  
  
  
  bubble_loss = loss_bind %>% 
    ggplot(aes(Significant, reorder(pathway.name,Significant), color = p.value, size = Annotated)) + 
    scale_color_gradient("p.value", low = loss_col, high = "black") + #theme_classic() +
    labs(#title = paste0(type_name," - ",gain_loss,"regulated transcripts"),
      x = "Significant", y = "") + theme_bw() + 
    theme(
      #plot.title=element_text(hjust=0.5),
      legend.key.size = unit(0.25, 'cm'),
      legend.title = element_text(size=9.5),
      legend.position = "right",
      text = element_text(family = "serif")) + 
    geom_point() +
    facet_grid(rows = vars(type), scales = "free_y", space = "free_y") + 
    guides(color = guide_colorbar(order = 1, barheight = 4))
  
  if (is.null(breaks_loss) == F & is.null(limits_loss) == F) {
    bubble_loss = bubble_loss + scale_x_continuous(breaks=breaks_loss,
                                                   limits=limits_loss)
  }
  
  Height = max(c(nrow(gain_bind),nrow(loss_bind)))/6.25
  if (Height < 3) {Height = 3}
  
  img_device(paste0(path_for_results,"/",treatment,".",annotation,".KEGG"),
      width = 9.90, h = Height)
  multiplot(bubble_gain, bubble_loss, cols=2)
  dev.off()
  
  #  print(paste(treatment,annotation, sep = "_"))
}
