run_GO <- function(comparison_name, genome_ann_path, GO_path, n.cores,
                   reference_bundle = NULL) {

  orgdb_package <- bundle_get(reference_bundle, 'functional.orgdb_package')
  orgdb_keytype <- bundle_get(reference_bundle, 'functional.orgdb_keytype')
  gene_to_go_path <- bundle_get(reference_bundle, 'functional.gene_to_go')
  gene2go <- NULL
  if (!is.null(gene_to_go_path)) {
    gene_to_go <- data.table::fread(gene_to_go_path, data.table = FALSE)
    if (!all(c('gene_id', 'go_id') %in% names(gene_to_go))) {
      stop('gene_to_go must contain gene_id and go_id columns.')
    }
    gene2go <- split(as.character(gene_to_go$go_id), as.character(gene_to_go$gene_id))
  }
  if (is.null(orgdb_package) && is.null(gene2go)) {
    stop('GO analysis requires functional.orgdb_package or functional.gene_to_go.')
  }
  if (!is.null(orgdb_package) && !requireNamespace(orgdb_package, quietly = TRUE)) {
    stop('GO annotation package is not installed: ', orgdb_package)
  }

  n.cores.ann = ifelse(n.cores >= 6, 6, n.cores)
  n.cores.cntx = ifelse(n.cores >= 18, 3, n.cores)

  mclapply(c("Genes", "Promoters", "CDS", "Introns", "fiveUTRs", "threeUTRs"), function(ann_loop) {
    mclapply(c("CG", "CHG", "CHH"), function(contx_loop) {
      for (gain_loss_loop in c("gain", "loss")) {
        for (GO_type_loop in c("BP", "MF", "CC")) {
          tryCatch(
            {
              suppressMessages(top.GO.fun(
                treatment = comparison_name,
                gain_OR_loss = gain_loss_loop,
                GO_Ontology_type = GO_type_loop,
                context = contx_loop,
                annotation = ann_loop,
                genome_ann_path = genome_ann_path,
                path_for_results = GO_path,
                orgdb_package = orgdb_package,
                orgdb_keytype = orgdb_keytype,
                gene2go = gene2go
              ))
            },
            error = function(cond) {
              message(paste0("Error in 'top.GO.fun': ", paste(ann_loop, contx_loop, gain_loss_loop, GO_type_loop, sep = "-")))
            }
          )
        }
      }
      tryCatch(
        {
          GO_one_plot(comparison_name, contx_loop, ann_loop, GO_path)
          # message(paste0("\tdone: ",ann_loop," in ",contx_loop," context\t"))
        },
        error = function(cond) {
          message(paste0("Error in GO plot: ", paste(ann_loop, contx_loop, sep = "-")))
        }
      )
    }, mc.cores = n.cores.cntx)
    message(paste0("\tdone: ", ann_loop))
  }, mc.cores = n.cores.ann)
}

########################################################

top.GO.fun = function(treatment,
                      gain_OR_loss, # "gain" OR "loss"
                      GO_Ontology_type, # BP, CC, MF
                      context, # "CG", "CHG", "CHH" or "all"
                      annotation, # "Genes" for example
                      pcutoff = 0.01,
                      GO_algorithm = "weight01", # classic, weight, weight01
                      genome_ann_path,
                      path_for_results,
                      n.nodes = NULL,
                      get_GO_term_genes = c(F, "GO_term"),
                      save.files = T,
                      orgdb_package = NULL,
                      orgdb_keytype = NULL,
                      gene2go = NULL) {
  
  ##########################################
  start_path = getwd()
  
  #####
  DMR_file = read.csv(paste0(genome_ann_path,"/",context,"/DMRs_",annotation,"_",context,"_genom_annotations.csv"))
  DMR_file = DMR_file[DMR_file$regionType == gain_OR_loss,
                      c("gene_id","pValue")]
  
  
  if (!is.null(gene2go)) {
    organism_ids <- data.frame(gene_id = names(gene2go))
  } else {
    orgdb <- getExportedValue(orgdb_package, orgdb_package)
    if (is.null(orgdb_keytype)) orgdb_keytype <- AnnotationDbi::keytypes(orgdb)[1]
    organism_ids <- data.frame(
      gene_id = AnnotationDbi::keys(orgdb, keytype = orgdb_keytype)
    )
  }
  
  all_genes = merge.data.frame(DMR_file, organism_ids, by = "gene_id", all.y = T)
  all_genes$pValue[is.na(all_genes$pValue)] <- 0.999
  all_genes$pValue[all_genes$pValue == 0] <- 1e-300
  
  geneList = ifelse(all_genes$pValue < pcutoff, 1, 0)
  names(geneList) = all_genes$gene_id
  
  ############################################################
  #### type Ontology  ####
  if (!is.null(gene2go)) {
    myGOdata <- new('topGOdata',
                    ontology = GO_Ontology_type,
                    allGenes = geneList,
                    geneSelectionFun = function(x)(x == 1),
                    annot = annFUN.gene2GO,
                    gene2GO = gene2go)
  } else {
    myGOdata <- new('topGOdata',
                    ontology = GO_Ontology_type,
                    allGenes = geneList,
                    geneSelectionFun = function(x)(x == 1),
                    annot = annFUN.org,
                    mapping = orgdb_package,
                    ID = orgdb_keytype)
  }
  
  sg <- sigGenes(myGOdata)
  str(sg)
  numSigGenes(myGOdata)
  
  resultFisher <- runTest(myGOdata, algorithm = GO_algorithm, statistic = "fisher")
  
  
  allRes <- GenTable(myGOdata, Fisher = resultFisher,
                     orderBy= "Fisher", ranksOf= "Fisher", topNodes=length(resultFisher@score))
  allRes$Fisher = as.numeric(allRes$Fisher)
  allRes$Term = gsub(",", ";", allRes$Term)
  
  allRes_view = allRes[allRes$Fisher <= pcutoff,]
  allRes_view = allRes_view[!is.na(allRes_view$GO.ID), ]
  
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
    write.csv(allRes_view, paste0(new_path,"/",GO_Ontology_type ,"." ,annotation,
                                  ".",context,".",gain_OR_loss,".",treatment,".topGO.csv"),
              quote = F, row.names = F)
  }
  
  
  if (is.null(n.nodes)) {
    if (nrow(allRes_view) > 10) {
      n.nodes = 10
    } else {
      n.nodes = nrow(allRes_view)
    }
  }
  
  setwd(new_path)
  if (nrow(allRes_view) > 0) {
    printGraph(myGOdata, resultFisher, 
               firstSigNodes = n.nodes, fn.prefix= paste(GO_Ontology_type,annotation,context,gain_OR_loss,treatment, sep = "_"),
               useInfo="all", pdfSW= T)
  }
  setwd(start_path)
}

########################################################

GO_one_plot <- function(treatment,
                        context,
                        annotation,
                        path_for__results,
                        breaks_gain = NULL, # c(2,4,6),
                        limits_gain = NULL, # c(0.5,6.5)
                        breaks_loss = NULL,
                        limits_loss = NULL
) {
  
  new_path = paste(path_for__results,context,annotation, sep = "/")
  
  bp_gain = read.csv(paste0(new_path,"/BP.",annotation,".",context,".gain.",treatment,".topGO.csv"))
  bp_loss = read.csv(paste0(new_path,"/BP.",annotation,".",context,".loss.",treatment,".topGO.csv"))
  cc_gain = read.csv(paste0(new_path,"/CC.",annotation,".",context,".gain.",treatment,".topGO.csv"))
  cc_loss = read.csv(paste0(new_path,"/CC.",annotation,".",context,".loss.",treatment,".topGO.csv"))
  mf_gain = read.csv(paste0(new_path,"/MF.",annotation,".",context,".gain.",treatment,".topGO.csv"))
  mf_loss = read.csv(paste0(new_path,"/MF.",annotation,".",context,".loss.",treatment,".topGO.csv"))
  
  if (length(bp_gain$GO.ID) != 0) {bp_gain$type = "BP"}
  if (length(bp_loss$GO.ID) != 0) {bp_loss$type = "BP"}
  if (length(cc_gain$GO.ID) != 0) {cc_gain$type = "CC"}
  if (length(cc_loss$GO.ID) != 0) {cc_loss$type = "CC"}
  if (length(mf_gain$GO.ID) != 0) {mf_gain$type = "MF"}
  if (length(mf_loss$GO.ID) != 0) {mf_loss$type = "MF"}
  
  gain_bind = rbind(bp_gain, cc_gain, mf_gain)
  loss_bind = rbind(bp_loss, cc_loss, mf_loss)
  
  gain_bind = gain_bind[!grepl("cellular_component|biological_process|molecular_function", gain_bind$Term),]
  loss_bind = loss_bind[!grepl("cellular_component|biological_process|molecular_function", loss_bind$Term),]
  
  gain_col = "#cf534c"
  loss_col = "#6397eb"
  
  bubble_gain = gain_bind %>% 
    ggplot(aes(Significant, reorder(Term,-Fisher), size = Annotated, color = Fisher)) + 
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
    ggplot(aes(Significant, reorder(Term,-Fisher), size = Annotated, color = Fisher)) + 
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
  #print(paste(treatment,context,annotation, sep = "_"))
  
  img_device(paste0(path_for__results,"/",treatment,"_",context,"_",annotation,"_GO"),
      width = 9.90, h = Height)
  multiplot(bubble_gain, bubble_loss, cols=2)
  dev.off()
}
