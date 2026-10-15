Revigo_plots <- function(GO_df_up, GO_df_down,
                         GO_type, #direction,
                         show_labels=TRUE,
                         reference_bundle,
                         semdata=NULL,
                         data_only = FALSE,
                         prepared_data = NULL,
                         specific_terms_lables_up=NULL, # specific by term column
                         specific_terms_lables_down=NULL, # specific by term column
                         specific_parents_lables_up=NULL, # specific by parents column
                         specific_parents_lables_down=NULL, # specific by parents column
                         specific_parents_lables_show=NULL # if null, will take the default terms label
                         ) {
  
  orgdb_package <- bundle_get(reference_bundle, "functional.orgdb_package")
  orgdb_keytype <- bundle_get(reference_bundle, "functional.orgdb_keytype")
  
  if (is.null(semdata)) {
    semdata <- suppressMessages(
      GOSemSim::godata(
        annoDb = orgdb_package,
        keytype = orgdb_keytype,
        ont = GO_type
      )
    )
  }

  # prepare revigo data frame
  revigo_df_fun <- function(GO_df) {
    
    # if (!is.numeric(GO_df$pValue)) {
    #   GO_df$pValue = gsub("< 1e-30", 1e-30, GO_df$pValue) %>% as.numeric()
    # }
    empty_result <- list(
      plot_df = data.frame(
        term_ID = character(),
        description = character(),
        plot_X = numeric(),
        plot_Y = numeric(),
        value = numeric(),
        parent = character(),
        parentTerm = character()
      ),
      reducedTerms = data.frame()
    )

    GO_df$pValue <- suppressWarnings(as.numeric(
      sub("^\\s*<\\s*", "", as.character(GO_df$pValue))
    ))

    GO_df <- GO_df[
      !is.na(GO_df$GO.ID) & nzchar(GO_df$GO.ID) &
        is.finite(GO_df$pValue) &
        GO_df$pValue >= 0 & GO_df$pValue <= 1,
      , drop = FALSE
    ]

    GO_df <- GO_df[order(GO_df$pValue), , drop = FALSE]
    GO_df <- GO_df[!duplicated(GO_df$GO.ID), , drop = FALSE]
    GO_df$pValue <- pmax(GO_df$pValue, .Machine$double.xmin)

    if (nrow(GO_df) < 3L) {
      message(time_msg(), "REVIGO ", GO_type,
              ": direction skipped (fewer than three usable GO terms)")
      return(empty_result)
    }
    
    ### rrvgo pipeline

    simMatrix <- calculateSimMatrix(
      GO_df$GO.ID,
      orgdb = orgdb_package,
      keytype = orgdb_keytype,
      ont = GO_type,
      semdata = semdata,
      method = "Rel"
    )

    if (!is.matrix(simMatrix) || nrow(simMatrix) < 3L) {
      message(time_msg(), "REVIGO ", GO_type,
              ": direction skipped (fewer than three terms retained)")
      return(empty_result)
    }

    scores <- setNames(-log10(GO_df$pValue), GO_df$GO.ID)

    reducedTerms <- suppressMessages(
      reduceSimMatrix(
        simMatrix,
        scores,
        threshold = 0.7,
        orgdb = orgdb_package,
        keytype = orgdb_keytype
      )
    )
    
    # tryCatch({
    #   x <<- scatterPlot(simMatrix, reducedTerms, algorithm = "umap")$data
    # }, error = function(e){
    #   message(time_msg(), "using 'PCA' instead of 'UMAP'")
    #   x <<- scatterPlot(simMatrix, reducedTerms, algorithm = "pca")$data
    # })
    x <- tryCatch(
      scatterPlot(simMatrix, reducedTerms, algorithm = "umap")$data,
      error = function(e) {
        message(time_msg(), "UMAP failed: ", conditionMessage(e),
                "; using the PCA fallback")
        scatterPlot(simMatrix, reducedTerms, algorithm = "pca")$data
      }
    )

    xx = data.frame(term_ID = row.names(x),
                           description = as.character(x$term),
                           plot_X = as.numeric(x$V1),
                           plot_Y = as.numeric(x$V2),
                           value = as.numeric(x$score),
                           parent = as.character(x$parent),
                           parentTerm = as.character(x$parentTerm))
    return(list(plot_df = xx,
                reducedTerms = reducedTerms))
  }

  # revigo_up_0 = revigo_df_fun(GO_df_up)
  # revigo_up = revigo_up_0$plot_df
  # revigo_up$direction <- rep("gain", nrow(revigo_up))
  # 
  # revigo_down_0 = revigo_df_fun(GO_df_down)
  # revigo_down = revigo_down_0$plot_df
  # revigo_down$direction <- rep("loss", nrow(revigo_down))
  if (is.null(prepared_data)) {
    revigo_up_0 <- revigo_df_fun(GO_df_up)
    revigo_down_0 <- revigo_df_fun(GO_df_down)
  } else {
    revigo_up_0 <- prepared_data$gain
    revigo_down_0 <- prepared_data$loss
  }

  if (data_only) {
    return(list(gain = revigo_up_0, loss = revigo_down_0))
  }

  revigo_up <- revigo_up_0$plot_df
  revigo_up$direction <- rep("gain", nrow(revigo_up))

  revigo_down <- revigo_down_0$plot_df
  revigo_down$direction <- rep("loss", nrow(revigo_down))
  
  ## bind 'up and 'down' data frames
  revigo_df = rbind(revigo_up,revigo_down) %>% arrange(parentTerm)# %>% arrange(desc(direction))
  if (!nrow(revigo_df)) return(NULL)
  
  ## set colors for different parent terms
  # n.unique_parents = length(unique(revigo_df$parentTerm))
  # if (n.unique_parents > 17) {
  #   color_palette = rainbow(n.unique_parents)
  # } else {
  #   color_palette <- suppressWarnings(brewer.pal(n = n.unique_parents, "Set1"))
  #   if (length(unique(revigo_df$parentTerm)) > 9) {
  #     color_palette = suppressWarnings(c(color_palette, brewer.pal(n = n.unique_parents-9, "Set2")))
  #   }
  # }
  parents <- unique(revigo_df$parentTerm)
  color_palette <- setNames(
    grDevices::hcl.colors(length(parents), "Dark 3"),
    parents
  )
  
  ## parents and color index
  # parent_col = data.frame(parentTerm = unique(revigo_df$parentTerm), hex_col = color_palette)
  parent_col <- data.frame(parentTerm = parents, hex_col = unname(color_palette))
  revigo_df = merge(revigo_df, parent_col, by = "parentTerm")
  revigo_df_up <- revigo_df[revigo_df$direction == "gain", ]
  revigo_df_down <- revigo_df[revigo_df$direction == "loss", ]
  
  
  ## max n.char for margin argument
  # terms_character = paste(specific_terms_lables_up,
  #                         specific_terms_lables_down,
  #                         specific_parents_lables_up,
  #                         specific_parents_lables_down, sep = "|")
  # if (length(terms_character) != 0) {
  #   lables_nchar = max(nchar(strsplit(terms_character, split = "\\|")[[1]]))
  #   lables_nchar_index = lables_nchar*0.15
  # } else {
  #   lables_nchar_index = ifelse(show_labels, 6, 0.1)
  # }
  label_patterns <- c(
    specific_terms_lables_up,
    specific_terms_lables_down,
    specific_parents_lables_up,
    specific_parents_lables_down
  )

  lables_nchar_index <- if (!show_labels) {
    0.1
  } else if (!length(label_patterns)) {
    6
  } else {
    max(nchar(label_patterns)) * 0.15
  }
  
  ##############################
  
  scatterPlot_function <- function(df, specific_terms_lables, specific_parents_lables, margin_index = lables_nchar_index) {
    if (!nrow(df)) return(NULL)

    one.x_range = (max(df$plot_X) - min(df$plot_X))
    one.y_range = (max(df$plot_Y) - min(df$plot_Y))
    
    #### labels - filter terms - 4 options
    # by specific terms
    if (!is.null(specific_terms_lables)) {
      filtered_data = df[grep(specific_terms_lables, df$description),]%>% # keep specific terms
        dplyr::distinct(description, .keep_all = TRUE) # remove duplicate terms
      filtered_data$parentTerm = filtered_data$description
      
    } else if (!is.null(specific_parents_lables)) {
      filtered_data = df[grep(specific_parents_lables, df$parentTerm),]%>% # keep specific terms
        dplyr::distinct(parentTerm, .keep_all = TRUE) # remove duplicate terms
      
      # change to costum label
      if (!is.null(specific_parents_lables_show)) {
        # Replace values in filtered_data$parentTerm
        lookup <- setNames(specific_parents_lables_show$b, specific_parents_lables_show$a)
        filtered_data$parentTerm <- ifelse(filtered_data$parentTerm %in% names(lookup),
                                           lookup[filtered_data$parentTerm],
                                           filtered_data$parentTerm)
      }
    } else {
      
      
      n.n = 8
      n.unique = length(unique(df$parentTerm))
      # less than n.n
      if (n.unique > n.n) {
        
        filtered_data <- df %>%
          dplyr::group_by(parentTerm) %>%
          dplyr::filter(dplyr::n() >= 2) %>% # filter term that appears less than 3 times
          dplyr::ungroup() %>%
          dplyr::distinct(parentTerm, .keep_all = TRUE) %>% # remove duplicate terms
          as.data.frame()
        
        # more than n.n
        if (nrow(filtered_data) < n.n) {
          filtered_data = df[!duplicated(df$parentTerm),]
        } else {
          filtered_data <- df %>%
            dplyr::group_by(parentTerm) %>%
            dplyr::summarize(count = dplyr::n()) %>% # keep top n.n terms
            dplyr::arrange(desc(count)) %>%
            dplyr::slice_max(order_by = count, n = n.n) %>%
            dplyr::slice_head(n = n.n) %>% 
            dplyr::select(parentTerm) %>%
            dplyr::inner_join(df, by = "parentTerm") %>%
            dplyr::ungroup()  %>%
            dplyr::distinct(parentTerm, .keep_all = TRUE) %>%
            as.data.frame() 
        }
        
      } else {
        filtered_data = df[!duplicated(df$parentTerm),] # keep all terms
      }
    }
    
    
    scatterPlot_vec = df %>% ggplot2::ggplot(ggplot2::aes(x = plot_X, y = plot_Y, color = parentTerm)) + 
      ggplot2::geom_point(ggplot2::aes(size = value), alpha = 0.5, show.legend = FALSE) + 
      
      # ggplot2::scale_color_manual(values = unique(df$hex_col), guide = "none") + 
      ggplot2::scale_color_manual(values = color_palette, guide = "none") +
      ggplot2::geom_point(ggplot2::aes(plot_X, plot_Y, size = value), 
                          shape = 21, fill = "transparent", colour = I(ggplot2::alpha("black", 0.6)),
                          show.legend = FALSE) + 
      
      ggplot2::theme_classic() + 
      ggplot2::geom_hline(yintercept=0, colour="grey50", linetype = "dashed") + 
      ggplot2::geom_vline(xintercept = 0, colour="grey50", linetype = "dashed") +
      
      ggplot2::theme(panel.border = element_rect(colour = "black", fill=NA, size=0.75),
                     axis.text.y = element_text(size = 9),
                     axis.text.x = element_text(size = 9),
                     axis.line = element_blank(),
                     plot.margin =  unit(c(0.1, margin_index, 0.1, 0.1), "cm")) +
      
      # ggrepel::geom_text_repel(ggplot2::aes(label = parentTerm), 
      #                          data = filtered_data,
      #                          direction = "y",
      #                          hjust = 0,
      #                          #alpha = 0.65,
      #                          #min.segment.length = 0,
      #                          #family = 'Times',
      #                          fontface = 'bold',
      #                          #box.padding = 1,
      #                          #box.padding = 0.35, 
      #                          #point.padding = 0.5,
      #                          #point.padding = NA,
      #                          #label.padding = .2,
      #                          size = 3,
      #                          #max.time = 10,
      #                          #max.iter = 1000000,
      #                          #force = 5,
      #                          colour="black",
      #                          xlim = c(max(df$plot_X)+one.x_range/7.5, Inf)) + 
      (if (show_labels) {
        ggrepel::geom_text_repel(
          ggplot2::aes(label = parentTerm),
          data = filtered_data,
          direction = "y",
          hjust = 0,
          fontface = "bold",
          size = 3,
          colour = "black",
          xlim = c(max(df$plot_X) + one.x_range / 7.5, Inf)
        )
      } else NULL) +
      
      ggplot2::labs (y = "semantic space y", x = "semantic space x") +#, size = "-log10(p)") + 
      #ggplot2::guides(color = guide_legend(override.aes = list(size = 4, shape = 21, colour = "black", fill = unique(df$hex_col)))) +
      #ggplot2::theme(legend.key = ggplot2::element_blank(), legend.box = "horizontal") +
      ggplot2::coord_cartesian(clip = "off",
                               xlim = c(min(df$plot_X)-one.x_range/10,max(df$plot_X)+one.x_range/10),
                               ylim = c(min(df$plot_Y)-one.y_range/10,max(df$plot_Y)+one.y_range/10))
    
    return(scatterPlot_vec)
  }
  scatterPlot_up = scatterPlot_function(revigo_df_up, specific_terms_lables_up, specific_parents_lables_up)
  scatterPlot_down = scatterPlot_function(revigo_df_down, specific_terms_lables_down, specific_parents_lables_down)
  
 
  return(list(scatterPlot_up = scatterPlot_up,
              scatterPlot_down = scatterPlot_down,
              scatterPlot_legend = parent_col,
              treemapPlot_up = revigo_up_0$reducedTerms,
              treemapPlot_down = revigo_down_0$reducedTerms,
              revigo_df = revigo_df))
}

#############################################################################

scatter_plot_run_function <- function(x,x.file, x.file.suffix) {
  save_plot <- function(name, w, h, code) {
    img_device(
      file.path(x.file, paste0(x.file.suffix, "_", name)),
      w = w, h = h
    )
    plot_device <- dev.cur()
    on.exit(dev.off(plot_device), add = TRUE)
    force(code)
  }

  ## save scatter plot
  # svg(paste0(x.file,"up_ScatterPlot",x.file.suffix,".svg"), width = 5, height = 2.38, family = "serif")
  # print(x$scatterPlot_up)
  # dev.off()

  # svg(paste0(x.file,"down_ScatterPlot",x.file.suffix,".svg"), width = 5, height = 2.38, family = "serif")
  # print(x$scatterPlot_down)
  # dev.off()
  if (!is.null(x$scatterPlot_up)) {
    save_plot("gain_ScatterPlot", 5, 2.38, print(x$scatterPlot_up))
  }

  if (!is.null(x$scatterPlot_down)) {
    save_plot("loss_ScatterPlot", 5, 2.38, print(x$scatterPlot_down))
  }

  ## save scatter plot legend
  # svg(paste0(x.file,"ScatterPlot_legend",x.file.suffix,".svg"), width = 4.32, height = 10, family = "serif")
  save_plot("ScatterPlot_legend", 4.32, ifelse(nrow(x$scatterPlot_legend) > 40, 10, 6), {
  par(mar = c(0.1, 0.1, 1.1, 0.1))
  plot.new()
  legend("topleft",
         legend = x$scatterPlot_legend$parentTerm,
         pch = 21,
         pt.cex = 1.5,
         pt.bg = paste0(x$scatterPlot_legend$hex_col,75),
         col = "black", 
         cex = 0.75, 
         bty = "n",
         x.intersp = 0.8,
         y.intersp = 1)
  title(main = "Parent Term", cex.main = 1, adj = 0.01)
  })

  ## save tree map
  # svg(paste0(x.file,"up_treeMap",x.file.suffix,".svg"), width = 5, height = 5, family = "serif")
  # treemapPlot(x$treemapPlot_up)
  # dev.off()

  # svg(paste0(x.file,"down_treeMap",x.file.suffix,".svg"), width = 5, height = 5, family = "serif")
  # treemapPlot(x$treemapPlot_down)
  # dev.off()
  if (nrow(x$treemapPlot_up)) {
    save_plot("gain_treeMap", 5, 5, rrvgo::treemapPlot(x$treemapPlot_up))
  }

  if (nrow(x$treemapPlot_down)) {
    save_plot("loss_treeMap", 5, 5, rrvgo::treemapPlot(x$treemapPlot_down))
  }

  # save scatter plot data frame
  write.csv(x$revigo_df, file.path(x.file, paste0(x.file.suffix, "_ScatterData.csv")), row.names = F) 
}

#############################################################################

run_REVIGO <- function(go_results_path, revigo_path, comparison_name, reference_bundle, n.cores = 1L) {
  statuses <- character()

  orgdb_package <- bundle_get(reference_bundle, "functional.orgdb_package")
  orgdb_keytype <- bundle_get(reference_bundle, "functional.orgdb_keytype")

  if (is.null(orgdb_package) || is.null(orgdb_keytype)) {
    message(time_msg(),
            "REVIGO: skipped (OrgDb package/keytype absent from bundle)")
    return(list(completed = 0L, skipped = 1L, failed = 0L))
  }

  ontologies <- c("BP", "MF", "CC")
  semdata_by_ontology <- setNames(
    lapply(ontologies, function(ontology) {
      suppressMessages(
        GOSemSim::godata(
          annoDb = orgdb_package,
          keytype = orgdb_keytype,
          ont = ontology
        )
      )
    }),
    ontologies
  )

  task_table <- expand.grid(
    type = ontologies,
    ann = c("Genes", "Promoters", "CDS", "Introns"),
    context = c("CG", "CHG", "CHH", "all"),
    stringsAsFactors = FALSE
  )
  task_table$id <- seq_len(nrow(task_table))
  jobs <- lapply(seq_len(nrow(task_table)), function(i) {
    as.list(task_table[i, , drop = FALSE])
  })

  revigo_settings <- list(
    go_path = normalizePath(go_results_path, mustWork = TRUE),
    comparison = comparison_name,
    bundle = reference_bundle,
    semdata = semdata_by_ontology
  )

  prepare_job <- function(job, settings = revigo_settings) {
    tryCatch({
      contexts <- if (job$context == "all") {
        c("CG", "CHG", "CHH")
      } else {
        job$context
      }

      read_go <- function(direction) {
        tables <- lapply(contexts, function(context) {
          filename <- paste(
            job$type, job$ann, context, direction,
            settings$comparison, "topGO.csv",
            sep = "."
          )
          read.csv(file.path(
            settings$go_path, context, job$ann, filename
          ))
        })
        tab <- do.call(rbind, tables)
        names(tab)[names(tab) == "Fisher"] <- "pValue"
        tab
      }

      set.seed(job$id)
      data <- Revigo_plots(
        GO_df_up = read_go("gain"),
        GO_df_down = read_go("loss"),
        GO_type = job$type,
        reference_bundle = settings$bundle,
        semdata = settings$semdata[[job$type]],
        data_only = TRUE
      )

      list(data = data, error = NULL)
    }, error = function(e) {
      list(data = NULL, error = conditionMessage(e))
    })
  }

  # Avoid exporting the runner's entire environment.
  environment(prepare_job) <- .GlobalEnv

  workers <- min(length(jobs), max(1L, as.integer(n.cores)))
  message(time_msg(), "REVIGO: computing ", length(jobs),
          " groups using ", workers, " worker(s)")
  cat("computing ", length(jobs),
          " groups using ", workers, " worker(s)\n")

  if (workers > 1L) {
    cl <- parallel::makePSOCKcluster(workers)
    on.exit({
      if (!is.null(cl)) parallel::stopCluster(cl)
    }, add = TRUE)

    parallel::clusterCall(cl, .libPaths, .libPaths())
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages({
        library(dplyr)
        library(ggplot2)
        library(rrvgo)
        library(grid)
      })
      time_msg <- function(suffix = "\t") {
        paste0(format(Sys.time(), "[%H:%M]"), suffix)
      }
      NULL
    })

    parallel::clusterExport(
      cl,
      c("Revigo_plots", "bundle_get", "prepare_job", "revigo_settings"),
      envir = environment()
    )

    results <- parallel::parLapplyLB(cl, jobs, prepare_job)

    parallel::stopCluster(cl)
    cl <- NULL
  } else {
    results <- lapply(jobs, prepare_job, settings = revigo_settings)
  }

  # message(time_msg(), "REVIGO: saving plots sequentially...")
  cat("saving plots sequentially...\n")

  for (i in seq_along(jobs)) {
    job <- jobs[[i]]
    item <- results[[i]]
    label <- paste(job$ann, job$context, job$type, sep = "-")

    status <- tryCatch({
      if (!is.null(item$error)) stop(item$error, call. = FALSE)

      plots <- Revigo_plots(
        GO_df_up = NULL,
        GO_df_down = NULL,
        GO_type = job$type,
        reference_bundle = reference_bundle,
        semdata = semdata_by_ontology[[job$type]],
        prepared_data = item$data
      )

      if (is.null(plots)) {
        message(time_msg(), "\tskipped: ", label)
        "skipped"
      } else {
        folder <- file.path(
          revigo_path, job$context, job$ann, job$type
        )
        dir.create(folder, recursive = TRUE, showWarnings = FALSE)
        scatter_plot_run_function(plots, folder, comparison_name)
        message(time_msg(), "\tcompleted: ", label)
        "completed"
      }
    }, error = function(e) {
      message(time_msg(), "\tfailed: ", label,
              ": ", conditionMessage(e))
      "failed"
    })

    statuses <- c(statuses, status)
  }

  list(
    completed = sum(statuses == "completed"),
    skipped = sum(statuses == "skipped"),
    failed = sum(statuses == "failed")
  )
}

  # for (context.l in c("CG","CHG","CHH","all")) {
#     path2results_meth <- file.path(revigo_path, context.l)
#     dir.create(path2results_meth, recursive = TRUE, showWarnings = FALSE)
#   
#     for (ann.l in c("Genes","Promoters", "CDS", "Introns")) {
#       path2results_meth_2 <- file.path(path2results_meth, ann.l)
#       dir.create(path2results_meth_2, recursive = TRUE, showWarnings = FALSE)
#       
#       for (type.l in c("BP","MF","CC")) {
#         path2results_meth_type <- file.path(path2results_meth_2, type.l)
#         dir.create(path2results_meth_type, recursive = TRUE, showWarnings = FALSE)
#         
#         status <- tryCatch({
#           if (context.l != "all") {
#             GO_file_up = read.csv(paste0(go_results_path,"/",context.l,"/",ann.l,"/",type.l,".",ann.l,".",context.l,".gain.",comparison_name,".topGO.csv"))
#             GO_file_down = read.csv(paste0(go_results_path,"/",context.l,"/",ann.l,"/",type.l,".",ann.l,".",context.l,".loss.", comparison_name, ".topGO.csv"))
#           } else {
#             GO_file_up = rbind(read.csv(paste0(go_results_path,"/CG/",ann.l,"/",type.l,".",ann.l,".CG.gain.",comparison_name,".topGO.csv")),
#                                read.csv(paste0(go_results_path,"/CHG/",ann.l,"/",type.l,".",ann.l,".CHG.gain.",comparison_name,".topGO.csv")),
#                                read.csv(paste0(go_results_path,"/CHH/",ann.l,"/",type.l,".",ann.l,".CHH.gain.",comparison_name,".topGO.csv")))
#             GO_file_down = rbind(read.csv(paste0(go_results_path,"/CG/",ann.l,"/",type.l,".",ann.l,".CG.loss.",comparison_name,".topGO.csv")),
#                                  read.csv(paste0(go_results_path,"/CHG/",ann.l,"/",type.l,".",ann.l,".CHG.loss.",comparison_name,".topGO.csv")),
#                                  read.csv(paste0(go_results_path,"/CHH/",ann.l,"/",type.l,".",ann.l,".CHH.loss.",comparison_name,".topGO.csv")))
#           }
#           names(GO_file_up) = gsub("Fisher", "pValue", names(GO_file_up))
#           names(GO_file_down) = gsub("Fisher", "pValue", names(GO_file_down))
#           
#           #####################################
#           
#           scatterPlot_fun = Revigo_plots(GO_df_up = GO_file_up,
#                                          GO_df_down = GO_file_down,
#                                          GO_type = type.l,
#                                          reference_bundle = reference_bundle,
#                                          semdata = semdata_by_ontology[[type.l]])
#           
#           ## save scatter plot
#           # scatter_plot_run_function(scatterPlot_fun, path2results_meth_type, comparison_name)
#           if (is.null(scatterPlot_fun)) {
#             message(time_msg(), "REVIGO skipped: ",
#                     paste(ann.l, context.l, type.l, sep = "-"))
#             "skipped"
#           } else {
#             scatter_plot_run_function(
#               scatterPlot_fun, path2results_meth_type, comparison_name
#             )
#             "completed"
#           }
#           
#         # }, error = function(e){message(time_msg(), paste0("\n**\tfail: ","_",ann.l,"_",context.l,"_",type.l,"\t**\n"))
#         }, error = function(e) {
#           message(time_msg(), "REVIGO failed: ",
#                   paste(ann.l, context.l, type.l, sep = "-"),
#                   ": ", conditionMessage(e))
#           "failed"
#         })
# 
#         statuses <- c(statuses, status)
#       }
#     }
  # }
  # list(
#     completed = sum(statuses == "completed"),
#     skipped = sum(statuses == "skipped"),
#     failed = sum(statuses == "failed")
  # )
# }