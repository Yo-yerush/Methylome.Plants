#!/usr/bin/env Rscript

ui_started <- FALSE

ui_script_dir <- function() {
  args <- commandArgs(FALSE)
  file <- args[startsWith(args, "--file=")]
  if (length(file)) dirname(normalizePath(sub("^--file=", "", file[1]))) else getwd()
}

find_pipeline <- function() {
  here <- ui_script_dir()
  candidates <- c(getwd(), here, dirname(here),
                  file.path(dirname(dirname(here)), "Methylome.Plants"))
  for (path in candidates) {
    if (file.exists(file.path(path, "scripts", "Methylome.Plants.sh"))) {
      return(normalizePath(path))
    }
  }
  ""
}

ui_defaults <- function(pipeline = find_pipeline()) {
  c(
    pipeline = pipeline, samples = "", bundle = if (nzchar(pipeline))
      file.path(pipeline, "reference_bundles", "arabidopsis_thaliana_TAIR10.yaml") else "",
    annotation = "", description = "", tes = "", fasta = "TAIR10",
    align_output = if (nzchar(pipeline)) file.path(pipeline, "bismark_CX_reports") else "",
    file_type = "CX_report", image_type = "png", cores = "8", memory = "8G",
    cg = "0.4", chg = "0.2", chh = "0.1", bin_size = "100",
    cytosines = "4", reads = "6", pvalue = "0.05",
    feature_bins = "10", random = "10000",
    align = "0", methylome = "1", dmrs = "1", qc = "1", pca = "1",
    total = "1", chromosomes = "1", te_distance = "1", annotation_counts = "1",
    tf = "1", groups = "1", go = "0", kegg = "0", strand = "0", dmv = "0",
    delta_h = "0", te_mp = "1", gene_mp = "1", feature_mp = "0"
  )
}

ui_analysis_labels <- c(
  dmrs = "Main DMR analysis", qc = "Sample QC and conversion controls",
  pca = "PCA", total = "Total methylation", chromosomes = "Chromosome profiles",
  te_distance = "TE size and centromere distance", annotation_counts = "Annotation counts",
  tf = "TFBS density", groups = "Functional gene groups", go = "GO enrichment",
  kegg = "KEGG pathways", strand = "Strand-specific DMRs", dmv = "DMVs",
  delta_h = "Experimental delta-H", te_mp = "TE metaplots",
  gene_mp = "Gene-body metaplots", feature_mp = "Gene-feature metaplots"
)

is_checked <- function(s, key) identical(unname(s[key]), "1")
quote_shell <- function(x) shQuote(as.character(x), type = "sh")

pipeline_args <- function(s, sample_file) {
  args <- c(
    file.path(s["pipeline"], "scripts", "Methylome.Plants.sh"),
    "--samples_file", sample_file, "--pipeline_path", s["pipeline"],
    "--reference_bundle", s["bundle"], "--file_type", s["file_type"],
    "--image_type", s["image_type"], "--n_cores", s["cores"],
    "--minProportionDiff_CG", s["cg"], "--minProportionDiff_CHG", s["chg"],
    "--minProportionDiff_CHH", s["chh"], "--binSize", s["bin_size"],
    "--minCytosinesCount", s["cytosines"], "--minReadsPerCytosine", s["reads"],
    "--pValueThreshold", s["pvalue"], "--MP_features_bin_size", s["feature_bins"],
    "--metaPlot_random", s["random"]
  )
  overrides <- c(annotation = "--annotation_file", description = "--description_file",
                 tes = "--TEs_file")
  for (key in names(overrides)) {
    if (nzchar(s[key])) args <- c(args, overrides[key], s[key])
  }
  flags <- c(pca = "--pca", total = "--total_methylation", chromosomes = "--CX_ChrPlot",
             te_distance = "--TEs_distance_n_size", annotation_counts = "--total_meth_ann",
             tf = "--TF_motifs", groups = "--func_groups", go = "--GO_analysis",
             kegg = "--KEGG_pathways", strand = "--strand_DMRs", dmv = "--DMVs",
             delta_h = "--dH", te_mp = "--MP_TEs", gene_mp = "--MP_Genes",
             feature_mp = "--MP_Gene_features")
  for (key in names(flags)) if (is_checked(s, key)) args <- c(args, flags[key])
  if (!is_checked(s, "dmrs")) args <- c(args, "--DMRs_off")
  if (!is_checked(s, "qc")) args <- c(args, "--QC_off")
  unname(args)
}

validate_settings <- function(s) {
  s <- trimws(s)
  if (anyNA(s) || any(grepl("[\r\n\t]", s))) stop("Settings cannot contain tabs or newlines.")
  if (!is_checked(s, "align") && !is_checked(s, "methylome")) stop("Select at least one pipeline.")
  for (key in c("pipeline", "samples")) {
    if (!nzchar(s[key]) || !file.exists(s[key])) stop("Select an existing ", key, " path.")
    s[key] <- normalizePath(s[key], mustWork = TRUE)
  }
  if (!dir.exists(s["pipeline"]) || dir.exists(s["samples"])) stop("Invalid pipeline or samples path.")
  if (any(!nzchar(Sys.which(c("bash", "nohup", "tail"))))) {
    stop("The UI needs the standard Linux commands bash, nohup and tail.")
  }
  if (is_checked(s, "methylome")) {
    if (!file.exists(file.path(s["pipeline"], "scripts", "Methylome.Plants.sh"))) {
      stop("The selected folder has no scripts/Methylome.Plants.sh.")
    }
    for (key in c("bundle", "annotation", "description", "tes")) {
      if (key == "bundle" || nzchar(s[key])) {
        if (!file.exists(s[key]) || dir.exists(s[key])) stop("File not found: ", key)
        s[key] <- normalizePath(s[key], mustWork = TRUE)
      }
    }
    if (!s["file_type"] %in% c("CX_report", "CGmap", "bedMethyl")) stop("Invalid input format.")
    if (!s["image_type"] %in% c("pdf", "svg", "png", "tiff", "jpeg", "bmp")) stop("Invalid image format.")
    if (!is_checked(s, "dmrs") && (is_checked(s, "go") || is_checked(s, "kegg"))) {
      stop("GO/KEGG need the main DMR analysis. Enable DMRs or disable GO/KEGG.")
    }
  }
  integer_keys <- c("cores", "bin_size", "cytosines", "reads", "feature_bins")
  for (key in integer_keys) {
    value <- suppressWarnings(as.numeric(s[key]))
    if (!is.finite(value) || value < 1 || value != floor(value)) stop(key, " must be a positive integer.")
  }
  for (key in c("cg", "chg", "chh", "pvalue")) {
    value <- suppressWarnings(as.numeric(s[key]))
    if (!is.finite(value) || value <= 0 || value >= 1) stop(key, " must be between 0 and 1.")
  }
  if (s["random"] != "all") {
    value <- suppressWarnings(as.numeric(s["random"]))
    if (!is.finite(value) || value < 1 || value != floor(value)) stop("Random genes must be a positive integer or all.")
  }
  if (is_checked(s, "align")) {
    if (!file.exists(file.path(s["pipeline"], "scripts", "run_bismark.sh"))) stop("run_bismark.sh was not found.")
    if (s["fasta"] != "TAIR10") {
      if (!file.exists(s["fasta"]) || dir.exists(s["fasta"])) stop("Select a reference FASTA or enter TAIR10.")
      s["fasta"] <- normalizePath(s["fasta"], mustWork = TRUE)
    }
    if (!nzchar(s["align_output"])) stop("Choose an alignment output folder.")
    s["align_output"] <- normalizePath(s["align_output"], mustWork = FALSE)
    if (!grepl("^[1-9][0-9]*[KMG]$", s["memory"])) stop("Memory must look like 8G.")
    tools <- c("bismark", "bismark_genome_preparation", "bismark_methylation_extractor")
    if (any(!nzchar(Sys.which(tools)))) stop("Bismark tools are missing from the active environment.")
  }
  sep <- if (grepl("\\.csv$", s["samples"], ignore.case = TRUE)) "," else "\t"
  samples <- utils::read.table(s["samples"], header = FALSE, sep = sep,
                             quote = "", comment.char = "", colClasses = "character")
  if (ncol(samples) != 2L || !nrow(samples) || anyNA(samples)) {
    stop("Use a two-column sample table without a header: condition/sample, file path.")
  }
  samples[] <- lapply(samples, trimws)
  if (any(!nzchar(as.matrix(samples)))) stop("The sample table contains empty values.")
  if (any(!grepl("^[A-Za-z0-9_-]+$", samples[[1]]))) {
    stop("Use letters, numbers, underscores and hyphens for sample/condition names.")
  }
  for (i in seq_len(nrow(samples))) {
    path <- samples[i, 2]
    if (!startsWith(path, "/")) path <- file.path(dirname(s["samples"]), path)
    if (!file.exists(path) || dir.exists(path)) stop("Sample file not found: ", path)
    samples[i, 2] <- normalizePath(path, mustWork = TRUE)
  }
  if (is_checked(s, "align") && any(grepl("[[:space:]]", c(
      s[c("pipeline", "fasta", "align_output")], samples[[2]])))) {
    stop("The existing Bismark script needs pipeline, FASTA, output and FASTQ paths without spaces.")
  }
  conditions <- samples[[1]]
  if (is_checked(s, "align")) conditions <- sub("[._][0-9]*$", "", conditions)
  groups <- unique(conditions)
  if (is_checked(s, "methylome")) {
    if (length(groups) != 2L) stop("Methylome analysis needs exactly two conditions: control first, treatment second.")
    if (length(rle(conditions)$values) != 2L) stop("Group control rows first, then treatment rows.")
    if (grepl(groups[1], groups[2], fixed = TRUE) || grepl(groups[2], groups[1], fixed = TRUE)) {
      stop("The existing runner needs condition names that do not contain each other.")
    }
  }
  list(settings = s, samples = samples, groups = groups)
}

make_run_script <- function(s, sample_file, run_dir) {
  status_file <- file.path(run_dir, "exit_status.txt")
  lines <- c(
    "#!/usr/bin/env bash", "set -euo pipefail",
    paste("status_file=", quote_shell(status_file), sep = ""),
    "trap 'code=$?; printf \"%s\\n\" \"$code\" > \"$status_file\"' EXIT",
    paste("cd --", quote_shell(s["pipeline"])),
    paste("samples_file=", quote_shell(sample_file), sep = "")
  )
  if (is_checked(s, "align")) {
    args <- c(file.path(s["pipeline"], "scripts", "run_bismark.sh"),
              "-s", sample_file, "-g", s["fasta"], "-o", s["align_output"],
              "-n", s["cores"], "-m", s["memory"], "--cx", "--mat")
    alignment_log <- file.path(run_dir, "bismark_stdout.log")
    lines <- c(lines, "printf 'Starting Bismark alignment\\n'",
               paste("bash", paste(quote_shell(args), collapse = " "),
                     "| tee", quote_shell(alignment_log)))
    if (is_checked(s, "methylome")) {
      lines <- c(lines, paste0("samples_file=$(tail -n 1 -- ", quote_shell(alignment_log), ")"),
                 "[[ -f \"$samples_file\" ]] || { printf 'Bismark did not produce a sample table\\n' >&2; exit 1; }")
      s["file_type"] <- "CX_report"
    }
  }
  if (is_checked(s, "methylome")) {
    args <- pipeline_args(s, sample_file)
    quoted <- quote_shell(args)
    quoted[which(args == "--samples_file") + 1L] <- '"$samples_file"'
    lines <- c(lines, "printf 'Starting methylome analysis\\n'",
               paste("bash", paste(quoted, collapse = " ")))
  }
  c(lines, "printf 'Runner finished. Check pipeline logs for analysis errors.\\n'")
}

launch_ui <- function(pipeline = find_pipeline()) {
  if (!isTRUE(capabilities("tcltk"))) stop("This R build has no Tcl/Tk support.")
  if (!nzchar(Sys.getenv("DISPLAY"))) {
    stop("No DISPLAY is available. Use a Linux desktop or X11 forwarding; whiptail remains available.")
  }
  suppressPackageStartupMessages(library(tcltk))
  if (!isTRUE(tcltk::.TkUp)) stop("Tk could not open the graphical display.")

  ####################################################
  # theme
  tcltk::tcl("ttk::style", "theme", "use", "classic")
  ####################################################

  s <- ui_defaults(pipeline)
  vars <- lapply(s, tcltk::tclVar)
  state <- new.env(parent = emptyenv())
  state$running <- FALSE
  state$run_dir <- ""
  state$run_script <- NULL
  state$pipeline_log <- ""
  state$timer <- NULL
  state$bismark_widgets <- list()
  state$methylome_widgets <- list()
  window <- tcltk::tktoplevel()
  tcltk::tkwm.title(window, "Methylome.Plants")
  tcltk::tkwm.geometry(window, "1040x800")
  tcltk::tkwm.minsize(window, 800, 650)
  tcltk::tkgrid.columnconfigure(window, 0, weight = 1)
  tcltk::tkgrid.rowconfigure(window, 1, weight = 1)
  title <- tcltk::ttklabel(window, text = "Methylome.Plants", font = "TkHeadingFont")
  tcltk::tkgrid(title, padx = 12, pady = 8, sticky = "w")
  notebook <- tcltk::ttknotebook(window)
  tcltk::tkgrid(notebook, sticky = "nsew", padx = 12)
  tabs <- list()
  for (label in c("Inputs", "Analyses", "Parameters")) {
    tab <- tcltk::ttkframe(notebook, padding = 12)
    tcltk::tkadd(notebook, tab, text = label)
    tabs[[label]] <- tab
    tcltk::tkgrid.columnconfigure(tab, 1, weight = 1)
  }
  status <- tcltk::tclVar("Select the pipeline folder and sample table. Preview before running.")
  footer <- tcltk::ttkframe(window, padding = 10)
  tcltk::tkgrid(footer, sticky = "ew")
  tcltk::tkgrid.columnconfigure(footer, 0, weight = 1)
  tcltk::tkgrid(tcltk::ttklabel(footer, textvariable = status, wraplength = 680),
               row = 0, column = 0, sticky = "w")

  get_settings <- function() {
    values <- vapply(vars, function(x) as.character(tcltk::tclvalue(x)), character(1))
    trimws(values)
  }
  update_bismark_fields <- function() {
    enabled <- if (as.character(tcltk::tclvalue(vars$align)) == "1") "normal" else "disabled"
    for (widget in state$bismark_widgets) tcltk::tkconfigure(widget, state = enabled)
  }
  update_methylome_fields <- function() {
    enabled <- if (as.character(tcltk::tclvalue(vars$methylome)) == "1") "normal" else "disabled"
    for (widget in state$methylome_widgets) tcltk::tkconfigure(widget, state = enabled)
  }
  set_settings <- function(values) {
    for (key in intersect(names(values), names(vars))) tcltk::tclvalue(vars[[key]]) <- values[key]
    update_bismark_fields()
    update_methylome_fields()
  }
  error_box <- function(e) tcltk::tkmessageBox(
    parent = window, icon = "error", type = "ok", title = "Please check",
    message = conditionMessage(e)
  )
  protect <- function(action) function() tryCatch(action(), error = error_box)
  create_bundle <- protect(function() {
    pipeline_dir <- normalizePath(get_settings()["pipeline"], mustWork = TRUE)
    backend <- file.path(pipeline_dir, "scripts", "reference_wizard", "reference_wizard.R")
    if (!file.exists(backend)) stop("The selected pipeline folder has no reference wizard backend.")
    dialog <- tcltk::tktoplevel()
    on.exit(try(tcltk::tkdestroy(dialog), silent = TRUE), add = TRUE)
    tcltk::tkwm.title(dialog, "Create reference bundle")
    tcltk::tkwm.transient(dialog, window)
    tcltk::tkwm.geometry(dialog, "860x600")
    tcltk::tkgrid.columnconfigure(dialog, 0, weight = 1)
    tcltk::tkgrid.rowconfigure(dialog, 0, weight = 1)
    wizard_status <- tcltk::tclVar("Fill the required (*) fields. Leave optional resources empty.")
    safe <- function(action, parent = dialog) function() tryCatch(action(), error = function(e) {
      tcltk::tclvalue(wizard_status) <- "Please check the message and try again."
      tcltk::tkmessageBox(parent = parent, icon = "error", type = "ok",
                          title = "Reference setup", message = conditionMessage(e))
    })
    busy <- function(message) {
      tcltk::tclvalue(wizard_status) <- message
      tcltk::tcl("update", "idletasks")
    }
    run_r <- function(args) {
      output <- tempfile("bundle_output_")
      errors <- tempfile("bundle_errors_")
      on.exit(unlink(c(output, errors)), add = TRUE)
      code <- system2(file.path(R.home("bin"), "Rscript"), args = shQuote(args),
                      stdout = output, stderr = errors, wait = TRUE)
      if (code != 0L) {
        stop(paste(c(paste("Reference setup exited with code", code),
                     readLines(errors, warn = FALSE), readLines(output, warn = FALSE)),
                   collapse = "\n"), call. = FALSE)
      }
      readLines(output, warn = FALSE)
    }
    pages <- list(
      Reference = c(organism = "Organism *", assembly = "Assembly *",
        annotation = "Gene annotation (GTF/GFF3/CSV) *", fasta = "FASTA or .fai (optional)",
        chromosome_sizes = "Chromosome sizes (optional)",
        primary_seqlevels = "Primary chromosomes (comma-separated)",
        chloroplast_seqlevels = "Chloroplast sequences (comma-separated)",
        mitochondrial_seqlevels = "Mitochondrial sequences (comma-separated)"),
      `GO / KEGG` = c(orgdb_package = "GO OrgDb package (optional)",
        kegg_organism = "KEGG organism code (optional)",
        gene_to_go = "Custom gene-to-GO table (optional)",
        kegg_id_map = "KEGG gene-ID mapping (optional)"),
      `Optional resources` = c(descriptions = "Gene descriptions", te = "TE annotation",
        stable_gbm = "Stable gbM gene list", dynamic_gbm = "Dynamic gbM gene list",
        centromeres = "Centromere BED", heterochromatin = "Heterochromatin BED",
        tfbs = "TFBS BED", gene_sets = "Functional gene sets (TSV/GMT)",
        seqname_aliases = "Sequence aliases (alias/canonical)")
    )
    field_keys <- unlist(lapply(pages, names), use.names = FALSE)
    fields <- lapply(setNames(rep("", length(field_keys)), field_keys), tcltk::tclVar)
    file_keys <- setdiff(field_keys, c("organism", "assembly", "primary_seqlevels",
      "chloroplast_seqlevels", "mitochondrial_seqlevels", "orgdb_package", "kegg_organism"))
    find_organism <- function(key) {
      chooser <- tcltk::tktoplevel()
      on.exit(try(tcltk::tkdestroy(chooser), silent = TRUE), add = TRUE)
      tcltk::tkwm.title(chooser, if (key == "orgdb_package") "Find GO / OrgDb" else "Find KEGG organism")
      tcltk::tkwm.transient(chooser, dialog)
      query <- tcltk::tclVar(as.character(tcltk::tclvalue(fields$organism)))
      rows <- data.frame()
      tcltk::tkgrid(tcltk::ttklabel(chooser, text = "Search by organism name or package/code:"),
                   row = 0, column = 0, columnspan = 2, padx = 8, pady = 5, sticky = "w")
      tcltk::tkgrid(tcltk::ttkentry(chooser, textvariable = query, width = 65),
                   row = 1, column = 0, padx = 8, pady = 5, sticky = "ew")
      results <- tcltk::tklistbox(chooser, width = 100, height = 14, selectmode = "browse", exportselection = FALSE)
      scroll <- tcltk::ttkscrollbar(chooser, command = function(...) tcltk::tkyview(results, ...))
      tcltk::tkconfigure(results, yscrollcommand = function(...) tcltk::tkset(scroll, ...))
      tcltk::tkgrid(results, row = 2, column = 0, padx = 8, sticky = "nsew")
      tcltk::tkgrid(scroll, row = 2, column = 1, sticky = "ns")
      search <- safe(function() {
        busy("Searching organism databases...")
        cache_dir <- file.path(pipeline_dir, "reference_bundles", "cache")
        dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
        args <- if (key == "orgdb_package") {
          c(backend, "list-orgdb", "--include-available", "--cache", file.path(cache_dir, "orgdb_packages.tsv"))
        } else {
          c(backend, "list-kegg", "--cache", file.path(cache_dir, "kegg_organisms.tsv"))
        }
        term <- trimws(as.character(tcltk::tclvalue(query)))
        if (nzchar(term)) args <- c(args, "--query", term)
        lines <- run_r(args)
        rows <<- if (length(lines)) utils::read.delim(text = paste(lines, collapse = "\n"),
          header = FALSE, quote = "", comment.char = "", colClasses = "character") else data.frame()
        tcltk::tkdelete(results, 0, "end")
        if (nrow(rows)) {
          for (label in apply(rows, 1, paste, collapse = " | ")) tcltk::tkinsert(results, "end", label)
        }
        tcltk::tclvalue(wizard_status) <- paste(nrow(rows), "organism database matches found.")
      }, parent = chooser)
      tcltk::tkgrid(tcltk::ttkbutton(chooser, text = "Search", command = search),
                   row = 1, column = 1, padx = 8)
      select <- safe(function() {
        selected <- as.character(tcltk::tclvalue(tcltk::tkcurselection(results)))
        if (!nzchar(selected)) stop("Select a result first.")
        tcltk::tclvalue(fields[[key]]) <- rows[as.integer(selected) + 1L, 1]
        tcltk::tkdestroy(chooser)
      }, parent = chooser)
      tcltk::tkgrid(tcltk::ttkbutton(chooser, text = "Select", command = select),
                   row = 3, column = 0, padx = 8, pady = 8, sticky = "w")
      tcltk::tkgrid(tcltk::ttkbutton(chooser, text = "Cancel", command = function() tcltk::tkdestroy(chooser)),
                   row = 3, column = 1, padx = 8, pady = 8)
      tcltk::tcl("update", "idletasks")
      tcltk::tkgrab.set(chooser)
      search()
      tcltk::tkwait.window(chooser)
      tcltk::tkgrab.set(dialog)
    }
    field <- function(parent, row, key, label) {
      tcltk::tkgrid(tcltk::ttklabel(parent, text = label), row = row, column = 0,
                   padx = 5, pady = 5, sticky = "w")
      tcltk::tkgrid(tcltk::ttkentry(parent, textvariable = fields[[key]], width = 48),
                   row = row, column = 1, padx = 5, pady = 5, sticky = "ew")
      if (key %in% file_keys) {
        button <- tcltk::ttkbutton(parent, text = "Browse…", command = safe(function() {
          path <- as.character(tcltk::tclvalue(tcltk::tkgetOpenFile(parent = dialog, title = label)))
          if (nzchar(path)) tcltk::tclvalue(fields[[key]]) <- path
        }))
      } else if (key %in% c("orgdb_package", "kegg_organism")) {
        button <- tcltk::ttkbutton(parent, text = "Find…", command = safe(function() find_organism(key)))
      } else return(invisible(NULL))
      tcltk::tkgrid(button, row = row, column = 2, padx = 5)
    }
    notebook <- tcltk::ttknotebook(dialog)
    tcltk::tkgrid(notebook, row = 0, column = 0, padx = 10, pady = 10, sticky = "nsew")
    tabs <- list()
    for (label in names(pages)) {
      tab <- tcltk::ttkframe(notebook, padding = 10)
      tabs[[label]] <- tab
      tcltk::tkadd(notebook, tab, text = label)
      tcltk::tkgrid.columnconfigure(tab, 1, weight = 1)
      for (i in seq_along(pages[[label]])) field(tab, i - 1L, names(pages[[label]])[i], pages[[label]][i])
    }
    orgdb_lengths <- tcltk::tclVar("1")
    tcltk::tkgrid(tcltk::ttkcheckbutton(tabs$Reference,
      text = "Use OrgDb chromosome lengths when no FASTA/size table is supplied", variable = orgdb_lengths),
      row = length(pages$Reference), column = 0, columnspan = 3, pady = 10, sticky = "w")
    tcltk::tkgrid(tcltk::ttklabel(tabs[["GO / KEGG"]], wraplength = 740,
      text = "GO and KEGG are independent. GO gene-ID compatibility is checked when saving. Leave either selection empty to skip it."),
      row = length(pages[["GO / KEGG"]]), column = 0, columnspan = 3, pady = 10, sticky = "w")
    tcltk::tkgrid(tcltk::ttklabel(dialog, textvariable = wizard_status, wraplength = 800),
                 row = 1, column = 0, padx = 15, pady = 5, sticky = "w")
    generate <- safe(function() {
      values <- trimws(vapply(fields, function(x) as.character(tcltk::tclvalue(x)), character(1)))
      if (any(!nzchar(values[c("organism", "assembly", "annotation")]))) {
        stop("Organism, assembly, and gene annotation are required.")
      }
      package <- unname(values["orgdb_package"])
      if (nzchar(package) && !requireNamespace(package, quietly = TRUE)) {
        answer <- tcltk::tkmessageBox(parent = dialog, icon = "question", type = "yesno", default = "no",
          title = "Install OrgDb?", message = paste(package, "is not installed. Install it with BiocManager?"))
        if (as.character(tcltk::tclvalue(answer)) != "yes") return(invisible(NULL))
        busy(paste("Installing", package, "..."))
        run_r(c("-e", paste(
          "pkg <- commandArgs(TRUE)[1];",
          "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org');",
          "BiocManager::install(pkg, ask=FALSE, update=FALSE); loadNamespace(pkg)"
        ), package))
      }
      output_dir <- file.path(pipeline_dir, "reference_bundles", "generated")
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      filename <- paste0(gsub("[^A-Za-z0-9._-]", "_", paste(values["organism"], values["assembly"], sep = "_")), ".yaml")
      output <- as.character(tcltk::tclvalue(tcltk::tkgetSaveFile(parent = dialog,
        title = "Save reference bundle", initialdir = output_dir, initialfile = filename, defaultextension = ".yaml")))
      if (!nzchar(output)) return(invisible(NULL))
      args <- c(backend, "generate", "--output", output)
      for (key in names(values)[nzchar(values)]) {
        args <- c(args, paste0("--", gsub("_", "-", key, fixed = TRUE)), unname(values[key]))
      }
      if (as.character(tcltk::tclvalue(orgdb_lengths)) != "1") args <- c(args, "--disable-orgdb-lengths")
      busy("Validating reference resources and saving the bundle...")
      summary <- run_r(args)
      bundle <- sub("^Bundle\t", "", grep("^Bundle\t", summary, value = TRUE))
      if (length(bundle) != 1L || !file.exists(bundle)) stop("The reference backend did not return a bundle file.")
      set_settings(c(bundle = normalizePath(bundle), annotation = "", description = "", tes = ""))
      tcltk::tclvalue(status) <- paste("Reference bundle created:", bundle)
      tcltk::tkmessageBox(parent = dialog, icon = "info", type = "ok", title = "Reference ready",
        message = paste(gsub("\t", ": ", summary, fixed = TRUE), collapse = "\n"))
      tcltk::tkdestroy(dialog)
    })
    buttons <- tcltk::ttkframe(dialog)
    tcltk::tkgrid(buttons, row = 2, column = 0, padx = 15, pady = 10, sticky = "w")
    tcltk::tkpack(tcltk::ttkbutton(buttons, text = "Save bundle", command = generate), side = "left", padx = 5)
    tcltk::tkpack(tcltk::ttkbutton(buttons, text = "Cancel", command = function() tcltk::tkdestroy(dialog)), side = "left", padx = 5)
    tcltk::tcl("update", "idletasks")
    tcltk::tkgrab.set(dialog)
    tcltk::tkwait.window(dialog)
  })
  add_field <- function(parent, row, label, key, browse = NULL, choices = NULL) {
    tcltk::tkgrid(tcltk::ttklabel(parent, text = label), row = row, column = 0,
                 padx = 5, pady = 5, sticky = "w")
    widget <- if (is.null(choices)) {
      tcltk::ttkentry(parent, textvariable = vars[[key]], width = 58)
    } else {
      tcltk::ttkcombobox(parent, textvariable = vars[[key]], values = choices,
                        state = "readonly", width = 20)
    }
    tcltk::tkgrid(widget, row = row, column = 1, padx = 5, pady = 5, sticky = "ew")
    widgets <- list(widget)
    if (!is.null(browse)) {
      button <- tcltk::ttkbutton(parent, text = "Browse…", command = protect(function() {
        path <- if (browse == "directory") {
          tcltk::tkchooseDirectory(parent = window, title = label, mustexist = TRUE)
        } else {
          tcltk::tkgetOpenFile(parent = window, title = label)
        }
        path <- as.character(tcltk::tclvalue(path))
        if (nzchar(path)) {
          tcltk::tclvalue(vars[[key]]) <- path
          if (key == "pipeline") {
            tcltk::tclvalue(vars$bundle) <- file.path(path, "reference_bundles", "arabidopsis_thaliana_TAIR10.yaml")
            tcltk::tclvalue(vars$align_output) <- file.path(path, "bismark_CX_reports")
          }
        }
      }))
      tcltk::tkgrid(button, row = row, column = 2, padx = 5)
      widgets <- c(widgets, list(button))
    }
    if (key %in% c("fasta", "align_output", "memory")) {
      state$bismark_widgets <- c(state$bismark_widgets, widgets)
    }
    if (key %in% c("cg", "chg", "chh", "bin_size", "cytosines", "reads",
                   "pvalue", "feature_bins", "random")) {
      state$methylome_widgets <- c(state$methylome_widgets, widgets)
    }
  }
  input <- tabs$Inputs
  add_field(input, 0, "Pipeline folder", "pipeline", "directory")
  add_field(input, 1, "Sample table (two columns, no header)", "samples", "file")
  add_field(input, 2, "Reference bundle YAML", "bundle", "file")
  tcltk::tkgrid(tcltk::ttkbutton(input, text = "Create bundle…", command = create_bundle),
               row = 2, column = 3, padx = 5, sticky = "w")
  add_field(input, 3, "Gene annotation override (optional)", "annotation", "file")
  add_field(input, 4, "Description override (optional)", "description", "file")
  add_field(input, 5, "TE annotation override (optional)", "tes", "file")
  add_field(input, 6, "Methylation input format", "file_type", choices = c("CX_report", "CGmap", "bedMethyl"))
  add_field(input, 7, "Image format", "image_type", choices = c("png", "pdf", "svg", "tiff", "jpeg", "bmp"))
  tcltk::tkgrid(tcltk::ttkcheckbutton(input, text = "Run Bismark: FASTQ → CX reports",
                                    variable = vars$align, command = update_bismark_fields),
               row = 8, column = 0, columnspan = 3, sticky = "w", pady = c(24, 6))
  add_field(input, 9, "Bismark FASTA (or TAIR10)", "fasta", "file")
  add_field(input, 10, "Bismark output folder", "align_output", "directory")
  tcltk::tkgrid(tcltk::ttklabel(input, wraplength = 870, text = paste(
    "Reference resources come from the selected YAML; optional overrides may be left empty.",
    "For methylation input, put control rows first and treatment rows second.",
    "For FASTQ alignment, use names such as wt_1 and mutant_1; Bismark removes the replicate suffix.",
    "Use Create bundle to prepare a new reference."
  )), row = 11, column = 0, columnspan = 3, pady = 12, sticky = "w")

  analyses <- tabs$Analyses
  tcltk::tkgrid(tcltk::ttkcheckbutton(analyses, text = "Run methylome analysis",
                                    variable = vars$methylome, command = update_methylome_fields),
               row = 0, column = 0, sticky = "w", pady = 6)
  keys <- names(ui_analysis_labels)
  for (i in seq_along(keys)) {
    checkbox <- tcltk::ttkcheckbutton(analyses, text = ui_analysis_labels[i],
                                     variable = vars[[keys[i]]])
    tcltk::tkgrid(checkbox,
                 row = 1 + (i - 1L) %% 9L, column = (i - 1L) %/% 9L,
                 padx = 6, pady = 7, sticky = "w")
    state$methylome_widgets <- c(state$methylome_widgets, list(checkbox))
  }
  preset <- function(mode) function() {
    values <- setNames(rep("0", length(keys)), keys)
    if (mode == "default") values <- ui_defaults()[keys]
    if (mode == "basic") values[c("dmrs", "qc", "pca", "total")] <- "1"
    set_settings(values)
  }
  presets <- tcltk::ttkframe(analyses)
  tcltk::tkgrid(presets, row = 11, column = 0, columnspan = 2, pady = 14, sticky = "w")
  preset_buttons <- list(
    tcltk::ttkbutton(presets, text = "Basic analyses", command = preset("basic")),
    tcltk::ttkbutton(presets, text = "All analyses off", command = preset("off")),
    tcltk::ttkbutton(presets, text = "Original UI defaults", command = preset("default"))
  )
  for (button in preset_buttons) tcltk::tkpack(button, side = "left", padx = 5)
  state$methylome_widgets <- c(state$methylome_widgets, preset_buttons)

  parameters <- tabs$Parameters
  labels <- c(cores = "CPU cores", memory = "Bismark extractor memory",
              cg = "CG minimum methylation difference", chg = "CHG minimum difference",
              chh = "CHH minimum difference", bin_size = "DMR bin size (bp)",
              cytosines = "Minimum cytosines per bin", reads = "Minimum reads per cytosine",
              pvalue = "DMR adjusted-score cutoff", feature_bins = "Gene-feature bin size",
              random = "Metaplot sample count (or all)")
  for (i in seq_along(labels)) add_field(parameters, i - 1L, labels[i], names(labels)[i])
  update_bismark_fields()
  update_methylome_fields()

  preview <- protect(function() {
    checked <- validate_settings(get_settings())
    s <- checked$settings
    commands <- make_run_script(s, "<validated sample table>", "<run directory>")
    cat(paste(c(paste("Conditions / sample groups:", paste(checked$groups, collapse = ", ")),
                paste("Sample-table rows:", nrow(checked$samples)), "",
                "Commands that will run:", commands, "",
                "The runner uses the existing pipeline. Outputs follow its existing folder rules."),
              collapse = "\n"), "\n")
    tcltk::tclvalue(status) <- "Inputs validated. Review commands in the console before pressing Run."
  })
  start <- protect(function() {
    if (state$running) stop("A pipeline run is already active in this window.")
    checked <- validate_settings(get_settings())
    s <- checked$settings
    log_root <- file.path(s["pipeline"], "results", "ui_runs")
    dir.create(log_root, recursive = TRUE, showWarnings = FALSE)
    run_dir <- tempfile(paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = log_root)
    if (!dir.create(run_dir)) stop("Could not create the run directory.")
    sample_file <- file.path(run_dir, "samples_validated.tsv")
    utils::write.table(checked$samples, sample_file, sep = "\t", quote = FALSE,
                       row.names = FALSE, col.names = FALSE)
    utils::write.table(data.frame(setting = names(s), value = unname(s)),
                       file.path(run_dir, "settings.tsv"), sep = "\t",
                       row.names = FALSE, quote = TRUE)
    script <- file.path(run_dir, "run_pipeline.sh")
    writeLines(make_run_script(s, sample_file, run_dir), script)
    state$run_script <- script
    tcltk::tkdestroy(window)
  })
  save_settings <- protect(function() {
    path <- as.character(tcltk::tclvalue(tcltk::tkgetSaveFile(
      parent = window, title = "Save UI settings", defaultextension = ".tsv"
    )))
    if (nzchar(path)) {
      values <- get_settings()
      utils::write.table(data.frame(setting = names(values), value = unname(values)),
                         path, sep = "\t", quote = TRUE, row.names = FALSE)
      tcltk::tclvalue(status) <- paste("Settings saved:", path)
    }
  })
  load_settings <- protect(function() {
    path <- as.character(tcltk::tclvalue(tcltk::tkgetOpenFile(parent = window, title = "Load UI settings")))
    if (nzchar(path)) {
      values <- utils::read.delim(path, colClasses = "character", check.names = FALSE)
      if (!identical(names(values), c("setting", "value")) || anyDuplicated(values$setting)) stop("Invalid settings table.")
      set_settings(setNames(values$value, values$setting))
      tcltk::tclvalue(status) <- paste("Settings loaded:", path)
    }
  })
  buttons <- tcltk::ttkframe(footer)
  tcltk::tkgrid(buttons, row = 1, column = 0, columnspan = 2, pady = 8, sticky = "w")
  tcltk::tkpack(tcltk::ttkbutton(buttons, text = "Load settings", command = load_settings), side = "left", padx = 4)
  tcltk::tkpack(tcltk::ttkbutton(buttons, text = "Save settings", command = save_settings), side = "left", padx = 4)
  tcltk::tkpack(tcltk::ttkbutton(buttons, text = "Validate / preview", command = preview), side = "left", padx = 4)
  run_button <- tcltk::ttkbutton(buttons, text = "Run pipeline", command = start)
  tcltk::tkpack(run_button, side = "left", padx = 12)
  close_window <- function() {
    if (state$running) {
      answer <- tcltk::tkmessageBox(parent = window, type = "yesno", default = "no",
                                   icon = "question", title = "Close UI?",
                                   message = paste("The pipeline will KEEP RUNNING. Its log is in:", state$run_dir,
                                                   "\nClose the UI?"))
      if (as.character(tcltk::tclvalue(answer)) != "yes") return(invisible(NULL))
    }
    if (!is.null(state$timer)) try(tcltk::tcl("after", "cancel", state$timer), silent = TRUE)
    tcltk::tkdestroy(window)
  }
  tcltk::tkwm.protocol(window, "WM_DELETE_WINDOW", close_window)
  
  ui_started <<- TRUE

  tcltk::tkwait.window(window)
  if (!is.null(state$run_script)) {
    exit_code <- system2("bash", args = shQuote(state$run_script), wait = TRUE)

    # Reserve 78 for UI startup failures.
    quit(status = if (exit_code == 78L) 1L else exit_code)
  }
}

ui_main <- function(args = commandArgs(TRUE)) {
  pipeline <- find_pipeline()
  if ("--help" %in% args || "-h" %in% args) {
    cat("Usage: Rscript Methylome.Plants_UI_tcltk.R [--pipeline-path /path/to/Methylome.Plants] [--check]\n")
    cat("--check reports Tcl/Tk and pipeline availability without opening a window.\n")
    return(invisible(NULL))
  }
  i <- 1L
  check <- FALSE
  while (i <= length(args)) {
    if (args[i] == "--pipeline-path" && i < length(args)) {
      pipeline <- args[i + 1L]
      i <- i + 2L
    } else if (args[i] == "--check") {
      check <- TRUE
      i <- i + 1L
    } else stop("Unknown/incomplete UI argument: ", args[i])
  }
  if (check) {
    cat("R:", as.character(getRversion()), "\n")
    cat("Conda environment:", Sys.getenv("CONDA_DEFAULT_ENV"), "\n")
    cat("Tcl/Tk support:", capabilities("tcltk"), "\n")
    cat("DISPLAY:", Sys.getenv("DISPLAY", unset = "<not set>"), "\n")
    cat("Pipeline folder:", pipeline, "\n")
    cat("Pipeline script found:", file.exists(file.path(pipeline, "scripts", "Methylome.Plants.sh")), "\n")
    return(invisible(NULL))
  }
  launch_ui(pipeline)
}

if (sys.nframe() == 0L) {
  tryCatch(ui_main(), error = function(e) {
    message("UI error: ", conditionMessage(e))
    quit(status = if (ui_started) 1L else 78L)
  })
}
