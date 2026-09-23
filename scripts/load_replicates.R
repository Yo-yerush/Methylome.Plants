load_replicates <- function(var_path, n.cores, var_name, metaPlot_script = F, file_type = "CX_report", convert_2_CX = F) {
  n.cores.f <- ifelse(n.cores > 2, length(var_path), 1)
  ### load the data for replicates
  methylationDataList <- mclapply(var_path,
    function(meth_file_path) {
      tryCatch(
        {
          cat(paste0("read ", basename(meth_file_path), " file...   \t[", var_name, "]\n"))
          if (file_type == "CX_report") {
            cx_4_read <- read_CX(meth_file_path, var_name)
          } else if (convert_2_CX) {
            cx_4_read <- read_CX(convert2cx_fun(meth_file_path, file_type))
          } else if (file_type == "bedMethyl") {
            cx_4_read <- read_bedMethyl(meth_file_path, var_name)
          } else if (file_type == "CGmap") {
            cx_4_read <- read_CGmap(meth_file_path, var_name)
          }
        },
        error = function(cond) {
          cat(paste0("\n* Error reading ", meth_file_path, " file: ", cond$message, "\n"))
          cx_4_read <- read_bismark(meth_file_path)
        }
      )
    },
    mc.cores = n.cores.f
  )

  # Pooling Methylation Datasets
  cat("pooling datasets\t[", var_name, "]\n", sep = "")
  methylationData_pool <- poolMethylationDatasets(GRangesList(methylationDataList))

  # Joining Replicates
  if (!metaPlot_script) {
    if (length(var_path) == 1) {
      methylationDataReplicates <- methylationDataList[[1]] # 1 sample
      names(mcols(methylationDataReplicates))[names(mcols(methylationDataReplicates)) %in% c("readsM", "readsN")] <- c("readsM1", "readsN1")
    } else {
      cat("joining ", length(var_path), " replicates\t[", var_name, "]\n", sep = "")
      # message(paste0(time_msg(), "\t----"))
      message(paste0(time_msg(), "\tjoin replicates   [", var_name, "]"))
      methylationDataReplicates <- joinReplicates(methylationDataList[[1]], methylationDataList[[2]]) # 2 samples
      if (length(var_path) > 2) {
        for (i in 3:length(var_path)) {
          methylationDataReplicates <- joinReplicates(methylationDataReplicates, methylationDataList[[i]]) # >2 samples
        }
      }
    }

    return(list(
      methylationData_pool = methylationData_pool,
      methylationDataReplicates = methylationDataReplicates
    ))
  } else {
    return(methylationData_pool)
  }
}

#################################################################

convert2cx_fun <- function(meth_file, file_type) {
  dir.create("CX_report_files", showWarnings = F)
  prefix_4_cx <- sub(paste0("\\.", file_type, "$"), "", basename(meth_file))

  # convert to 'CX_report' if the files type is 'bedMeyhl' or 'CGmap'
  if (file_type == "bedMethyl") {
    cat("convert 'bedMethyl' to 'CX_report'...\n")
    system(paste0("scripts/bedmethyl_2_cx.sh -i ", meth_file, " -o CX_report_files/", prefix_4_cx))
  } else if (file_type == "CGmap") {
    cat("convert 'CGmap' to 'CX_report'...\n")
    system(paste0("scripts/cgmap_2_cx.sh -i ", meth_file, " -o CX_report_files/", prefix_4_cx))
  }
  cx_name_out <- paste0("CX_report_files/", prefix_4_cx, ".CX_report.txt")
  cat("done:", cx_name_out, "\n")
  cx_name_out
}

#################################################################

read_bismark <- function(cx_path, name) {
  msg_0 <- capture.output(rb_out <- readBismark(cx_path), type = "message")
  if (length(msg_0) > 0 && length(strsplit(msg_0, " ")[[1]]) >= 2) {
    msg <- strsplit(msg_0, " ")[[1]][2]
  } else {
    msg <- "unknown"
  }
  message(paste0(time_msg(), "\tread ", msg, " C's [", name, "]"))

  rb_out
}

#################################################################

read_CX <- function(cx_path, name) {
  cx_df <- fread(
    cmd = ifelse(grepl("\\.gz$", cx_path),
      sprintf("gunzip -c %s", shQuote(cx_path)),
      sprintf("cat %s", shQuote(cx_path))
    ),
    col.names = c(
      "seqnames", "pos", "strand", "readsM", "readsUnM",
      "context", "trinucleotide_context"
    ),
    showProgress = FALSE
  )
  setorder(cx_df, seqnames, pos)

  cx_gr <- GRanges(
    seqnames = factor(cx_df$seqnames, levels = unique(cx_df$seqnames)),
    ranges = IRanges(as.integer(cx_df$pos), width = 1),
    strand = factor(cx_df$strand, levels = c("+", "-")),
    context = factor(cx_df$context, levels = c("CG", "CHG", "CHH")),
    readsM = as.integer(cx_df$readsM),
    readsN = as.integer(cx_df$readsM) + as.integer(cx_df$readsUnM),
    trinucleotide_context = as.character(cx_df$trinucleotide_context)
  )
  message(paste0(time_msg(), "\tread ", nrow(cx_df), " cytosine sites\t[", name, "]"))

  cx_gr
}

#################################################################

read_CGmap <- function(cgmap_path, name) {
  cgmap_df <- fread(
    cmd = ifelse(grepl("\\.gz$", cgmap_path),
      sprintf("gunzip -c %s", shQuote(cgmap_path)),
      sprintf("cat %s", shQuote(cgmap_path))
    ),
    col.names = c(
      "seqnames", "pos", "strand", "context", "ratio", "readsM", "readsN"
    ),
    showProgress = FALSE
  )
  setorder(cgmap_df, seqnames, pos)

  cgmap_gr <- GRanges(
    seqnames = factor(cgmap_df$seqnames, levels = unique(cgmap_df$seqnames)),
    ranges = IRanges(as.integer(cgmap_df$pos), width = 1),
    strand = factor(cgmap_df$strand, levels = c("+", "-")),
    context = factor(cgmap_df$context, levels = c("CG", "CHG", "CHH")),
    readsM = as.integer(cgmap_df$readsM),
    readsN = as.integer(cgmap_df$readsN),
    trinucleotide_context = as.character(rep("-", nrow(cgmap_df)))
  )
  message(paste0(time_msg(), "\tread ", nrow(cgmap_df), " cytosine sites\t[", name, "]"))

  cgmap_gr
}

#################################################################
############ continue!!!!!!
read_bedMethyl <- function(bed_path, name) {
  bed_df <- fread(
    cmd = ifelse(grepl("\\.gz$", bed_path),
      sprintf("gunzip -c %s", shQuote(bed_path)),
      sprintf("cat %s", shQuote(bed_path))
    ),
    col.names = c(
      "seqnames", "start0", "end", "strand",
      "context", "ratio", "readsM", "readsN", "trinuc"
    ),
    colClasses = list(
      character = c("seqnames", "strand", "context", "trinuc"),
      integer   = c("start0", "end", "readsM", "readsN"),
      numeric   = "ratio"
    ),
    showProgress = FALSE
  )

  ## Convert 0-based BED start → 1-based genomic coordinate
  bed_df[, pos := start0 + 1L]
  setorder(bed_df, seqnames, pos)

  bed_gr <- GRanges(
    seqnames = factor(bed_df$seqnames, levels = unique(bed_df$seqnames)),
    ranges = IRanges(bed_df$pos, width = 1),
    strand = factor(bed_df$strand, levels = c("+", "-")),
    context = factor(bed_df$context, levels = c("CG", "CHG", "CHH")),
    readsM = bed_df$readsM,
    readsN = bed_df$readsN,
    ratio = bed_df$ratio,
    trinucleotide_context = if ("trinuc" %in% names(bed_df)) {
      bed_df$trinuc
    } else {
      "-"
    }
  )
  message(paste0(time_msg(), "\tread ", nrow(bed_df), " cytosine sites\t[", name, "]"))

  bed_gr
}
