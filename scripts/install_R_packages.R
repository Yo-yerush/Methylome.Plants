check_if_installed <- as.logical(commandArgs(trailingOnly = TRUE))

pkg_name <- c("dplyr", "tidyr", "ggplot2", "data.table", "lattice", "PeakSegDisk", "geomtextpath", "parallel", "BiocManager", "RColorBrewer", "circlize", "cowplot", "knitr", "kableExtra", "yaml")
pkg_biocond <- c("DMRcaller", "rtracklayer", "topGO", "KEGGREST", "Rgraphviz", "GenomicFeatures", "plyranges", "AnnotationDbi", "Biostrings")

if (!check_if_installed) {
  # Update 'textshaping' version to '0.4.1' if required
  if (packageVersion("textshaping") < "0.4.1") {
    tryCatch(
      {
        suppressMessages(devtools::install_github("r-lib/textshaping", ref = "v0.4.1", quiet = T))
        message(paste0("Update 'textshaping' to v", packageVersion("textshaping")))
      },
      error = function(e) {
        message("Error installing 'textshaping'")
      }
    )
  }

  # install packages if required
  cat(paste0("\rInstall R packages...\t    "))
  installed_packages <- rownames(installed.packages())
  i <- 1
  for (pkg in pkg_name) {
    if (!(pkg %in% installed_packages)) {
      try(
        {
          suppressMessages(install.packages(pkg, repos = "http://cran.r-project.org", quiet = T))
        },
        silent = T
      )
    }
    perc_val <- (i / length(c(pkg_name, pkg_biocond))) * 100
    cat(paste0("\rInstall R packages...\t", round(perc_val, 1), "% "))
    i <- i + 1
  }

  for (pkg in pkg_biocond) {
    if (!(pkg %in% installed_packages)) {
      try(
        {
          suppressMessages(BiocManager::install(pkg, quiet = T, update = F, ask = F))
        },
        silent = T
      )
    }
    perc_val <- (i / length(c(pkg_name, pkg_biocond))) * 100
    cat(paste0("\rInstall R packages...\t", round(perc_val, 1), "% "))
    i <- i + 1
  }

  cat("\n\n")
}

# Check each package if installed
c.pkg <- .packages(all.available = TRUE)
if (!check_if_installed) message("\n")
cat("\n \tChecking installed R packages:\n")
for (i in c("textshaping", pkg_name, pkg_biocond)) {
  tryCatch(
    {
      if (c.pkg[grep(paste0("^", i, "$"), c.pkg)] == i) {
        cat(paste0("*\tinstalled ", i, ": yes\n"))
        if (!check_if_installed) message(paste0("*\tinstalled ", i, ": yes"))
      }
    },
    error = function(cond) {
      cat(paste0("*\tinstalled ", i, ": no\n"))
      if (!check_if_installed) message(paste0("*\tinstalled ", i, ": no"))
    }
  )
}
