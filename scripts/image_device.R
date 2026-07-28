  img_device <- function(filename, w, h, img_type = "pdf") {
    dev_call <- ifelse(img_type == "pdf", "cairo_pdf", img_type)
    img_env <- get(dev_call, envir = asNamespace("grDevices"))
    full_file_name <- paste0(filename, ".", img_type)

    if (img_type == "svg" | img_type == "pdf") {
      img_env(full_file_name, width = w, height = h, family = "serif")
    } else if (img_type == "tiff") {
      img_env(full_file_name, width = w, height = h, units = "in", res = 900, family = "serif", compression = "lzw")
    } else {
      img_env(full_file_name, width = w, height = h, units = "in", res = 900, family = "serif")
    }
  }