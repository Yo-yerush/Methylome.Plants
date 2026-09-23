chromosome_plot <- function(DMRs_trimmed,var1_name,var2_name,chr.n,ymax,ymin,pch,
                            chr_label = chr.n) {
  chr.vec = DMRs_trimmed[DMRs_trimmed@seqnames == chr.n]

  if (length(chr.vec) == 0L) {
    plot(NA, NA,
         xlim = c(0, 1), ylim = c(ymin, ymax),
         axes = FALSE, type = "n", xlab = "", ylab = "", main = "")
    mtext(chr_label, side = 1, line = 0, cex = 1)
    return(invisible(NULL))
  }

  methylationProfiles <- GRangesList("var" = chr.vec)
  legend_text = paste0(var2_name," vs ",var1_name)
  
  col = c("#302c2c","#f54545")
  #pch = rep(20,length(methylationProfiles))
  lty = rep(0,length(methylationProfiles))
  col_legend = col[1:length(methylationProfiles)]
  lty_legend = lty[1:length(methylationProfiles)]
  
  pos = (start(methylationProfiles[[1]]) + end(methylationProfiles[[1]]))/2
  plot(pos, methylationProfiles[[1]]$plog, type = "o",
       ylim = c(ymin, ymax), 
       col = col[1], pch = pch, lty = lty[1], 
       yaxt = "n", xaxt = "n", 
       #axes = F,
       main = ""
  )
  mtext(chr_label, side = 1, line = 0, cex = 1)#adj = 0,
}

##################################################

ChrPlots_DMRs <- function(comparison_name, DMRs_gr_pl, var1, var2, context, scripts_dir) {
  chr_labels <- as.character(seqlevels(DMRs_gr_pl))
  chr_amount = length(chr_labels)
  if (chr_amount == 0L || length(DMRs_gr_pl) == 0L) return(invisible(NULL))
  
  DMRs_gr_pl$plog = -log10(DMRs_gr_pl$pValue)
  DMRs_gr_pl$plog[which(DMRs_gr_pl$plog == Inf)] = 300
  DMRs_gr_pl$plog[which(DMRs_gr_pl$plog > 300)] = 300
  DMRs_gr_pl$plog[DMRs_gr_pl$regionType == "loss"] = -DMRs_gr_pl$plog[DMRs_gr_pl$regionType == "loss"]
  gain_plog <- DMRs_gr_pl[DMRs_gr_pl$regionType == "gain",]$plog
  loss_plog <- DMRs_gr_pl[DMRs_gr_pl$regionType == "loss",]$plog
  p.max = if (length(gain_plog)) round(max(gain_plog, na.rm = TRUE)) else 0
  p.min = if (length(loss_plog)) round(min(loss_plog, na.rm = TRUE)) else 0
  if (!is.finite(p.max)) p.max <- 0
  if (!is.finite(p.min)) p.min <- 0
  if (p.max == p.min) {
    p.max <- p.max + 1
    p.min <- p.min - 1
  }
  
  #### change seqnames for plots
  DMRs_gr_pl_trimmed = renameSeqlevels(
    DMRs_gr_pl,
    stats::setNames(as.character(seq_along(chr_labels)), chr_labels)
  )



    chr_widths <- tapply(end(DMRs_gr_pl), seqnames(DMRs_gr_pl), max)
    chr_widths <- chr_widths[chr_labels]
    chr_widths[!is.finite(chr_widths)] <- 1
    chr_widths <- as.numeric(chr_widths) / max(chr_widths, na.rm = TRUE)
    fig_width <- max(8, min(14, 2.1 + 0.55 * chr_amount))

    img_device(paste0("ChrPlot_DMRs_",context,"_",comparison_name), w = fig_width, h = 2)
    layout(matrix(seq_len(chr_amount + 2), nrow = 1),
           widths = c(1.2, chr_widths, 1.6))
    
    #left axis
    par(mar = c(1,4,0.5,0))
    plot(runif(10),runif(10),xlim=c(0, 0.01), 
         ylim=c(p.min,p.max),axes=FALSE,type="n",ylab=paste(context," DMRs (-log[p])"),xlab="")
    axis(2, c(p.min, 0, p.max), labels = c(-p.min, 0, p.max))
    
    # main ChrPlot
    for (chr_number in seq(chr_amount)) {
      par(mar = c(1,0.1,0.5,0.1))
      chromosome_plot(DMRs_gr_pl_trimmed,var1,var2,chr_number,p.max,p.min,16, chr_labels[chr_number])
    }
    
    # right labels and legend
    par(mar = c(1,0,0.5,3))
    plot(NA, NA, xlim = c(0, 1), ylim = c(p.min, p.max),
         axes = FALSE, type = "n", xlab = "", ylab = "")
    axis(4, c(p.min, 0, p.max), labels = F)
    axis(4, c(p.min/2, p.max/2), labels = c("Loss", "Gain"), tick = F)
    legend("topright", legend = c(paste0(var2," vs ",var1)), 
           lty = 0, col = "#302c2c", 
           pch = 16, bty = "n")
    dev.off()
}
