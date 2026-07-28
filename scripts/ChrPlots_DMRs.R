chromosome_plot <- function(DMRs_trimmed,var1_name,var2_name,chr.n,ymax,ymin,pch,
                            chr_label = chr.n) {
  chr.vec = DMRs_trimmed[DMRs_trimmed@seqnames == chr.n]
  
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
  
  DMRs_gr_pl$plog = -log10(DMRs_gr_pl$pValue)
  DMRs_gr_pl$plog[which(DMRs_gr_pl$plog == Inf)] = 300
  DMRs_gr_pl$plog[which(DMRs_gr_pl$plog > 300)] = 300
  DMRs_gr_pl$plog[DMRs_gr_pl$regionType == "loss"] = -DMRs_gr_pl$plog[DMRs_gr_pl$regionType == "loss"]
  p.max = round(max(DMRs_gr_pl[DMRs_gr_pl$regionType == "gain",]$plog))
  p.min = round(min(DMRs_gr_pl[DMRs_gr_pl$regionType == "loss",]$plog))
  
  #### change seqnames for plots
  DMRs_gr_pl_trimmed = renameSeqlevels(
    DMRs_gr_pl,
    stats::setNames(as.character(seq_along(chr_labels)), chr_labels)
  )



    img_device(paste0("ChrPlot_DMRs_",context,"_",comparison_name), w = 8, h = 2)
    
    #left axis
    par(mar = c(1,4,2,0))
    par(fig=c(0,2,0,10)/10)
    plot(runif(10),runif(10),xlim=c(0, 0.01), 
         ylim=c(p.min,p.max),axes=FALSE,type="n",ylab=paste(context," DMRs (-log[p])"),xlab="")
    axis(2, c(p.min, 0, p.max), labels = c(-p.min, 0, p.max))
    par(new=T)
    
    # main ChrPlot
    i=1
    i.0=1.15
    u = (10-i.0)/6
    for (chr_number in seq(chr_amount)) {
      par(mar = c(1,0,2,0))
      par(fig=c(i,i+u,0,10)/10)
      chromosome_plot(DMRs_gr_pl_trimmed,var1,var2,chr_number,p.max,p.min,"•", chr_labels[chr_number])
      par(new=T)
      i=i+u
    }
    
    # right axis
    par(mar = c(1,0,2,2))
    par(fig=c(u*4+1-0.025,u*4+3-0.025,0,10)/10)
    axis(4, c(p.min, 0, p.max), labels = F, )
    par(new=T)
    
    # right axis - gain or loss
    par(mar = c(1,0,2,2))
    par(fig=c(6.6,8.6,0,10)/10) # c(u*4+1,u*4+3,0,10)/10
    axis(4, c(p.min/2, p.max/2), labels = c("Loss", "Gain"), tick = F)
    par(new=T)
    
    # legend
    par(mar = c(1,0,2,0))
    par(fig=c(8,10,0,10)/10) # c(u*4+3,10,0,10)/10)
    legend("topright", legend = c(paste0(var2," vs ",var1)), 
           lty = 0, col = "#302c2c", 
           pch = 16, bty = "n")
    dev.off()
}

