#!/usr/bin/env Rscript

## take the arguments from environment
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=4) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)}


setwd(getwd())

require(ggplot2)
require(reshape2) # use to convert data.frame from wide to long
#detach(package:reshape2)

dir.create("AT_count")

# open CBI reference
file1 = paste("ref/", args[2], "_CBI_CAIavg.txt", sep = "")

## function to draw correlation plot with polyA internal/UTR signal
r = read.table(file1, sep = "\t", header = T)

cor.graph = function(sam)    {
  file = paste(sam, "/CDS_UTR_PolyA_counts.txt", sep = "")
  data = read.table(file, header = T, sep = "\t", quote = "")
  df = merge(data, r, by.x = "Gene", by.y = "transcription_id")
  df$ratio = df$CDS_counts/df$UTR_counts
  #df$ratio_utr5 = df$UTR5_counts/df$UTR_counts
  df$ratio_norm = df$CDS_counts*1000/df$CDS_size/df$UTR_counts
  df$rpm = df$UTR_counts*1000000/sum(df$UTR_counts)
  df$cds_dens = log10(df$CDS_counts*1000/df$CDS_size+1)
  df$intron_dens = log10(df$intron_counts*1000/df$CDS_intron_size+1)

  df = df[df$rpm > 10,]
  df = na.omit(df)
  
  setEPS()
  postscript(paste(sam, "/", sam, " CBI vs ratio dotplot.eps", sep=""), pointsize = 8, width = 8, height = 6)
  
  par(mfrow=c(2,3))

  cor1 = cor.test(df$CBI[df$ratio>0], log10(df$ratio[df$ratio>0]))
  cor2 = cor.test(df$CBI[df$UTR_counts > 0], log10(df$UTR_counts[df$UTR_counts > 0]))
  cor3 = cor.test(df$GC[df$ratio>0], log10(df$ratio_norm[df$ratio_norm>0]))
  cor4 = cor.test(df$avg_CAI[df$ratio_norm>0], log10(df$ratio_norm[df$ratio_norm>0]))
  cor5 = cor.test(log(df$ratio_norm[df$ratio_norm>0]), log(df$rpm[df$ratio_norm>0]))
  cor6 = cor.test(df$GC[df$rpm>0], log(df$rpm[df$rpm>1]))
  
  
  plot(df$CBI[df$ratio>0], df$ratio[df$ratio>0], pch=16, col="red", cex=0.5, main= paste(sam, "dotplot"),
       xlab="CBI ", ylab="polyA internal to UTR signal ratio", log = "y")
  abline(lm(log10(df$ratio[df$ratio>0])~df$CBI[df$ratio>0]), lty = 2, lwd = 2, col = "grey50")
  text(min(df$CBI)+0.2, max(df$ratio)*0.15, paste(" r =", round(cor1$estimate, 3), "\n n =", 
                                                  cor1$parameter, "\n P=", formatC(cor1$p.value, digits = 3)))
  
  plot(df$CBI[df$UTR_counts > 1], df$UTR_counts[df$UTR_counts > 1], pch=16, col="red", cex=0.5, main= paste(sam, "dotplot"),
       xlab="CBI ", ylab="3' UTR reads", log = "y")
  abline(lm(log10(df$UTR_counts[df$ratio>0])~df$CBI[df$ratio>0]), lty = 2, lwd = 2, col = "grey50")
  text(min(df$CBI)+0.4, max(df$UTR_counts)*0.15, paste(" r =", round(cor2$estimate, 3), "\n n =", 
                                                       cor2$parameter, "\n P=", formatC(cor2$p.value, digits = 3)))
  
  plot(df$GC[df$ratio_norm>0], df$ratio_norm[df$ratio_norm>0], pch=16, col="blue", cex=0.5, main= paste(sam, "dotplot"),
       xlab="gene ORF GC% ", ylab="polyA normalized internal (by CDS size) to UTR signal ratio", log = "y")
  abline(lm(log10(df$ratio_norm[df$ratio_norm>0])~df$GC[df$ratio_norm>0]), lty = 2, lwd = 2, col = "grey50")
  text(max(df$GC)-0.15, max(df$ratio_norm)*0.15, paste(" r =", round(cor3$estimate, 3), "\n n =", 
                                                       cor3$parameter, "\n P=", formatC(cor3$p.value, digits = 3)))
  
  plot(df$avg_CAI[df$ratio_norm>0], df$ratio_norm[df$ratio_norm>0], pch=16, col="blue", cex=0.5, main= paste(sam, "dotplot"),
       xlab="CAI ", ylab="polyA normalized internal (by CDS size) to UTR signal ratio", log = "y")
  abline(lm(log10(df$ratio_norm[df$ratio_norm>0])~df$avg_CAI[df$ratio_norm>0]), lty = 2, lwd = 2, col = "grey50")
  text(min(df$avg_CAI)*1.15, max(df$ratio_norm)*0.15, paste(" r =", round(cor4$estimate, 3), "\n n =", 
                                                            cor4$parameter, "\n P=", formatC(cor4$p.value, digits = 3)))
 
  plot(df$GC[df$rpm>0], df$rpm[df$rpm>1], pch=16, col="chartreuse3", cex=0.5, main= paste(sam, "RPM to CDS cleavage ratio"),
       xlab="gene ORF GC%", ylab="3' UTR counts", log = "y")
  abline(lm(log10(df$rpm[df$rpm>0])~df$GC[df$rpm>0]), lty = 2, lwd = 2, col = "grey50")
  text(min(df$GC[df$rpm>0])+0.05, max(df$rpm[df$rpm>0])*0.15, paste(" r =", round(cor6$estimate, 3), "\n n =", 
                                                                    cor6$parameter, "\n P=", formatC(cor6$p.value, digits = 3)))
  
  dev.off()
  
  jpeg(filename = paste(sam, "/", sam, " CDS vs CDS intron density.jpeg", sep=""),
       width = 300, height = 350, units = "px", quality = 100)
  par(mfrow=c(1,1))
  d1 = df[df$CDS_counts>0,]
  d1 = d1[,c("Gene", "cds_dens", "intron_dens")]
  names(d1) = c("Gene", "CDS", "Intron")
  d2 = melt(d1)
  boxplot(value~variable, d2, main = "Poly(A) signal density (per kb)",
          ylab = "signal density (log10)", col = c("red", "green"), notch = T )

  
  dev.off()
  
}


## function to draw line plots for cleavage sites

gcline = function(sam, pos)    {
  file = paste(sam, "/gc_", pos, ".txt", sep = "")
  data = read.table(file, header = T, sep = "\t")
  
  data = data[order(data$pos),]
  
  df =apply(data[,2:5], 1, function(x) x/sum(x))
  colnames(df) = data$pos
  
  df1 = melt(df)
  names(df1) = c("base", "position", "frequency")
  
  
  write.csv(df1[df1$base == "A",], paste("AT_count/", sam, "_", pos, "_A_.csv", sep = ""), row.names = F)
  write.csv(df1[df1$base == "T",], paste("AT_count/", sam, "_", pos, "_T_.csv", sep = ""), row.names = F)
  
  p = ggplot(df1, aes(position, frequency, group = base)) + 
    geom_line(aes(color=base)) + 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  p+xlab('position')+ylab('base frequency')+
    ggtitle(paste(sam, pos, "PAS GC content"))+
    scale_x_continuous(breaks = seq(-(ncol(df)-1)/2, (ncol(df)-1)/2, by = 10))+
    scale_y_continuous(breaks = seq(0,1,0.25))
  
  
  ggsave(filename = paste(sam,"/gc_", pos, ".eps",sep = ""), plot = last_plot(), device = "eps",
         width = 10, height = 4, units = "in", dpi = 300)
  
}

# function 3 draw stacked barchart, compare PAS signal sequence
stackbar = function(sam, pos="UTR", col = "steelblue")    {
  cds.file = paste(sam, "/polyA_CDS_", args[3],"-", args[4], "_sig_motifs.txt", sep = "")
  utr.file = paste(sam, "/polyA_UTR3_", args[3],"-", args[4], "_sig_motifs.txt", sep = "")
  intron.file = paste(sam, "/polyA_CDS_intron_", args[3],"-", args[4], "_sig_motifs.txt", sep = "")
  cds = read.table(cds.file, sep = "\t", header = T)
  utr = read.table(utr.file, sep = "\t", header = T)
  intron = read.table(intron.file, sep = "\t", header = T)
  
  cds$CDS = cds$frqency/cds[1,"total"]*100
  utr$UTR = utr$frqency/utr[1,"total"]*100
  intron$CDS_intron = intron$frqency/intron[1,"total"]*100
  
  cds = cds[,c(1,9)]
  utr = utr[,c(1,9)]
  intron = intron[,c(1,9)]
  
  df = merge(utr, cds, by.x = "motif", by.y = "motif", all = T)
  df = merge(df, intron, by.x = "motif", by.y = "motif", all = T)
  df[is.na(df)] = 0
  
  df = df[order(-df[[pos]]),]
  
  df1 = df[1:20,]
  df1 = t(df1)
  colnames(df1) = df1[1,]
  df1 = df1[-1,]
  df2 = melt(df1)
  names(df2) = c("position", "motif", "ratio")
  df2$ratio = as.numeric(as.character(df2$ratio))
  #df2$motif = as.character(df2$motif)
  
  ggplot(df2, aes(x=motif, y=ratio, fill=position)) + 
    geom_bar(stat="identity") +
    xlab("motif") +
    ylab("ratio (% over all sequences)") +
    theme(axis.title.x = element_text(face="bold", colour="#990000", size=8),
          axis.text.x  = element_text(angle=90, vjust=0.5, size=16))
  
  ggsave(filename = paste(sam,"/",pos, "_motifs.eps",sep = ""), plot = last_plot(), 
         width = 10, height = 6, units = "in", dpi = 300)
  
  df3 = df2[df2$position == pos, 2:3]
  df3 = df3[1:15,]

  ggplot(df3, aes(motif, ratio)) + 
              geom_bar(stat = "identity", fill=col)+
              theme_minimal()+
              xlab("motif") +
              ylab("ratio (% over all sequences)") +
              theme(axis.title.x = element_text(face="bold", colour="#990000", size=8),
                    axis.text.x  = element_text(angle=90, vjust=0.5, size=16)) 
    
  ggsave(filename = paste(sam,"/",pos, " ONLY_motifs.eps",sep = ""), plot = last_plot(), 
                     width = 6, height = 6, units = "in", dpi = 300)
              
  if (pos == "UTR")    {
    df4 = df2[df2$position != "CDS_intron",]
    df4$position = factor(df4$position, levels = c("CDS", "UTR"))
    #df4 = df4[order(-df4$ratio),]
    ggplot(df4, aes(x=motif, y=ratio, fill= position)) + 
      geom_bar(stat="identity") +
      xlab("motif") +
      ylab("ratio (% over all sequences)") +
      theme(axis.title.x = element_text(face="bold", colour="#990000", size=8),
            axis.text.x  = element_text(angle=90, vjust=0.5, size=16))
    
    ggsave(filename = paste(sam,"/",pos, "_motifs UTR CDS.eps",sep = ""), plot = last_plot(), 
           width = 10, height = 6, units = "in", dpi = 300)    
  }
  
  
  write.csv(df, paste(sam,"/combined motif.csv", sep = ""), row.names = F)
}



s = args[1]
cor.graph(s)
gcline(s, "UTR")
gcline(s, "CDS")
gcline(s, "CDS_intron")
stackbar(s, "UTR")
stackbar(s, "CDS", "firebrick")
stackbar(s, "CDS_intron", "darkolivegreen")

