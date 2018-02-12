#!/usr/bin/env Rscript

## this script gonna generate something like line graph with confidence interval 
## output file is PAS_average.jepg

## take the arguments from environment
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default reference
  args[2] = "NC12"
}

setwd(getwd())


library(reshape2)
library(ggplot2)


intup = function(x) {
  d = t.test(x)
  d$conf.int[2]
}

intdn = function(x)  {
  d = t.test(x)
  d$conf.int[1]
}


#### end of basic function part

## ==========================================================================
## draw lines with ggplot2 to show the difference of average CAI at 10 position (or more depending on input) 
## between putative polyA adding signals and random region 
## ==========================================================================

PAS.line = function(sam)  {

  pcodon = read.table(paste(sam, "/PAS_CAI.txt", sep = ""), sep = "\t", fill = T) 
  
  pcodon = na.omit(pcodon)
  
  ## calculate the average, upper and lower confidence interval, store in d with row bind
  len = ncol(pcodon)
  PAS_avg = apply(pcodon[,4:len], 2, mean)
  PAS_up = apply(pcodon[,4:len], 2, intup)
  PAS_dn = apply(pcodon[,4:len], 2, intdn)
  
  rcodon = read.table(paste(sam, "/PAS_random_CAI.txt", sep = ""), sep = "\t", fill = T)
  
  rcodon = na.omit(rcodon)
  len = ncol(rcodon)
  ran_avg = apply(rcodon[,2:len], 2, mean)
  ran_up = apply(rcodon[,2:len], 2, intup)
  ran_dn = apply(rcodon[,2:len], 2, intdn)
  
  d = rbind(PAS_avg, PAS_up, PAS_dn, ran_avg, ran_up, ran_dn) 
  d = as.data.frame(d)
  colnames(d) = paste("AA", 1:ncol(d), sep = "")
  d$cat1 = gsub("_.+", "", rownames(d), perl = T)
  d$cat2 = gsub("^.+_", "", rownames(d), perl = T)
  
  # use melt from reshape2 package
  d = melt(d)
  dv = d[d$cat2 == "avg", ]
  names(dv)[4] = "avg"
  dv = dv[,-2]
  
  du = d[d$cat2 == "up", ]
  names(du)[4] = "up"
  du = du[,-2]
  
  dn = d[d$cat2 == "dn", ]
  names(dn)[4] = "low"
  dn = dn[,-2]
  
  d = merge(dv, du, by = c("cat1", "variable"))
  d = merge(d, dn, by = c("cat1", "variable"))
  
  d$variable = gsub("AA", "", d$variable)
  d$variable = as.numeric(as.character(d$variable))
  
  names(d)[1:2] = c("sample", "position")
  
  
  
  p = ggplot(data = d, aes(x = position, y = avg, color =sample))+ geom_line()
  p + geom_ribbon(aes(ymin = low, ymax = up, fill = sample), alpha = 0.4, color = NA) +
    scale_y_continuous(limits = c(min(d$low)-0.1, max(d$up)+0.1)) +
    ggtitle(sam) + ylab("CAI")
  
  ggsave(filename = paste(sam,"/PAS_averge.jpeg",sep = ""), plot = last_plot(), 
         width = 5, height = 4, units = "in", dpi = 300)  
  }



### ========================================
##  PAS Codon frequency and HEATMAP
### ========================================
pas.heatmap = function(sam, st, ed, out)  {
  
  cf = read.table(paste(sam, "/PAS_codons.txt", sep = ""), sep = "\t", fill = T)
  rpm = read.table(paste(sam, "/CDS_UTR_PolyA_counts.txt", sep = ""), sep = "\t", header = T)
  rpm = rpm[,1:3]
  cf = merge(rpm, cf, by.x = "Gene", by.y = "V1")
  cf$ratio = cf$V2*100/cf$UTR_counts
  #df = cf[cf$V3 >= 5 & cf$ratio > 0.1 & cf$ratio < 100, ]
  #df = cf[cf$V3 > 1,]
  df = cf
  
  #make a table to compare the composition of codons 
  x = as.data.frame(table(df$V4))
  names(x) = c("codon", "P1")
  column = 5:(ncol(df)-3)
  for (i in column)  {
    cl = paste("V", i, sep = "")
    cl2 = paste("P", i-3, sep = "")
    d1 = as.data.frame(table(df[[cl]]))
    names(d1) = c("codon", cl2)
    x = merge(x, d1, by.x = "codon", by.y = "codon", all = T)
  }
  x[is.na(x)] = 0
  x = x[!(x$codon %in% c("TGA", "TAA", "TAG", "")),]
  x$sAll = apply(x[,2:ncol(x)], 1, sum)
  
  
  
  ## draw control
  ck = read.table(paste(sam, "/PAS_random_codons.txt", sep = ""), sep = "\t", fill = T)
  ck = na.omit(ck)
  ck = ck[ck$V1 %in% df$Gene,]  # gene used for analyses must be the same group as having internal cleavage
  x1 = as.data.frame(table(ck$V2))
  names(x1) = c("codon", "P1")
  column = 3:ncol(ck)
  for (i in column)  {
    cl = paste("V", i, sep = "")
    cl2 = paste("P", i-1, sep = "")
    d1 = as.data.frame(table(ck[[cl]]))
    names(d1) = c("codon", cl2)
    x1 = merge(x1, d1, by.x = "codon", by.y = "codon")
  }
  x1$ckAll = apply(x1[,2:ncol(x1)], 1, sum)
  x1 = x1[x1$codon %in% x$codon,]
  
  
  m1 = as.matrix(x[,2:(ncol(x)-1)]) 
  row.names(m1) = x$codon
  m1 = m1/max(m1)
  m2 = as.matrix(x1[,2:(ncol(x)-1)])
  row.names(m2) = x1$codon
  m2 = m2/max(m2)
  
  mat = log2(m1/m2+1)
  
  
  ### ---------------------------------------------
  ### This part draw codon frequency by amino acid level
  ### ---------------------------------------------
  
  df = merge(x, x1, by.x = "codon", by.y = "codon")
  df$sAll = df$sAll/sum(df$sAll) 
  ref = read.table("ref/CBI_ref.txt", sep = "\t", header = T)
  df = merge(ref, df, by.x = "codon", by.y = "codon")
  ref = read.table("ref/codon_frequency.txt", sep = "\t")
  ref = ref[,-2]
  df = merge(ref, df, by.x = "V1", by.y = "codon")
  
  names(df)[1:2] = c("codon", "frequency")
  v = st:ed # region, i.e. number of first n codon (from st to end, for example 1 to 6 codon from left)
  ve = paste("P", v, ".x", sep = "")
  df$sum = apply(df[,ve], 1, sum)
  df = df[df$sum > 0, ]
  df$sum = df$sum/sum(df$sum)
  df$frequency = df$frequency/sum(df$frequency)
  df$ratio = log2(df$sum/df$frequency)
 

  df = df[order(df$AA, -df$CAI),]
 
  ## -------------------------
  # draw dotplot
  
  setEPS()
  if (out == "eps")   {
    postscript(paste("codon_ICS/", sam, " P", st, "-", ed, " ICS dot.eps", sep=""), pointsize =15)
  }else{
    jpeg(filename = paste("codon_ICS/", sam, " P", st, "-", ed, " ICS dot.jpeg", sep=""),
       width = 600, height = 1200, units = "px", quality = 100)   }
  
  par(mfrow = c(1,1))
  
  cor = cor.test(df$CAI[df$CAI<1], df$ratio[df$CAI<1])
  plot(df$CAI, df$ratio, main = sam, xlab = "CAI", ylab = "normalized codon frequency (log2)", pch = 16, col = "dodgerblue3")
  l = lm(df$ratio~df$CAI)
  abline(l, col = "red", lty = 3, lwd = 3)
  text(max(df$CAI)-0.2, max(df$ratio)-0.5, paste(" r =", round(cor$estimate, 3), "\n n =", 
                                                       cor$parameter, "\n P=", formatC(cor$p.value, digits = 3)))
  
  dev.off()
  
  ##-------------------
  ## draw multiple tiles 
  
  if (out == "eps")   {
    postscript(paste("codon_ICS/", sam, " P", st, "-", ed, " ICS atio codon.eps", sep=""), pointsize =5,
               width = 4, height = 8)
  }else{
    jpeg(filename = paste("codon_ICS/", sam, " P", st, "-", ed, " ICS ratio codon.jpeg", sep=""), 
       width = 600, height = 1200, units = "px", quality = 100)      }
  
  par(mar=c(4,2,2,2))
  par(mfrow=c(6,3),tcl=-0.5,mai=c(0.4,0.2,0.1,0.2))  #mar set the margin from bottom, left, top, right
  
  aa <- c("A","C","D","E","F","G","H","I","K","L","N","P","Q","R","S","T","V","Y")
  
  for (i in aa) {
 
    n = nrow(df[df$AA == i,])
    n = ifelse(n = 2, 2.5, n) 
    min = min(df$ratio[df$AA == i])
    min = ifelse(min>0, -0.2, min-0.2)
    max = max(df$ratio[df$AA == i])
    max = ifelse(max>0, max+0.2, 0.1)
    
    mp<-barplot(t(df[df$AA == i, "ratio"]),ann=FALSE, names.arg = df$codon[df$AA == i], 
                cex.names = 1.2, main = i, las = 2, #axes = F,
                ylim = c(min, max), 
                col = "grey",xlim = c(0,n), width=c(0.5,0.5,0.5), space = 0.8)
    par(new = T)
    abline(h = 0, col = "cyan", lwd = 2, lty = 2)
    par(new = T)
    plot(mp,df$CAI[df$AA == i],type="b",lwd = 2,col="red", xlim = c(0,n), ylim = c(0,1.2), ann=FALSE, axes = F)
    axis(4,at=seq(0,1.2,0.2))
    
    box()   
  }
  dev.off() 
  


  ### ------------------------------
  #report relative changes to control in heatmap 
  ### ------------------------------
  library(pheatmap)
  jpeg(filename = paste(sam, "/codon relative frequency.jpg", sep = ""),
       width = 600, height = 1200, units = "px", quality = 100)
  
  pheatmap(mat, clustering_distance_rows = 'euclidean',  #option is correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'
           border_color = NA, cluster_rows = T, annotation_names_row = T, annotation_names_col = T, show_rownames = T, cluster_cols = F, silent = F,
           annotation_row = NA, fontsize_row = 17, fontsize_col = 22, fontsize = 14)
 
  dev.off() 
  
  # report codon frequency changes 
  jpeg(filename = paste(sam, "/codon frequency.jpg", sep = ""),
       width = 600, height = 1200, units = "px", quality = 100)
  
  pheatmap(m1, clustering_distance_rows = 'euclidean',  #option is correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'
           border_color = NA, cluster_rows = T, annotation_names_row = T, annotation_names_col = T, show_rownames = T, cluster_cols = F, silent = F,
           annotation_row = NA, fontsize_row = 17, fontsize_col = 22, fontsize = 14)
  dev.off() 
  
  # report random codon frequency changes 
  jpeg(filename = paste(sam, "/codon random frequency.jpg", sep = ""),
       width = 600, height = 1200, units = "px", quality = 100)
  
  pheatmap(m2, clustering_distance_rows = 'euclidean',  #option is correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'
           border_color = NA, cluster_rows = T, annotation_names_row = T, annotation_names_col = T, show_rownames = T, cluster_cols = F, silent = F,
           annotation_row = NA, fontsize_row = 17, fontsize_col = 22, fontsize = 14)
  dev.off() 
  
}
  

####==========================================================
## diCodon heatmap
####==========================================================

pas.diCod = function(sam, out="jpeg")  {
  
  dc = read.table(paste(sam, "/PAS_dicodons.txt", sep = ""), sep = "\t", fill = T)
  rpm = read.table(paste(sam, "/CDS_UTR_PolyA_counts.txt", sep = ""), sep = "\t", header = T)
  rpm = rpm[,1:3]
  dc = merge(rpm, dc, by.x = "Gene", by.y = "V1")
  dc$ratio = dc$V2*100/dc$UTR_counts
  df = dc[dc$V3 >= 5 & dc$ratio > 0.01 & dc$ratio < 100, ]
  #make a table to compare the composition of codons 
  x = as.data.frame(table(df$V4))
  names(x) = c("codon", "P1")
  column = 5:12
  for (i in column)  {
    cl = paste("V", i, sep = "")
    cl2 = paste("P", i-3, sep = "")
    d1 = as.data.frame(table(df[[cl]]))
    names(d1) = c("codon", cl2)
    x = merge(x, d1, by.x = "codon", by.y = "codon")
  }
  x$All = apply(x[,2:10], 1, sum)
  dx = sweep(x[,2:11],2,colSums(x[,2:11]),`/`)
  dx$codon = x$codon
  
  ## draw control
  ck = read.table(paste(sam, "/PAS_random_dicodons.txt", sep = ""), sep = "\t")
  #ck = ck[ck$V1 %in% df$Gene,]
  y = as.data.frame(table(ck$V2))
  names(y) = c("codon", "P1")
  column = 3:10
  for (i in column)  {
    cl = paste("V", i, sep = "")
    cl2 = paste("P", i-1, sep = "")
    d1 = as.data.frame(table(ck[[cl]]))
    names(d1) = c("codon", cl2)
    y = merge(y, d1, by.x = "codon", by.y = "codon")
  }
  y$All = apply(y[,2:10], 1, sum)
  y = y[y$codon %in% x$codon,]
  dy = sweep(y[,2:11],2,colSums(y[,2:11]),`/`)
  dy$codon = y$codon
  dy = dy[dy$codon %in% dx$codon,]
  
  m1 = as.matrix(dx[,1:10]) 
  row.names(m1) = dx$codon
  m2 = as.matrix(dy[,1:10])
  row.names(m2) = dy$codon
  mat = m1/m2
  #this two lines is used in case Na, NaN or Inf occur
  mat[is.na(mat)] = 0
  mat[!is.finite(mat)] = max(mat[is.finite(mat)])
  
  d = as.data.frame(mat)
  d = d[order(-d$All),]
  d = d[d$All > 0,]
  write.csv(d, paste(sam, "/diCodon frequency.csv", sep = ""))
  d$name = rownames(d)
  d$order = 1:nrow(d)
  par(mfrow=c(1,1))   
  
  setEPS()
  if (out == "eps")   {
    postscript(paste(sam, "/dicodon frequency.eps", sep=""), pointsize =15, width = 10, height = 5)
  }else{
    jpeg(filename = paste(sam, "/dicodon frequency.jpg", sep=""),
       width = 500, height = 500, units = "px", quality = 100)}
  plot(d$order, d$All, xaxt = "n", xlab = "di-codons", ylab = "relative frequency over background ",
       main = "di-codon frequency", pch = 15, col = "grey40")
  # add label to first a few data
  x = max(d$order)/3
  y = floor(max(d$All)*0.8)
  y = y:(y-4)
  with(d[1:5,],text(x,y,name))
  
  dev.off()
  
}
 



setwd("~/Desktop/ZP_polyA/mouse_SRP063051")
# mouse
sample = c("SRR2225336", "SRR2225337", "SRR2225338", "SRR2225339",
           "SRR2225340", "SRR2225341", "SRR2225342", "SRR2225343", 
           "SRR2225344", "SRR2225345", "SRR2225346", "SRR2225347", 
           "SRR2225348", "SRR2225349", "SRR2225350", "SRR2225351")


dir.create("codon_ICS")
for (s in sample)  {
  PAS.line(s)
  pas.heatmap(s, 6, 11, "eps")
  }









