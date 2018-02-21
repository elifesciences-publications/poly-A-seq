# this code is modified from N crassa and add special conditions for mouse

# define the position where dicodon frequency is started to count, and end at where

start = 1
end = 5
# 1 to 5 means position of dicodons from about -30 to -10 region. which will have 10 codons upstream of cleavage
# 1 to 5 approximately represent the region where PAS motif located.


pas.diCod = function(sam, out="jpeg")  {
  
  dc = read.table(paste(sam, "/PAS_dicodons.txt", sep = ""), sep = "\t", fill = T, stringsAsFactors=FALSE)
  rpm = read.table(paste(sam, "/CDS_UTR_PolyA_counts.txt", sep = ""), sep = "\t", header = T, stringsAsFactors=FALSE)
  rpm = rpm[,1:3]
  dc = merge(rpm, dc, by.x = "Gene", by.y = "V1")
  dc$ratio = dc$V2*100/dc$UTR_counts
  dc$rpkm = dc$UTR_counts*1000000/sum(dc$UTR_counts)
  # remove the gene has very low expression (rpkm < 1)
  dc = dc[dc$rpkm>1,]
  
  # first set those peak that only have 1 reads as a unique group
  d = dc[dc$V3 == 1, ]
  x = as.data.frame(table(d$V4), stringsAsFactors=FALSE)
  names(x) = c("dicodon", "P1")
  column = (start+4):(end+3)
  for (i in column)  {
    cl = paste("V", i, sep = "")
    cl2 = paste("P", i-3, sep = "")
    d1 = as.data.frame(table(d[[cl]]), stringsAsFactors=FALSE)
    names(d1) = c("dicodon", cl2)
    x = merge(x, d1, by.x = "dicodon", by.y = "dicodon", all = T)
  }
  x[is.na(x)] = 0
  x[["P=1"]] = apply(x[,(start+1):(end+1)], 1, sum)
  df = x[,c(1,(end+2))]
  
  
  ##--------------------------------
  ##analyze the rest 
  ## -------------------------------
  d = dc[dc$V3 > 1, ]
  
  # set a quantile to see the enrichment of different ratio of motif.
  qu = quantile(d$ratio)
  
  
  for (j in 1:4)  {
    
    # first set those peak that only have 1 reads as a unique group
    d2 = d[d$ratio >= qu[j] & d$ratio < qu[j+1], ]
    x = as.data.frame(table(d2$V4), stringsAsFactors=FALSE)
    names(x) = c("dicodon", "P1")
    column = (start+4):(end+3)
    for (i in column)  {
      cl = paste("V", i, sep = "")
      cl2 = paste("P", i-3, sep = "")
      d1 = as.data.frame(table(d[[cl]]), stringsAsFactors=FALSE)
      names(d1) = c("dicodon", cl2)
      x = merge(x, d1, by.x = "dicodon", by.y = "dicodon", all = T)
    }
    x[is.na(x)] = 0
    x[[paste("Q", j, sep = "")]] = apply(x[,(start+1):(end+1)], 1, sum)
    d2 = x[,c(1,(end+2))]
    
    df = merge(df, d2, by.x = "dicodon", by.y = "dicodon", all = T)
    df[is.na(df)] = 0
    
  }
  
  df$Qall = apply(df[,3:6], 1, sum)
  df1 = sweep(df[,2:7],2,colSums(df[,2:7]),`/`)
  df1$dicodon = df$dicodon
  
  
  # read the reference 
  gc = read.table(paste(sam, "/dicodon_genome_frequency.txt", sep = ""), sep = "\t", header = T, stringsAsFactors=FALSE)
  gc[is.na(gc)] = 0
  gc$cod1 = substr(gc$dicodon, 1, 3)
  gc$cod2 = substr(gc$dicodon, 4, 6)
  
  cod = read.table("ref/codon_frequency.txt", sep = "\t", stringsAsFactors = F)
  cod$V3 = cod$V3/sum(cod$V3)
  cod = cod[,c(1,3)]
  
  gc = merge(gc, cod, by.x = "cod1", by.y = "V1")
  names(gc)[7] = "cai1"
  gc = merge(gc, cod, by.x = "cod2", by.y = "V1")
  names(gc)[8] = "cai2"
  gc$theory = gc$cai1*gc$cai2
  gc$genome = gc$genome/sum(gc$genome)
  gc$ORF_cleavage = gc$ORF_cleavage/sum(gc$ORF_cleavage)
  gc$none = gc$none/sum(gc$none)
  
  gc = gc[, c("dicodon","theory", "genome", "ORF_cleavage", "none")]
  #sum(gc$theory) = this does not comes to 0 became not all combination exists
  
  
  
  df = merge(gc, df1, by.x = "dicodon", by.y = "dicodon")
  
  df[-1] = sweep(df[-1],1,df$genome,`/`)
  
  
  df$col = "black"
  
  motif = read.table(paste(sam, "/polyA_UTR3_-30--10_sig_motifs.txt", sep = ""), sep = "\t", header = T,, stringsAsFactors=FALSE)
  motif = motif[1:15,]
  motif$motif = gsub("U", "T", motif$motif)
  df$col[df$dicodon %in% motif$motif] = "red"
  df = df[-2]
  
  for (k in 2:10)   {
    dt = df[,c(1,k,11)]
    dt = dt[dt[,2]>0,]
    dt = dt[order(-dt[,2]),]
    dt[,2] = log2(dt[,2])
    dt$order = 1:nrow(dt)
    
    setEPS()
    if (out == "eps")   {
      postscript(paste(sam, "/dicodon frequency ", names(df)[k], ".eps", sep=""), pointsize =10, width = 4, height = 4)
    }else{
      jpeg(filename = paste(sam, "/dicodon frequency ", names(df)[k], ".eps", sep=""),
           width = 500, height = 500, units = "px", quality = 100)}
    plot(dt$order, dt[,2], xaxt = "n", xlab = "di-codons", ylab = "relative frequency over genome ",
         main = paste(sam, "dicodon frequency", names(df)[k]), pch = 16, col = "grey40", cex = 0.25)
    #par(new=T)
    #points(dt$order[dt$col=="red"], dt[,2][dt$col=="red"], col = "red", cex = 1, pch = 1)
    
    par(new=T)
    points(dt$order[dt$dicodon=="AATAAA"], dt[,2][dt$dicodon=="AATAAA"], col = "red", cex = 1, pch = 1)
    #text(dt$order[dt$dicodon=="AATAAA"]+400, dt[,2][dt$dicodon=="AATAAA"], 
         #col = "red", labels = "AAUAAA")
    
    par(new=T)
    points(dt$order[dt$dicodon=="ATTAAA"], dt[,2][dt$dicodon=="ATTAAA"], col = "green", cex = 1, pch = 1)
    #text(dt$order[dt$dicodon=="ATTAAA"]+400, dt[,2][dt$dicodon=="ATTAAA"], 
         #col = "green", labels = "AUUAAA")
 
    dev.off()    
  }
}



setwd("<your working directory>")  

#sample name should be the same as your used in previous step, from step1 to now
sample = c("SRR2225336", "SRR2225337", "SRR2225338", "SRR2225339",
           "SRR2225340", "SRR2225341", "SRR2225342", "SRR2225343", 
           "SRR2225344", "SRR2225345", "SRR2225346", "SRR2225347", 
           "SRR2225348", "SRR2225349", "SRR2225350", "SRR2225351")

for (s in sample)  {
  pas.diCod(s, "eps")
}
