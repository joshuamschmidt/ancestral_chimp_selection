# data
#---- Notes ----#
# pryr is used because it is a way of wrapping a call to plot as a function.
# helps with complicated plot layouts.
# will send on the appropriate colours (hex) for the different lineages.
# could try thinning the PBS and Diff plots by randomly sampling positions with 
# -log10(p) < 1.
# still want to check all coordinates...
library(data.table)
library(pryr)
cd4regionGenes <- readRDS("cd4regionGenes.rds")
cd4.pclr <- readRDS("cd4.windows.rds")
windows <- readRDS("windows.rds")
cd4.dafs <- readRDS("cd4.dafs.rds")
cd4.pbsnj  <- readRDS("cd4.pbsnj.rds")
nonsyn <- readRDS("cd4.nonsyn.rds")


setkey(cd4regionGenes,chr,start,end)
setkey(cd4.pclr,Chr,Physpos)
setkey(cd4.dafs,chr,start,end)
setkey(cd4.pbsnj,chr,start,end)


# updating gene tracks...
ens <- fread("ens_genes_chr12_6799607-7332026.txt",select = c("name","chrom","strand","txStart","txEnd","name2"))
ens_2_name <- fread("ens_2_gene-name.txt")
ens <- ens_2_name[ens,on=.(name)]
ens <- ens[!value %in% NA]
ens[,strand:=ifelse(strand=="+",1,-1)]
ens <- ens[,.(chrom,start=txStart,stop=txEnd,gene=value,score=0,strand)]
ens <- ens[gene %in% genes$gene]

xmin <- 6875000
xmax <- 7100000

ens[,new.start:=ifelse(start < xmin, xmin, start)]
ens[,new.stop:=ifelse(stop > xmax, xmax, stop)]
# remove some gene names for plotting clarity - list them in the figure legend or methods?
thin <- c("ING4","ZNF384","PIANP","PTMS","LAG3","GPR162","P3H3","GNB3","TPI1","C12orf57","LRRC23")

# plot without p-value cutoffs; annotated AA
# 3p-clr plot

# make the gene layout plot
# xmin <- cd4.dafs[,min(start)]
# xmax <- cd4.dafs[,max(end)]

genes.pryr %<a-% {
  plot(0, type="n", xlab="", ylab="", xlim=c(xmin, xmax), ylim=c(-5, 5),axes=FALSE)
  abline(h = 0)
  for (g in ens[,gene]){
    if(ens[gene==g]$strand==1){
      rect(ens[gene==g]$new.start, 0, ens[gene==g]$new.stop, 3.5,
           col="grey", border=par("fg"), lty=NULL, lwd=par("lwd"), xpd=FALSE)
      if(!g %in% thin){
        text(x=ens[gene==g]$new.start + ((ens[gene==g]$new.stop-ens[gene==g]$new.start)/2),
             y=4.5,
             labels=g)
      }
    }
    if(ens[gene==g]$strand== -1){
      rect(ens[gene==g]$new.start, 0, ens[gene==g]$new.stop, -3.5,
           col="grey", border=par("fg"), lty=NULL, lwd=par("lwd"), xpd=FALSE)
      if(!g %in% thin){
        text(x=ens[gene==g]$new.start + ((ens[gene==g]$new.stop-ens[gene==g]$new.start)/2),
           y= -4.5,
           labels = g)
      }
    }
  }
  # #patch in two labels
  # text(x=cd4.pclr[,min(Physpos)]-3e3,
  #      y= -4.5,
  #      labels="ACRBP")
  # text(x=cd4.pclr[,max(Physpos)]+4e3,
  #      y= -4.5,
  #      labels="C1R")
}
# make the 3p-clr scores plot
pclr.pryr %<a-% {
  plot(windows[Chr==12 & Physpos >= xmin & Physpos <= xmax][order(Physpos)]$Physpos,
       -log10(windows[Chr==12 & Physpos >= xmin & Physpos <= xmax][order(Physpos)]$anc.p),
       type='b',
       pch=19,
       col='brown',
       lwd=2,
       xaxt = "n",
       yaxt = "n",
       ylim=c(0,5),
       xlim=c(xmin,xmax),
       ylab = "",
       cex.lab=1.4,
       bty="n"
  )
  # add axis on y
  ytick<-seq(0, 5, by=2.5)
  axis(2, at = ytick,las=1,cex=1)
  legend("topright",
         title=NULL,
         legend="central - eastern\nancestor 3P-CLR",
         col="brown",
         pch = 19,
         pt.bg = 'white',
         lty = 1,
         lwd=2,
         bty="n")
  
}

clade.DAF.pryr %<a-% {
  cd4.dafs[order(start)][cladeDiffP < 0.3][,plot(start,-log10(cladeDiffP),type='b',col='chocolate',ylim=c(0,3),xaxt = "n",
                                                 yaxt = "n",
                                                 ylab = "",
                                                 cex.lab=1.4,
                                                 lwd=2,
                                                 bty="n",
                                                 pch=19,
                                                 xlim=c(xmin,xmax))]
  legend("topright",
         title=NULL,
         legend="internal-DAF",
         col="chocolate",
         pch = 19,
         lty = 1,
         lwd=2,
         bty="n")
  ytick<-seq(0, 5, by=2.5)
  axis(2, at = ytick,las=1,cex=1)
  xtick<-seq(round(xmin/1e6,1), round(xmax/1e6,1), by=0.05)
  #axis(1, at = xtick,cex=1)
  axis(side=1,at=xtick,tcl=0.4,lwd.ticks=3,mgp=c(0,0.5,0),outer=F,pos = 2)
  cd4.dafs[start==nonsyn[1]$start][,points(start,-log10(cladeDiffP),pch=19,col='purple',cex=2)]
  cd4.dafs[start==nonsyn[3]$start][,points(start,-log10(cladeDiffP),pch=19,col='seagreen',cex=2)]
  cd4.dafs[start==nonsyn[2]$start][,points(start,-log10(cladeDiffP),pch=19,col='pink',cex=2)]
  #c
}

# central pbsnj
centralpbsnj.pryr %<a-% {
  cd4.pbsnj[order(start)][c.p < 0.3][,plot(start,-log10(c.p),type='b',col='green',ylim=c(0,5),xaxt = "n",
                                           yaxt = "n",
                                           ylab = "",
                                           cex.lab=1.4,
                                           lwd=2,
                                           bty="n",
                                           pch=19,
                                           xlim=c(xmin,xmax))]
  legend("topright",
         title=NULL,
         legend="central PBSnj",
         col="green",
         pch = 19,
         lty = 1,
         lwd=2,
         bty="n")
  ytick<-seq(0, 5, by=2.5)
  axis(2, at = ytick,las=1,cex=1)
  cd4.pbsnj[start==nonsyn[1]$start][,points(start,-log10(c.p),pch=19,col='purple',cex=2)]
  cd4.pbsnj[start==nonsyn[3]$start][,points(start,-log10(c.p),pch=19,col='seagreen',cex=2)]
  cd4.pbsnj[start==nonsyn[2]$start][,points(start,-log10(c.p),pch=19,col='pink',cex=2)]
}
# eastern pbsnj
easternpbsnj.pryr %<a-% {
  cd4.pbsnj[order(start)][e.p < 0.3][,plot(start,-log10(e.p),type='b',col='orange',ylim=c(0,3), xlab="Position on Chr12",
                                           yaxt = "n",
                                           ylab = "",
                                           cex.lab=1.4,
                                           lwd=2,
                                           bty="n",
                                           pch=19,
                                           xlim=c(xmin,xmax))]
  legend("topright",
         title=NULL,
         legend="eastern PBSnj",
         col="orange",
         pch = 19,
         lty = 1,
         lwd=2,
         bty="n")
  ytick<-seq(0, 5, by=2.5)
  axis(2, at = ytick,las=1,cex=1)
  cd4.pbsnj[start==nonsyn[1]$start][,points(start,-log10(e.p),pch=19,col='purple',cex=2)]
  cd4.pbsnj[start==nonsyn[3]$start][,points(start,-log10(e.p),pch=19,col='seagreen',cex=2)]
  cd4.pbsnj[start==nonsyn[2]$start][,points(start,-log10(e.p),pch=19,col='pink',cex=2)]
}

# internal pbsnj
intpbsnj.pryr %<a-% {
  cd4.pbsnj[order(start)][i.p < 0.3][,plot(start,-log10(i.p),type='b',col='brown',ylim=c(0,3),xaxt = "n",
                                           yaxt = "n",
                                           ylab = "",
                                           cex.lab=1.4,
                                           lwd=2,
                                           pch=19,
                                           bty="n",
                                           xlim=c(xmin,xmax))]
  legend("topright",
         title=NULL,
         legend="internal PBSnj",
         col="brown",
         pch = 19,
         lty = 1,
         lwd=2,
         bty="n")
  ytick<-seq(0, 5, by=2.5)
  axis(2, at = ytick,las=1,cex=1)
  cd4.pbsnj[start==nonsyn[1]$start][,points(start,-log10(i.p),pch=19,col='purple',cex=2)]
  cd4.pbsnj[start==nonsyn[3]$start][,points(start,-log10(i.p),pch=19,col='seagreen',cex=2)]
  cd4.pbsnj[start==nonsyn[2]$start][,points(start,-log10(i.p),pch=19,col='pink',cex=2)]
}

pdf(file="CD4_region.sweep.stats_nopvaluelines_annotatedAA.pdf", pointsize=5, width=16.9/2.54, height=8/2.54)
par(oma=c(3,1,1,1),mar = c(0.5,4.5,0,0))
layout(matrix(c(1,1,2,2,3,3,4,4,5,5,5,6,6), 13, 1, byrow = TRUE))
genes.pryr
pclr.pryr
intpbsnj.pryr
clade.DAF.pryr
centralpbsnj.pryr
easternpbsnj.pryr
mtext("-log10(p-value)", outer = TRUE, cex = 1,side = 2,line = -1)
mtext("Position along Chr12", outer = TRUE, cex = 1,side = 1,line = 1.5)
dev.off()

