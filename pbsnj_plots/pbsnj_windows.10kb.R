library(data.table)
library(vioplot)
setwd("~/Projects/ancestral_chimp_selection")
# combine all
allPBSnj <- data.table()
for ( c in c(1,"2A","2B",3:22)){
  load(paste0("pbsnj_plots/pbsnj_windows/chr",c,"_WindowSize_5KB_pbs4_bahtia_nj.sliding.per.snp.Rdata"))
  allPBSnj <- rbindlist(list(allPBSnj,PBSNJ))
  rm(PBSNJ)
}

allPBSnj[,Internal.s:=Internal/(1+Internal+Western + Nigeria + Central + Eastern)]
allPBSnj[,Internal.n:=(Internal.s-min(Internal.s))/(max(Internal.s)-min(Internal.s))]
# candidates
cands_5 <- fread("pbsnj_plots/aligned_3pclr_out_genome_wide_0.005percentsubsetCANDSNP.txt",col.names = c("chr","pos"))
cands_5 <- allPBSnj[cands_5,on=.(chr,pos),roll="nearest"]


cands_05 <- fread("pbsnj_plots/aligned_3pclr_out_genome_wide_0.0005percentsubsetCANDSNP.txt",col.names = c("chr","pos"))
cands_05 <- allPBSnj[cands_05,on=.(chr,pos),roll="nearest"]

cands_1 <- fread("pbsnj_plots/aligned_3pclr_out_genome_wide_0.1percentsubsetCANDSNP.txt",col.names = c("chr","pos"))
cands_1 <- allPBSnj[cands_1,on=.(chr,pos),roll="nearest"]

# background
back <- fread("pbsnj_plots/all.windows.3pclr.ce.ncBACKGROUNDSNPint.txt",col.names = c("chr","pos"))
back <- allPBSnj[back,on=.(chr,pos),roll="nearest"]
# back[,median(Internal.n,na.rm = T)]
# back[,mean(Internal.n,na.rm = T)]
opar <- par()
pdf(file="pbsnj_plots/vioplots-PBSnjInt-5kb.windows-3pclr.Candidates.pdf",
    width = 1.75,
    height = 1.75,
    pointsize = 6,fonts = "Helvetica")
par(mfrow=c(1,1))
par(las=1)
par(mgp=c(1.8,0.6,0))
par(mar = c(3, 2.5, 0.1, 0.1))
vioplot(back[!Internal.n %in% NA]$Internal.n,
        cands_5[!Internal.n %in% NA]$Internal.n,
        cands_1[!Internal.n %in% NA]$Internal.n,
        cands_05[!Internal.n %in% NA]$Internal.n,
        names=c("BG","0.5","0.1","0.05"),
        col=c("#659e78","#dec4a3","#dc9e76","#da735a"),
        axes=FALSE,
        ann=FALSE,
        horizontal = T,
        ylim=c(0,1))
mtext(side = 1.2,text = "PBSnjINT",line = 1.5,cex = 1.25,at = 0.5)
dev.off()
# # also try ggplot style
# library(ggplot2)
# longDT <- data.table()
# 
# longDT <- rbindlist(list(longDT,data.table(class="background",PBSnjInt=back[!Internal.n %in% NA]$Internal.n)))
# longDT <- rbindlist(list(longDT,data.table(class="0.5%",PBSnjInt=cands[!Internal.n %in% NA]$Internal.n)))
# longDT <- rbindlist(list(longDT,data.table(class="0.1%",PBSnjInt=cands2[!Internal.n %in% NA]$Internal.n)))
# longDT <- rbindlist(list(longDT,data.table(class="0.05%",PBSnjInt=cands3[!Internal.n %in% NA]$Internal.n)))
# 
# 
# median.quartile <- function(x){
#   out <- quantile(x, probs = c(0.25,0.5,0.75))
#   names(out) <- c("ymin","y","ymax")
#   return(out) 
# }
# 
# 
# p <- ggplot(longDT, aes(PBSnjInt, factor(class)))
# p + geom_violin() + stat_summary(fun.y=median.quartile,geom='point')


t.test(cands[!Internal.n %in% NA]$Internal.n,back[!Internal.n %in% NA]$Internal.n)