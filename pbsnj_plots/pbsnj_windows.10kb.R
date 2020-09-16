library(data.table)
setwd("~/Projects/ancestral_chimpanzee_selection")
# combine all
allPBSnj <- data.table()
for ( c in c(1,"2A","2B",3:22)){
  load(paste0("chr",c,"_WindowSize_5KB_pbs4_bahtia_nj.sliding.per.snp.Rdata"))
  allPBSnj <- rbindlist(list(allPBSnj,PBSNJ))
  rm(PBSNJ)
}

allPBSnj[,Internal.s:=Internal/(1+Internal+Western + Nigeria + Central + Eastern)]
allPBSnj[,Internal.n:=(Internal.s-min(Internal.s))/(max(Internal.s)-min(Internal.s))]
# candidates
cands <- fread("pbsnj_plots/aligned_3pclr_out_genome_wide_0.005percentsubsetCANDSNP.txt",col.names = c("chr","pos"))
cands <- allPBSnj[cands,on=.(chr,pos),roll="nearest"]


cands2 <- fread("pbsnj_plots/aligned_3pclr_out_genome_wide_0.0005percentsubsetCANDSNP.txt",col.names = c("chr","pos"))
cands2 <- allPBSnj[cands2,on=.(chr,pos),roll="nearest"]

cands3 <- fread("pbsnj_plots/aligned_3pclr_out_genome_wide_0.1percentsubsetCANDSNP.txt",col.names = c("chr","pos"))
cands3 <- allPBSnj[cands3,on=.(chr,pos),roll="nearest"]

# background
back <- fread("pbsnj_plots/all.windows.3pclr.ce.ncBACKGROUNDSNPint.txt",col.names = c("chr","pos"))
back <- allPBSnj[back,on=.(chr,pos),roll="nearest"]
back[,median(Internal.n,na.rm = T)]
back[,mean(Internal.n,na.rm = T)]

library(vioplot)
pdf(file="vioplots-PBSnjInt-5kb.windows-3pclr.Candidates.pdf",
    pointsize = 10,
    width = 5,
    height = 5)
vioplot(back[!Internal.n %in% NA]$Internal.n,
        cands[!Internal.n %in% NA]$Internal.n,
        cands3[!Internal.n %in% NA]$Internal.n,
        cands2[!Internal.n %in% NA]$Internal.n,
        names=c("background","0.5%","0.1%","0.05%"),
        col=c("#659e78","#dec4a3","#dc9e76","#da735a"),
        horizontal = T,xlab="PBSnjInt")
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