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


# stat test?
# need the 3pclr pvalues or LL
ll <- fread("sed -e '/Chr/d' /Users/joshuaschmidt/Projects/ancestral_chimp/all.windows.3pclr.ce.nc.txt",select = c("V1","V2","V9"),col.names = c("chr","pos","LL"))
setkey(ll,chr,pos)
ll <- ll[back,on=.(chr,pos)][!Internal.n %in% NA]
model <- lm(Internal.n ~ LL , data = ll)
summary(model)


# Call:
#   lm(formula = Internal.n ~ LL, data = ll)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.29463 -0.06789 -0.02794  0.04016  0.75704 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 1.072e-01  1.221e-04   877.9   <2e-16 ***
#   LL          3.206e-04  2.981e-06   107.5   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.09568 on 824808 degrees of freedom
# Multiple R-squared:  0.01383,	Adjusted R-squared:  0.01383 
# F-statistic: 1.157e+04 on 1 and 824808 DF,  p-value: < 2.2e-16