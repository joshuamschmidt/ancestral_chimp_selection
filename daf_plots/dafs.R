library(data.table)
# dafs
load("~/Box/ancestral_chimp/chimp.bonobo.dafs.conditoned.human.Rdata")
# dafs.human <- dafs.human[e.daf+c.daf+n.daf+w.daf!=0]
# dafs.human <- dafs.human[e.daf+c.daf+n.daf+w.daf!=4]
dafs.human <- dafs.human[n.daf!=0 & n.daf!=1]
setkey(dafs.human,chr,start,end)


# windows
all <- fread("~/Projects/ancestral_chimp/all.windows.3pclr.ce.nc.txt",select = c("V1","V2","V4","V5"),col.names = c("chr","focal","start","end"))
all[,`:=`(focal=as.integer(focal),start=as.integer(start),end=as.integer(end))]
all <- all[chr!="Chr"]
all[,id:=1:.N]
setkey(all,chr,start,end)


# # randomly sample 200000 windows
# sample.win <- all[all[,sample(.N,size = 200e3,replace = F)]]
# setkey(sample.win,chr,start,end)
# sample.daf <- unique(foverlaps(dafs.human,sample.win,nomatch = 0)[,.(chr,start=i.start,end=i.end,e.daf,c.daf)])


# DT hist function

dtHist <- function(dt,stat,maxC) {
  dt[,.(count=round(get(stat)*maxC))][,.N,by=count][order(count),.(count,frequency=N/sum(N))]
}



# candidate 0.5% windows
cands_0.5 <- fread("~/Projects/ancestral_chimp/aligned_3pclr_out_genome_wide_0.005percentsubsetCANDSNP.txt",col.names = c("chr","focal"))
cands_0.5.win <- all[cands_0.5,on=.(chr,focal)]
setkey(cands_0.5.win,chr,start,end)
cands_0.5.daf <- unique(foverlaps(dafs.human,cands_0.5.win,nomatch = 0)[,.(chr,start=i.start,end=i.end,e.daf,c.daf)])


setwd("~/Projects/ancestral_chimp/candidate.windows.dafs")
pdf(file="candidates.0.5.daf.pdf",height = 8, width=8)
par(mfrow=c(2,1))
# eastern chimp
e1_0.5 <- dtHist(dafs.human,"e.daf",38)
e2_0.5 <- dtHist(cands_0.5.daf,"e.daf",38)

# central chimp
c1_0.5 <- dtHist(dafs.human,"c.daf",36)
c2_0.5 <- dtHist(cands_0.5.daf,"c.daf",36)

# plots
yscale <- round(range(e1_0.5[,range(frequency)],e2_0.5[,range(frequency)],c1_0.5[,range(frequency)],c2_0.5[,range(frequency)]),digits=1)

# eastern plot
barplot(height=e1_0.5$frequency,
        xlab = "SFS",
        ylab="Proportion of sites",
        ylim=yscale,
        names.arg=0:38,
        col=rgb(0,0,1,1/4),
        main="Eastern SFS: background vs anc candidates")
barplot(height=e2_0.5$frequency,
        add = T,
        col=rgb(1,0,0,1/4))
legend("topright",legend = c("background","0.5% candidates"),fill = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)))

# central plot
barplot(height=c1_0.5$frequency,
        xlab = "SFS",
        ylab="Proportion of sites",
        ylim=yscale,
        names.arg=0:36,
        col=rgb(0,0,1,1/4),
        main="Central SFS: background vs anc candidates")
barplot(height=c2_0.5$frequency,
        add = T,
        col=rgb(1,0,0,1/4))
legend("topright",legend = c("background","0.5% candidates"),fill = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)))
dev.off()

# candidate 0.1% windows
cands_0.1 <- fread("~/Projects/ancestral_chimp/aligned_3pclr_out_genome_wide_0.1percentsubsetCANDSNP.txt",col.names = c("chr","focal"))
cands_0.1.win <- all[cands_0.1,on=.(chr,focal)]
setkey(cands_0.1.win,chr,start,end)
cands_0.1.daf <- unique(foverlaps(dafs.human,cands_0.1.win,nomatch = 0)[,.(chr,start=i.start,end=i.end,e.daf,c.daf)])


setwd("~/Projects/ancestral_chimp/candidate.windows.dafs")
pdf(file="candidates.0.1.daf.pdf",height = 8, width=8)
par(mfrow=c(2,1))
# eastern chimp
e1_0.1 <- dtHist(dafs.human,"e.daf",38)
e2_0.1 <- dtHist(cands_0.1.daf,"e.daf",38)

# central chimp
c1_0.1 <- dtHist(dafs.human,"c.daf",36)
c2_0.1 <- dtHist(cands_0.1.daf,"c.daf",36)

# plots
yscale <- round(range(e1_0.1[,range(frequency)],e2_0.1[,range(frequency)],c1_0.1[,range(frequency)],c2_0.1[,range(frequency)]),digits=1)

# eastern plot
barplot(height=e1_0.1$frequency,
        xlab = "SFS",
        ylab="Proportion of sites",
        ylim=yscale,
        names.arg=0:38,
        col=rgb(0,0,1,1/4),
        main="Eastern SFS: background vs anc candidates")
barplot(height=e2_0.1$frequency,
        add = T,
        col=rgb(1,0,0,1/4))
legend("topright",legend = c("background","0.1% candidates"),fill = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)))

# central plot
barplot(height=c1_0.1$frequency,
        xlab = "SFS",
        ylab="Proportion of sites",
        ylim=yscale,
        names.arg=0:36,
        col=rgb(0,0,1,1/4),
        main="Central SFS: background vs anc candidates")
barplot(height=c2_0.1$frequency,
        add = T,
        col=rgb(1,0,0,1/4))
legend("topright",legend = c("background","0.1% candidates"),fill = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)))
dev.off()



# candidate 0.05% windows
cands_0.05 <- fread("~/Projects/ancestral_chimp/aligned_3pclr_out_genome_wide_0.0005percentsubsetCANDSNP.txt",col.names = c("chr","focal"))
cands_0.05.win <- all[cands_0.05,on=.(chr,focal)]
setkey(cands_0.05.win,chr,start,end)
cands_0.05.daf <- unique(foverlaps(dafs.human,cands_0.05.win,nomatch = 0)[,.(chr,start=i.start,end=i.end,e.daf,c.daf)])


setwd("~/Projects/ancestral_chimp/candidate.windows.dafs")
pdf(file="candidates.0.05.daf.pdf",height = 8, width=8)
par(mfrow=c(2,1))
# eastern chimp
e1_0.05 <- dtHist(dafs.human,"e.daf",38)
e2_0.05 <- dtHist(cands_0.05.daf,"e.daf",38)

# central chimp
c1_0.05 <- dtHist(dafs.human,"c.daf",36)
c2_0.05 <- dtHist(cands_0.05.daf,"c.daf",36)

# plots
yscale <- round(range(e1_0.05[,range(frequency)],e2_0.05[,range(frequency)],c1_0.05[,range(frequency)],c2_0.05[,range(frequency)]),digits=1)

# eastern plot
barplot(height=e1_0.05$frequency,
        xlab = "SFS",
        ylab="Proportion of sites",
        ylim=yscale,
        names.arg=0:38,
        col=rgb(0,0,1,1/4),
        main="Eastern SFS: background vs anc candidates")
barplot(height=e2_0.05$frequency,
        add = T,
        col=rgb(1,0,0,1/4))
legend("topright",legend = c("background","0.05% candidates"),fill = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)))

# central plot
barplot(height=c1_0.05$frequency,
        xlab = "SFS",
        ylab="Proportion of sites",
        ylim=yscale,
        names.arg=0:36,
        col=rgb(0,0,1,1/4),
        main="Central SFS: background vs anc candidates")
barplot(height=c2_0.05$frequency,
        add = T,
        col=rgb(1,0,0,1/4))
legend("topright",legend = c("background","0.05% candidates"),fill = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)))
dev.off()

# MWU tests on the FOLDED SFSs (MAF)

etest_05 <- wilcox.test(x=round(dafs.human[,.(e.maf=ifelse(e.daf < 0.5,e.daf,1-e.daf))]$e.maf*38),
                       y=round(cands_0.5.daf[,.(e.maf=ifelse(e.daf < 0.5,e.daf,1-e.daf))]$e.maf*38),
                       alternative = "two.sided")

etest_01 <- wilcox.test(x=round(dafs.human[,.(e.maf=ifelse(e.daf < 0.5,e.daf,1-e.daf))]$e.maf*38),
                       y=round(cands_0.1.daf[,.(e.maf=ifelse(e.daf < 0.5,e.daf,1-e.daf))]$e.maf*38),
                       alternative = "two.sided")

etest_005 <- wilcox.test(x=round(dafs.human[,.(e.maf=ifelse(e.daf < 0.5,e.daf,1-e.daf))]$e.maf*38),
                       y=round(cands_0.05.daf[,.(e.maf=ifelse(e.daf < 0.5,e.daf,1-e.daf))]$e.maf*38),
                       alternative = "two.sided")


ctest_05 <- wilcox.test(x=round(dafs.human[,.(c.maf=ifelse(c.daf < 0.5,c.daf,1-c.daf))]$c.maf*36),
                        y=round(cands_0.5.daf[,.(c.maf=ifelse(c.daf < 0.5,c.daf,1-c.daf))]$c.maf*36),
                        alternative = "two.sided")

ctest_01 <- wilcox.test(x=round(dafs.human[,.(c.maf=ifelse(c.daf < 0.5,c.daf,1-c.daf))]$c.maf*36),
                        y=round(cands_0.1.daf[,.(c.maf=ifelse(c.daf < 0.5,c.daf,1-c.daf))]$c.maf*36),
                        alternative = "two.sided")

ctest_005 <- wilcox.test(x=round(dafs.human[,.(c.maf=ifelse(c.daf < 0.5,c.daf,1-c.daf))]$c.maf*36),
                         y=round(cands_0.05.daf[,.(c.maf=ifelse(c.daf < 0.5,c.daf,1-c.daf))]$c.maf*36),
                         alternative = "two.sided")

# plot the scaled SFSs
escaled05 <- e2_0.5$density/e1_0.5$density
escaled01 <- e2_0.1$density/e1_0.1$density
escaled005 <- e2_0.05$density/e1_0.05$density
cscaled05 <- c2_0.5$density/c1_0.5$density
cscaled01 <- c2_0.1$density/c1_0.1$density
cscaled005 <- c2_0.05$density/c1_0.05$density

yscale <- range(c(cscaled005,escaled005,cscaled01,escaled01,cscaled05,escaled05))
par(mfrow=c(2,1))
plot(x=0:38,y=escaled005,type='b',ylim=c(yscale[1],yscale[2]))


# # nigeria chimp
# p1 <- hist(round(sample.daf$n.daf*20),breaks = 20,probability = T,plot = F)
# p2 <- hist(round(cands.daf$n.daf*20),breaks = 20,probability = T,plot = F)
# plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,20),freq = F,ylim=c(0,0.8),main="Nigeria DAF: background vs anc candidates",xlab="DAF")  # first histogram
# plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,20), add=T,freq = F)  # second
# legend("topright",legend = c("background","0.5% candidates"),fill = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)))
# 
# 
# # western chimp
# p1 <- hist(round(sample.daf$w.daf*22),breaks = 22,probability = T,plot = F)
# p2 <- hist(round(cands.daf$w.daf*22),breaks = 22,probability = T,plot = F)
# plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,22),freq = F,ylim=c(0,0.8),main="Western DAF: background vs anc candidates",xlab="DAF")  # first histogram
# plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,22), add=T,freq = F)  # second
# legend("topright",legend = c("background","0.5% candidates"),fill = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)))
dev.off()




