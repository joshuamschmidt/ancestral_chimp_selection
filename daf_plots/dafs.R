
library(data.table)
# dafs
load("/Users/joshuaschmidt/Box/ancestral_chimp/chimp.bonobo.dafs.conditoned.human.Rdata")
dafs.human <- dafs.human[e.daf+c.daf+n.daf+w.daf!=0]
dafs.human <- dafs.human[e.daf+c.daf+n.daf+w.daf!=4]
dafs.human <- dafs.human[n.daf!=0 & n.daf!=1]
setkey(dafs.human,chr,start,end)


# windows
all <- fread("/Users/joshuaschmidt/Projects/ancestral_chimp/all.windows.3pclr.ce.nc.txt",select = c("V1","V2","V4","V5"),col.names = c("chr","focal","start","end"))
all[,`:=`(focal=as.integer(focal),start=as.integer(start),end=as.integer(end))]
all <- all[chr!="Chr"]
all[,id:=1:.N]
setkey(all,chr,start,end)


# randomly sample 50000 windows
sample.win <- all[all[,sample(.N,size = 50000,replace = F)]]
setkey(sample.win,chr,start,end)
sample.daf <- foverlaps(dafs.human,sample.win,nomatch = 0)

# candidate windows
cands <- fread("/Users/joshuaschmidt/Projects/ancestral_chimp/aligned_3pclr_out_genome_wide_0.005percentsubsetCANDSNP.txt",col.names = c("chr","focal"))
cands.win <- all[cands,on=.(chr,focal)]
setkey(cands.win,chr,start,end)
cands.daf <- foverlaps(dafs.human,cands.win,nomatch = 0)


setwd("/Users/joshuaschmidt/Projects/ancestral_chimp/candidate.windows.dafs")
pdf(file="candidates.0.5.daf.pdf",height = 8, width=8)
par(mfrow=c(2,2))
# eastern chimp
p1 <- hist(round(sample.daf$e.daf*38),breaks = 38,probability = T,plot = F)
p2 <- hist(round(cands.daf$e.daf*38),breaks = 38,probability = T,plot = F)
plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,38),freq = F,ylim=c(0,0.8),main="Eastern DAF: background vs anc candidates",xlab="DAF")  # first histogram
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,38), add=T,freq = F)  # second
legend("topright",legend = c("background","0.5% candidates"),fill = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)))
# central chimp
p1 <- hist(round(sample.daf$c.daf*36),breaks = 36,probability = T,plot = F)
p2 <- hist(round(cands.daf$c.daf*36),breaks = 36,probability = T,plot = F)
plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,36),freq = F,ylim=c(0,0.8),main="Central DAF: background vs anc candidates",xlab="DAF")  # first histogram
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,36), add=T,freq = F)  # second
legend("topright",legend = c("background","0.5% candidates"),fill = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)))

# nigeria chimp
p1 <- hist(round(sample.daf$n.daf*20),breaks = 20,probability = T,plot = F)
p2 <- hist(round(cands.daf$n.daf*20),breaks = 20,probability = T,plot = F)
plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,20),freq = F,ylim=c(0,0.8),main="Nigeria DAF: background vs anc candidates",xlab="DAF")  # first histogram
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,20), add=T,freq = F)  # second
legend("topright",legend = c("background","0.5% candidates"),fill = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)))


# western chimp
p1 <- hist(round(sample.daf$w.daf*22),breaks = 22,probability = T,plot = F)
p2 <- hist(round(cands.daf$w.daf*22),breaks = 22,probability = T,plot = F)
plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,22),freq = F,ylim=c(0,0.8),main="Western DAF: background vs anc candidates",xlab="DAF")  # first histogram
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,22), add=T,freq = F)  # second
legend("topright",legend = c("background","0.5% candidates"),fill = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)))
dev.off()




# 0.05% candidates
# candidate windows
cands <- fread("/Users/joshuaschmidt/Projects/ancestral_chimp/aligned_3pclr_out_genome_wide_0.0005percentsubsetCANDSNP.txt",col.names = c("chr","focal"))
cands.win <- all[cands,on=.(chr,focal)]
setkey(cands.win,chr,start,end)
cands.daf <- foverlaps(dafs.human,cands.win,nomatch = 0)


pdf(file="candidates.0.05.daf.pdf",height = 8, width=8)
par(mfrow=c(2,2))
# eastern chimp
p1 <- hist(round(sample.daf$e.daf*38),breaks = 38,probability = T,plot = F)
p2 <- hist(round(cands.daf$e.daf*38),breaks = 38,probability = T,plot = F)
plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,38),freq = F,ylim=c(0,0.8),main="Eastern DAF: background vs anc candidates",xlab="DAF")  # first histogram
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,38), add=T,freq = F)  # second
legend("topright",legend = c("background","0.05% candidates"),fill = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)))
# central chimp
p1 <- hist(round(sample.daf$c.daf*36),breaks = 36,probability = T,plot = F)
p2 <- hist(round(cands.daf$c.daf*36),breaks = 36,probability = T,plot = F)
plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,36),freq = F,ylim=c(0,0.8),main="Central DAF: background vs anc candidates",xlab="DAF")  # first histogram
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,36), add=T,freq = F)  # second
legend("topright",legend = c("background","0.05% candidates"),fill = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)))

# nigeria chimp
p1 <- hist(round(sample.daf$n.daf*20),breaks = 20,probability = T,plot = F)
p2 <- hist(round(cands.daf$n.daf*20),breaks = 20,probability = T,plot = F)
plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,20),freq = F,ylim=c(0,0.8),main="Nigeria DAF: background vs anc candidates",xlab="DAF")  # first histogram
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,20), add=T,freq = F)  # second
legend("topright",legend = c("background","0.05% candidates"),fill = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)))


# western chimp
p1 <- hist(round(sample.daf$w.daf*22),breaks = 22,probability = T,plot = F)
p2 <- hist(round(cands.daf$w.daf*22),breaks = 22,probability = T,plot = F)
plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,22),freq = F,ylim=c(0,0.8),main="Western DAF: background vs anc candidates",xlab="DAF")  # first histogram
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,22), add=T,freq = F)  # second
legend("topright",legend = c("background","0.05% candidates"),fill = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)))
dev.off()