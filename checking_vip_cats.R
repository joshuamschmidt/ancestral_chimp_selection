
library(data.table)
library(ggplot2)
library(hrbrthemes)
library(scales)
setwd("/Users/joshuaschmidt/Projects/bat_exomes/xmas_island_bams/cov_stats/")

files <- list.files(pattern="dedup.rg.fixed.bamcov.txt")
all.cov <- data.table()
for (f in files){
  tmp <- fread(f,header=TRUE,select = c("#rname","meandepth"))
  names(tmp) <- c("target","meandepth")
  sample <- strsplit(f,split = "_")[[1]][2]
  tmp[,sample:=sample]
  all.cov <- rbindlist(list(all.cov,tmp))
}

# y-axis is sample, x axis is target. sort both in ascending order
y_order <- all.cov[,mean(meandepth),by=sample][order(V1)]$sample
x_order <- all.cov[,mean(meandepth),by=target][order(V1)]$target
all.cov$sample <- factor(all.cov$sample, levels = y_order)
all.cov$target <- factor(all.cov$target, levels = x_order)
p1 <- ggplot(all.cov, aes(target, sample, fill=meandepth)) + geom_tile() +
  scale_fill_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0, trans = 'log10') + theme_ipsum() + theme(axis.title.x=element_blank(),
                                                                                                                      axis.text.x=element_blank(),
                                                                                                                      axis.ticks.x=element_blank())

pdf(file="coverage_target.by.sample.pdf",height=10,width=12)
p1
dev.off()

png(file="coverage_target.by.sample.png",width = 12, height = 10, units = "in", pointsize = 12,res=600)
p1
dev.off()