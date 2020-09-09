library(rtracklayer)
library(data.table)
#
setwd("/Users/joshuaschmidt/Projects/ancestral_chimpanzee_selection/ens_reg")
hg19_elements <- fread("unique_elements.chr.bed",col.names = c("chr","start","end","id") )
panTro_elements <- fread("pantro4.unique_elements.chr.bed",col.names = c("chr","start","end","id") )
remaphg19_elements <- fread("pantro4.unique_elements_remaphg19.chr.bed",col.names = c("chr","start","end","id") )
setkey(hg19_elements,id)
setkey(panTro_elements,id)
setkey(remaphg19_elements,id)

# filter for map to homologous chr in panTro, same chr and equal to or contained within the original location
hc_map_ids <- hg19_elements[panTro_elements][i.chr==chr][,.(chr,start,end,id)][remaphg19_elements][i.chr==chr][i.start>=start & i.end <= end]$id

panTro_elements <- panTro_elements[id %in% hc_map_ids]
setkey(panTro_elements,id)
# tie these to feature type
info <- unique(fread("ens.reg.elements.txt.gz",sep="\t",header = T, fill = TRUE,select = c("Feature type","Regulatory stable ID")))
names(info) <- c("type","id")
setkey(info,id)
panTro_elements <- panTro_elements[info][!chr %in% NA]