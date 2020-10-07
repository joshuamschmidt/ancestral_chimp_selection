library(rtracklayer)
library(data.table)
library(R.utils)
#
setwd("~/Projects/ancestral_chimp_selection/ens_reg")
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
info <- unique(fread("ens.reg.elements.txt.gz",sep="\t",header = T, fill = TRUE,select = c("Feature type","Regulatory stable ID","SO term accession")))
names(info) <- c("type","id","SO")
setkey(info,id)
panTro_elements <- panTro_elements[info][!chr %in% NA]

# create a basic go lookup file for GOWINDA
reg_types <- panTro_elements[,unique(type)]
if (file.exists("all_ens.reg.elements.SO.txt")){
  file.remove("all_ens.reg.elements.SO.txt")
}
for (rt in reg_types){
  so_term <- panTro_elements[type==eval(rt)][1]$SO
  elements <- paste(sort(panTro_elements[type==eval(rt)][,unique(id)]),collapse = ' ')
  fwrite(x = data.table(so_term,rt,elements),
         file = "all_ens.reg.elements.SO.txt",
         append = T,
         quote = F,
         row.names = F,
         col.names = F,
         sep = "\t")
}

# create a gtf file for go
panTro_elements[, chr.s := tstrsplit(chr, "r",keep = 2,fixed=T)]
if (file.exists("all_ens.SO.definitions.gtf")){
  file.remove("all_ens.SO.definitions.gtf")
}

fwrite(x = panTro_elements[,.(chr=chr.s,
                              source='josh',
                              type='exon',
                              start,
                              end,
                              filter=".",
                              strand="+",
                              phase="0",
                              desc=paste('gene_id \"',id,'\";',sep=""))],
       file = "all_ens.SO.definitions.gtf",
       quote = F,
       col.names = F,
       row.names = F,
       sep = "\t")



