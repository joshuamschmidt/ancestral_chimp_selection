library(data.table)
th_candidates <- fread("/Users/joshuaschmidt/Projects/ancestral_chimp_selection/th_set_test/th_candidate_classification.txt")

names(th_candidates) <- c("gene","th_class","sweep_candidate")
# ensure overlap with my KEGG testing
th_candidates <- th_candidates[gene %in% KEGG[id=="hsa04658"][,strsplit(genes," ")]$V1]
# exlude genes in both pathways
th_candidates <- th_candidates[th_class!=3]
th_candidates <- th_candidates[order(-sweep_candidate)]
# ordered so non-candidates are row 1-n, candisates n+1 -> m
n.genes <- dim(th_candidates)[1]
n.candidates <- th_candidates[sweep_candidate==1,.N]
th_candidates[1:n.candidates]

th_candidates[,.N,by=.(sweep_candidate,th_class)]
# proportion of candiates which are th1 specific.
# simplfy to th1secific vs non
#th_candidates[th_class==3,th_class:=2]
obs <- th_candidates[1:n.candidates][th_class==1,.N]

# samples
n.perm <- 1e3
perms <- rep(0.00, n.perm)
for(i in 1:n.perm){
  perms[i] <- th_candidates[][sample(.N,n.genes)][1:n.candidates][th_class==1,.N]
}

hist(perms)

# test that Th2 is underrepresented in candidates (NULL candidate gene equally likely to be Th1 specific vs other)
pvalue <- (1+length(perms[perms >= obs]))/(n.perm+1)
# pvalue <- (1+length(perms[perms <= obs]))/(n.perm+1)
# KEGG[id=="hsa04658"]


# Chi-square test

M <- as.table(rbind(c(7, 18), c(8, 47)))
dimnames(M) <- list(candidate = c("Y", "N"),
                    th = c("1","Other"))
