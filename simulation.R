#meta-analysis 

#simulation of data
# 0: reference marker; 1: alternate marker
mark <- setClass("mark", representation(type = "numeric"))
ref <- new("mark", type = 0)
alt <- new("mark", type = 1)

#marker
marker <- setClass("marker", representation(i = "mark", j = "mark"))
m0 <- new("marker", i = ref, j = ref)
m1 <- new("marker", i = ref, j = alt)
m2 <- new("marker", i = alt, j = alt)

#set seed
set.seed(126)

#function returns percentage of case
genegen <- function(prob){
  return(mean(rpois(100, prob)))
}

genegen(5)

genome <- setClass("genome", representation(genes = "matrix"))
num <- 1000*genegen(15)
genes <- new("genome", genes = matrix(rep(0, num), rep(1, 1000-num)))


