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

rep(ref,10)

genome <- setClass("genome", representation(genes = "matrix"), prototype = list(genes = (matrix(rep(ref, 1000)))))
num <- 1000*genegen(15)
genes <- new("genome", genes = matrix(rep(0, num), rep(1, 1000-num)))



#function returns genome
genegen <- function(af){
  hap1=rbinom(1000,1,af)
  hap2=rbinom(1000,1,af)
  gene = hap1+hap2
  return(gene)
}

genegen(0.4)




