#meta-analysis 

#function returns a simulated marker
#takes in number of samples and allele frequency
genegen <- function(samples, af){
  hap1=rbinom(samples,1,af)
  hap2=rbinom(samples,1,af)
  gene = hap1+hap2
  return(gene)
}

samp_size = 100000
prop = 0.2
sims <- genegen(samp_size, prop)

#generate allele proportions and odds ratio
sim_prop <- table(sims)/samp_size
names(sim_prop)
sim_prop["1"]
odds <- sim_prop/sim_prop["0"]

alph = 0.00005

#funcition that generates probability of a phenotype to occur in a certain marker
#takes in aray of marker samples, alph, beta
phenprob <- function(alpha, beta, marker){
  return(exp(alpha+beta*marker)/(1+exp(alpha+beta*marker)))
}

marker_types = c(0, 1, 2)

phenprob(alph, odds, marker_types)




