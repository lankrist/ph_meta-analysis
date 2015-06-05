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
al_freq = 0.2
sims <- genegen(samp_size, al_freq)

#generate allele proportions and odds ratio
sim_prop <- table(sims)/samp_size
names(sim_prop)
sim_prop["1"]
odds <- sim_prop/sim_prop["0"]

alph = 2

#funcition that generates probability of a phenotype to occur in a certain marker
#takes in aray of marker samples, alph, beta
phenprob <- function(alpha, beta, marker){
  return(exp(alpha+beta*marker)/(1+exp(alpha+beta*marker)))
}

#generate conditional probability of illness given marker
marker_types = c(0, 1, 2)
logit_prob <- phenprob(alph, odds, marker_types)

#create phenotype probability
sim_phen <- c(
rbinom(length(which(sims ==0)), 1, logit_prob[1]),
rbinom(length(which(sims ==1)), 1, logit_prob[2]),
rbinom(length(which(sims ==2)), 1, logit_prob[3]))

sim_geno <- sort(sims)

table(sim_phen, sim_geno)/samp_size
geno <- data.frame(sim_geno, sim_phen)

#output data
write.table(geno, file = "/Users/kristine/Documents/Summer_2015/genomics/meta_analysis/ph_meta-analysis/sim_gene1.csv", sep=",")


