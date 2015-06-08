#meta-analysis simulation

#VARIABLES
samp_size = 100000
al_freq = 0.2
alpha = 2
marker_types = c(0, 1, 2)

#FUNCTION
#Generates markers for a given sample size
#takes in number of samples and allele frequency
genegen <- function(samples, af){
  hap1=rbinom(samples,1,af)
  hap2=rbinom(samples,1,af)
  gene = hap1+hap2
  return(gene)
}

#already know the phenotype
#control_num = number of non-infected( = 0)
#r is risk for allele
simul_geno <- function(control, case, p, r){
  #d is prevalance of disease
  d = case/(control+case)
  q = 1-p
  
  prob_disease = d*(p^2) + 2*p*1*d*r + d*(r^2)*(q^2)
  
  #prob_0 is prob that genotype is 0 given that diseased
  prob_0 = (d*p^2)/prob_disease
  #prob_1 is prob that genotype is 1 given that diseased
  prob_1 = (2*p*q*d*r)/prob_disease
  #prob_2 is prob that genotype is 2 given that diseased
  prob_2 = ((q^2)*d*r^2)/prob_disease
  
  return(prob_0)
}

simul_geno(100,100,.5,.1)




#Generates probability of a phenotype to occur in a certain marker
#takes in array of marker samples, alph, beta
phenprob <- function(alpha, beta, marker){
  return(exp(alpha+beta*marker)/(1+exp(alpha+beta*marker)))
}

#SIMULATION
#simulate genotypes
sims <- genegen(samp_size, al_freq)

#calculate allele proportions and odds ratio
sim_prop <- table(sims)/samp_size
odds <- sim_prop/sim_prop["0"]

#calculate conditional probability of illness given marker
logit_prob <- phenprob(alpha, odds, marker_types)

#simulate phenotypes
sim_phen <- c(
  rbinom(length(which(sims ==0)), 1, logit_prob[1]),
  rbinom(length(which(sims ==1)), 1, logit_prob[2]),
  rbinom(length(which(sims ==2)), 1, logit_prob[3]))

#sort genotype data simulation
sim_geno <- sort(sims)

#compare data and convert to dataframe
table(sim_phen, sim_geno)/samp_size
geno <- data.frame(sim_geno, sim_phen)

#output data
#write.table(geno, file = "/Users/kristine/Documents/Summer_2015/genomics/meta_analysis/ph_meta-analysis/sim_gene4.csv", row.names = F, sep=",")

