#meta-analysis simulation

#VARIABLES
control_samp = 1000
case_samp = 100


#FUNCTION
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



#SIMULATION
#simulate genotypes
sims <- simul_geno(100,100,.5,.1)

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

