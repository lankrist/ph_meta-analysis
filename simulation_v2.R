#meta-analysis simulation
library(mosaic)

#VARIABLES
cosamp <- 600
casamp <- 400

#FUNCTION
#Returns genotype proportions
#control_num = number of non-infected( = 0)
#r is risk for allele; p is allele frequency 
geno_prob <- function(control, case, p, r){
  #d is prevalance of disease
  d = case/(control+case)
  q = 1-p
  
  prob_disease = d*(p^2) + 2*p*q*d*r + d*(r^2)*(q^2)
  #prob_d0 is prob that genotype is 0 given that diseased
  prob_d0 = (d*p^2)/prob_disease
  #prob_d1 is prob that genotype is 1 given that diseased
  prob_d1 = (2*p*q*d*r)/prob_disease
  #prob_d2 is prob that genotype is 2 given that diseased
  prob_d2 = ((q^2)*d*r^2)/prob_disease
  #proportions of disease
  disease <- do(case)*sample(c(0,1,2), 1, prob=c(prob_d0, prob_d1, prob_d2))
  
  prob_healthy = (p^2)*(1-d)+2*p*q*(1-d*r)+(q^2)*(1-d*r^2)

  #prob_h0 is prob that genotype is 0 given that diseased
  prob_h0 = ((p^2)*(1-d))/prob_healthy
  #prob_h1 is prob that genotype is 1 given that diseased
  prob_h1 = (2*p*q*(1-d*r))/prob_healthy
  #prob_h2 is prob that genotype is 2 given that diseased
  prob_h2 = ((q^2)*(1-d*r^2))/prob_healthy

  #proportions of healthy
  healthy <- do(control)*sample(c(0,1,2), 1, prob=c(prob_h0, prob_h1, prob_h2))
  
  #geno = c(case*prob_d0 + control*prob_h0, case*prob_d1 + control*prob_h1, case*prob_d2+control*prob_h2)
  #return(geno)
  genotype = mapply(c, disease, healthy)
  return(genotype)  

}

#put all the simulations together
simulation <- function(control_samp, case_samp, dom_allele, risk_factor){
  gprob <- geno_prob(control_samp,case_samp, dom_allele, risk_factor)
  pprob <- c(rep(1, cosamp), rep(0, casamp))
  return(data.frame(sim_geno = gprob, sim_pheno = pprob))
}

#SIMULATION
#this portion should be written as function too
geno <- simulation(cosamp, casamp, 0.4, 0.1)
table(geno)

write.table(geno, 
            file = "/Users/kristine/Documents/Summer_2015/genomics/meta_analysis/ph_meta-analysis/sim_gene.csv", 
            row.names = F, sep=",")

