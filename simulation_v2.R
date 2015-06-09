#meta-analysis simulation 2
library(mosaic)

#FUNCTION
#Returns genotype proportions
#control_num = number of non-infected( = 0)
#r is risk for allele; p is ...
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
  disease <- do(control)*sample(c(0,1,2), 1, prob=c(prob_d0, prob_d1, prob_d2))
  
  prob_healthy = (p^2)*(1-d)+2*p*q*((1-d)*r)+(q^2)*((1-d)*r^2)
  
  #prob_h0 is prob that genotype is 0 given that diseased
  prob_h0 = ((p^2)*(1-d))/prob_healthy
  #prob_h1 is prob that genotype is 1 given that diseased
  prob_h1 = (2*p*q*((1-d)*r))/prob_healthy
  #prob_h2 is prob that genotype is 2 given that diseased
  prob_h2 = ((q^2)*((1-d)*r^2))/prob_healthy
  
  #proportions of healthy
  healthy <- do(control)*sample(c(0,1,2), 1, prob=c(prob_h0, prob_h1, prob_h2))
  
  genotype = mapply(c, disease, healthy)  

  return(genotype)  
}

#SIMULATION
gprob <- geno_prob(100,100, .5, .1)
pprob <- c(rep(1, 100), rep(0, 100))
geno <- data.frame(gprob, pprob)
write.table(geno, 
            file = "/Users/kristine/Documents/Summer_2015/genomics/meta_analysis/ph_meta-analysis/sim_gene4.csv", 
            row.names = F, sep=",")

