#meta-analysis simulation


#control_num = number of non-infected( = 0)
#r is risk for allele

simul_geno <- function(control, case, p, r){
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
  
  
  prob_healthy = (p^2)*(1-d)+2*p*q*((1-d)*r)+(q^2)*((1-d)*r^2)
  
  #prob_h0 is prob that genotype is 0 given that diseased
  prob_h0 = ((p^2)*(1-d))/prob_healthy
  #prob_h1 is prob that genotype is 1 given that diseased
  prob_h1 = (2*p*q*((1-d)*r))/prob_healthy
  #prob_h2 is prob that genotype is 2 given that diseased
  prob_h2 = ((q^2)*((1-d)*r^2))/prob_healthy
  
  genotype = c(prob_d0*case + prob_h0*control, prob_d1*case + prob_h1*control, prob_d2*case + prob_h2*control)  

  return(genotype)  
}

simul_geno(100,100, .5, .1)
