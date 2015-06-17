
rgs <- commandArgs(TRUE)
if (length(Args)!=5){
  print ("Usage: CaseNo ControlNo Allele_Freq Risk_Ratio Rounds")
  quit()
}

Cosamp=as.numeric(Args[1])
Casamp=as.numeric(Args[2])
Al_fr <- as.numeric(Args[3])
Risk<-as.numeric( Args[4])
Rounds <-as.numeric(Args[5])


#cosamp, casamp, al_fr, risk, roundscosamp

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
  disease <- sample(c(0,1,2), case, prob=c(prob_d0, prob_d1, prob_d2),replace=TRUE)
  
  prob_healthy = (p^2)*(1-d)+2*p*q*(1-d*r)+(q^2)*(1-d*r^2)
  
  #prob_h0 is prob that genotype is 0 given that diseased
  prob_h0 = ((p^2)*(1-d))/prob_healthy
  #prob_h1 is prob that genotype is 1 given that diseased
  prob_h1 = (2*p*q*(1-d*r))/prob_healthy
  #prob_h2 is prob that genotype is 2 given that diseased
  prob_h2 = ((q^2)*(1-d*r^2))/prob_healthy
  
  #proportions of healthy
  healthy <- sample(c(0,1,2), control, prob=c(prob_h0, prob_h1, prob_h2),replace=TRUE)
  
  #geno = c(case*prob_d0 + control*prob_h0, case*prob_d1 + control*prob_h1, case*prob_d2+control*prob_h2)
  #return(geno)
  genotype = c( disease, healthy)
  return(genotype)
  
}
#put all the simulations together
simulation <- function(control_samp, case_samp, dom_allele, risk_factor){
  gprob <- geno_prob(control_samp,case_samp, dom_allele, risk_factor)
  cosamp=control_samp
  casamp=case_samp
  pprob <- c(rep(1, cosamp), rep(0, casamp))
  return(data.frame(sim_geno = gprob, sim_pheno = pprob))
}



#generate parameters
#LOGISTIC REGRESSION                    
#returns log_reg beta_0, beta_1, log_reg error
log_reg <- function(v){
  lreg <- glm(formula = v[,2] ~ v[,1], family = binomial(logit), data = v)
  #coefficients
  lreg_beta <- summary(lreg)$coefficients[,1][2]
  #std_err
  lreg_err <- summary(lreg)$coefficients[,2][2]
  pVal <- summary(lreg)$coefficients[,4][2]
  return(c(lreg_beta, lreg_err,pVal))
}

#SCORE TEST
#returns z-score, standard error, p val
score_test <- function(vector){
  score <- crossprod( (vector[,1] - mean(vector[,1])), (vector[,2] - mean(vector[,2])) )
  StdErr <- sqrt( crossprod( (vector[,2] - mean(vector[,2]))^2, (vector[,1] - mean(vector[,1]))^2 ) )
  p = 2*pnorm(-abs(score/StdErr))
  return(c(score/StdErr, sqrt(StdErr), p))
}

#ALLELE FREQUENCY TEST
#risk allele is q
#input vector that is 2 x n where n is the number of samples(controls+cases)
#returns standard error, p-value, z score
allele_freq <- function(dat){
  
  Case=which(dat[,2]==1)
  Ctrl=which(dat[,2]==0)
  NCase=length(Case)
  NCtrl=length(Ctrl)
  
  PCase=sum(dat[Case,1])/(2*NCase)
  PCtrl=sum(dat[Ctrl,1])/(2*NCtrl)
  
  stand_error = sqrt(  (PCase*(1-PCase))/(2*NCase) + (PCtrl*(1-PCtrl))/(2*NCtrl) )
  
  z = (PCase - PCtrl)/stand_error
  p = 2*pnorm(-abs(z))
  result = c(stand_error, z,p)
  return(result)
}

#CHI SQUARE TEST
#returns p value
chi_sq <- function(v){
  return(chisq.test(table(v))$p.value)
}

PValInd = function(x)
{
  if(x<0.05)
  {
    return(1)
  }
  else
  {
    return(0)
  }
}

powerReturn <- function(cosamp, casamp, al_fr, risk, rounds){
  
  pValLog=NULL
  pValScore=NULL
  pValChi=NULL
  pValProp=NULL
  
  
  for (i in 1:rounds) {
    data=simulation(cosamp, casamp, al_fr, risk)
    
    temp = log_reg(data)
    pValLog=c(pValLog,PValInd(temp[3]))
    
    temp = score_test(data)
    pValScore=c(pValScore,PValInd(temp[3]))
    
    temp = allele_freq(data)
    pValProp=c(pValProp,PValInd(temp[3]))
    
    temp = chi_sq(data)
    pValChi=c(pValChi,PValInd(temp[1]))
  }
  
  return(c(mean(pValLog), mean(pValScore), mean(pValProp), mean(pValChi)))
}


print(powerReturn(Cosamp,Casamp,Al_fr,Risk,Rounds))
                                      