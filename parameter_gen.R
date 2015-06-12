#generate parameters
#load data
sims <- read.csv("/Users/kristine/Documents/Summer_2015/genomics/meta_analysis/ph_meta-analysis/sim_gene.csv", header = T)
head(sims); table(sims)

#LOGISTIC REGRESSION
log_reg <- function(v){
  lreg <- glm(formula = v[,2] ~ v[,1], family = binomial(logit), data = v)
  #coefficients
  lreg_beta <- lreg$coefficients[2]
  #std_err
  lreg_err <- summary(lreg)$coefficients[,2][2]
  return(c(lreg_beta, lreg_err))
}
log_reg(sims)

#SCORE TEST
score_test <- function(vector){
  score <- crossprod((vector[,1] - mean(vector[,1])), (vector[,2] - mean(vector[,2])))
  variance <- crossprod((vector[,1] - mean(vector[,1])), (vector[,1] - mean(vector[,1])))
  p = 2*pnorm(-abs(score/variance))
  return(c(score/variance, sqrt(variance), p))
}


#ALLELE FREQUENCY TEST
#risk allele is q
#input vector that is 2 x n where n is the number of samples(controls+cases)
allele_freq <- function(dat){
  
  #number of healthy individuals with respective genotype
  geno_h0 = length((intersect(which(dat[ , 2] == 0), which(dat[,1] ==0))))
  geno_h1 = length((intersect(which(dat[ , 2] == 0), which(dat[,1] ==1))))  
  geno_h2 = length((intersect(which(dat[ , 2] == 0), which(dat[,1] ==2))))
  
  #number of infected individuals with respective genotype
  geno_d0 = length((intersect(which(dat[ , 2] == 1), which(dat[,1] ==0))))
  geno_d1 = length((intersect(which(dat[ , 2] == 1), which(dat[,1] ==1))))
  geno_d2 = length((intersect(which(dat[ , 2] == 1), which(dat[,1] ==2))))
  
  #checked the geno_h0, geno_h1, etc and it worked
  
  tot_healthy = geno_h0+geno_h1+geno_h2
  tot_infect = geno_d0+geno_d1+geno_d2
  
  risk_control = (2*geno_h2+geno_h1)/(2*tot_healthy)
  
  risk_cases = (2*geno_d2+geno_d1)/(2*tot_infect)
  
  stand_error = sqrt((risk_cases*(1-risk_cases) + risk_control*(1-risk_control)) / (2*(tot_healthy+tot_infect)) )
  
  z = (risk_cases - risk_control)/stand_error
  
  p = 2*pnorm(-abs(z))
  
  result = c(stand_error, p, z)
  
  return(result)
  
}


#CHI SQUARE TEST
chi_sq <- function(v){
  return(chisq.test(table(v))$p.value)
}


