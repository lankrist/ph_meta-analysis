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
  score <- crossprod((vector[,1] - mean(vector[,1])), (vector[,2] - mean(vector[,2])))
  StdErr <- sqrt(var(vector[,2])*crossprod((vector[,1] - mean(vector[,1])), 
                                           (vector[,1] - mean(vector[,1]))))
  p = 2*pnorm(-abs(score/StdErr))
  return(c(score/variance, sqrt(variance), p))
}

#ALLELE FREQUENCY TEST
#function that takes in matrix, genotype, phenotype and returns the number of samples
samp_count <- function(v, geno, pheno){
  length(intersect(which(v[,1] == geno)), which(v[ , 2] == pheno))
}
#actual function 
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
  
  stand_error = sqrt(  (PCase*(1-PCase) + PCtrl*(1-PCtrl)) / (2*(NCase+NCtrl)) )
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




