#generate parameters
#load data
dat <- read.csv("/Users/kristine/Documents/Summer_2015/genomics/meta_analysis/ph_meta-analysis/sim_gene.csv", header = T)
head(sims); table(sims)

#LOGISTIC REGRESSION
#returns log_reg beta_0, beta_1, log_reg error
log_reg <- function(v){
  lreg <- glm(formula = v[,2] ~ v[,1], family = binomial(logit), data = v)
  #coefficients
  lreg_beta <- lreg$coefficients
  #std_err
  lreg_err <- summary(lreg)$coefficients[,2][2]
  return(c(lreg_beta, lreg_err))
}

#SCORE TEST
#returns z-score, standard error, p val
score_test <- function(vector){
  score <- crossprod((vector[,1] - mean(vector[,1])), (vector[,2] - mean(vector[,2])))
  variance <- var(vector[,2])*crossprod((vector[,1] - mean(vector[,1])), (vector[,1] - mean(vector[,1])))
  p = 2*pnorm(-abs(score/variance))
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
  gentbl <- table(dat)
  #number of healthy individuals with respective genotype
  geno_h0 = gentbl[1,1]; geno_h1 = gentbl[2,1]; geno_h2 = gentbl[3,1]
  #number of infected individuals with respective genotype
  geno_d0 = gentbl[1,2]; geno_d1 = gentbl[2,2]; geno_d2 = gentbl[3,2]
  
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
#returns p value
chi_sq <- function(v){
  return(chisq.test(table(v))$p.value)
}
