#generate parameters
<<<<<<< HEAD

#load data
sims <- read.csv("/Users/kristine/Documents/Summer_2015/genomics/meta_analysis/ph_meta-analysis/sim_gene1.csv", header = T)
head(sims); table(sims)

#LOGISTIC REGRESSION
lreg <- glm(formula = sim_phen ~ sim_geno, family = binomial(logit), data = sims)
names(lreg)
summary(lreg)
#coefficients
lreg_beta <- lreg$coefficients
#std_err
lreg_err <- mean(lreg$residuals)
#deviance
lreg_deviance <- lreg$deviance
=======
install.packages("statmod")
Rlibrary(statmod)

#load data
sims <- read.csv("/Users/kristine/Documents/Summer_2015/genomics/meta_analysis/ph_meta-analysis/sim_gene.csv", header = T)
head(sims); table(sims)

#LOGISTIC REGRESSION
lreg <- glm(formula = sim_pheno ~ result, family = binomial(logit), data = sims)
names(lreg);summary(lreg)

#coefficients
lreg_beta <- lreg$coefficients
#std_err
lreg_err <- summary(lreg)$coefficients[,2]

#SCORE TEST
score_test <- function(vector){
  gen_mean <- mean(vector$result)
  phen_mean <- mean(vector$sim_pheno)
  return(crossprod((vector[,1] - mean(vector[,1])), (vector[,2] - mean(vector[,2]))))
}

score_test(sims)
#sims <- data.frame(result = c(11, 15, 3, 7), sim_pheno = c(10, 12, 19, 4))
#crossprod((sims[,1] - mean(sims[,1])), (sims[,2] - mean(sims[,2])))

#ALLELE FREQUENCY TEST

#risk allele is q


#create a 2 x 10 vector

samp = data.frame(geno = c(0,1,2,0,1,1,0,2,1,0), pheno = c(0,0,1,0,0,0,0,1,0,0))

geno_h0 = length((intersect(which(samp[ , 2] == 0), which(samp[,1] ==0))))
geno_h1 = length((intersect(which(samp[ , 2] == 0), which(samp[,1] ==1))))

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
  
  result = c(stand_error, z)
  
  return(result)
  
}

allele_freq(samp)


#Chi-Squared Test
v <- NULL
gprob <- geno_prob(cosamp,casamp, .5, .1)
pprob <- c(rep(1, cosamp), rep(0, casamp))
geno <- data.frame(sim_geno = gprob, sim_pheno = pprob)
tbl <- table(geno)
Xsq <- chisq.test(tbl)
pval <- Xsq$p.value
v <- c(v, pval)


#CHI SQUARE TEST

>>>>>>> a649cb5496c86d825ec2e3d115e92d38cf1211ac


