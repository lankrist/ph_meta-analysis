#generate parameters
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

#ALLELE FREQUENCY COMPARISON


#CHI SQUARE TEST



