#generate parameters

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


