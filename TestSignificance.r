
#SIMULATION 


#takes in p-values and returns 1 if pval < 0.05
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

cosamp <- 1000
casamp <- 1000
al_fr <- 0.2
risk <- 1
rounds <- 200

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


sum(pValLog)
sum(pValScore)
sum(pValProp)
sum(pValChi)



