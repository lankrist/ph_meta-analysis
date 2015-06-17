#PowerTesting

cosamp <- 1000
casamp <- 200
al_fr <- 0.1
rounds <- 200

risk <- NULL
for (i in -3:3) {
  risk <- c(risk, 2**i)
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


print(powerReturn(1000, 1000, 0.2, 1, 100))


pwLog=NULL
pwScore=NULL
pwChi=NULL
pwProp=NULL

for (j in risk)
{

   
  pwLog <- c(pwLog, sum(pValLog)/rounds)
#   pwScore <- c(pwScore, sum(pValScore)/rounds)
#   pwProp <- c(pwProp, sum(pValProp)/rounds)
#   pwChi <- c(pwChi, sum(pValChi)/rounds)
  
}

