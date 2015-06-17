riskVal <- NULL
for (i in -3:3) {
  riskVal <- c(risk, 2**i)
}

for(risk in riskVal)
{
  
  system(paste("cat NewPowerGenerator.r | R --slave --no-save --args 2000 2000 0.2 ",risk," 10000 > TempOutput_for_Risk_Factor_",risk ,
               " &",sep=""))
  
}