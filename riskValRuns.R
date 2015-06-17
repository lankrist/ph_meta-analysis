riskVal <- NULL
for (i in -3:0) {
  riskVal <- c(riskVal, 2**i)
}

for(risk in riskVal)
{
  system(paste("cat powerTests.r | R --slave --no-save --args 2000 200 0.2 ",risk," 1000 > TempOutput_for_Risk_Factor_",risk ,
               " &",sep=""))
 
}
