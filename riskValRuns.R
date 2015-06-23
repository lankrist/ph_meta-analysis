#Goncala version riskvalRuns
source("powerTests.r")
#source("powerTests2.r")

#total sample size 10000
Samp <- 10000  #Casamp = n; Cosamp = Samp-n
rounds = 10000
disfr = 0.0001

sims <- data.frame()
#cat("Ncases Ncontrols AF Risk Power_Log Power_Score Power_Prop Power_Chi Disease_frequency\n")
for (Cosamp in c(9000, 8000, 7000, 6000, 5000)){
  Casamp <- Samp - Cosamp
  for (af in c(0.01, 0.05, 0.10, 0.20, 0.50)){
    for (risk in c(0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0, 1.125, 1.25, 1.375, 1.5, 1.75, 2.0, 5.0)){
      #cat(Casamp, Cosamp, af, risk, powerReturn(Cosamp,Casamp,af,risk,rounds, disfr), "\n")
      sims <- rbind(sims,c(Casamp, Cosamp, af, risk, powerReturn(Cosamp,Casamp,af,risk,rounds, disfr), disfr))
    }
  }
}
names(sims) <- c("Ncases","Ncontrols", "AF", "Risk", "Power_Log", "Power_Score", "Power_Prop", "Power_Chi", "Disease_frequency")

write.table(sims, file = "sim_fail.csv", sep = ",", col.names = T, row.names = F)