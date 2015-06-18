#Goncala version riskvalRuns

source("powerTests2.r")

Casamp = 2000
Cosamp = 500
rounds = 100
disfr = 0.0001

cat("Ncases Ncontrols AF Risk Power_Log Power_Score Power_Prop Power_Chi Disease_frequency\n")

for (af in c(0.01, 0.05, 0.10, 0.20, 0.50)){
  for (risk in c(0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0, 1.125, 1.25, 1.375, 1.5, 1.75, 2.0, 5.0)){
    cat(Casamp, Cosamp, af, risk, powerReturn(Cosamp,Casamp,af,risk,rounds, disfr), "\n")
  }
  cat("\n")
}