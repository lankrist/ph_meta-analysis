#Goncala version riskvalRuns
source("powerTests2.r")

Casamp = 10000
Cosamp = 10000
rounds = 100

cat("Ncases Ncontrols AF Risk Power_Log Power_Score Power_Prop Power_Chi\n")

for (af in c(0.01, 0.05, 0.10, 0.20, 0.50))
  for (risk in c(1.0, 1.5, 2.0, 5.0))
    cat(Casamp, Cosamp, af, risk, powerReturn(Cosamp,Casamp,af,risk,rounds), "\n")
