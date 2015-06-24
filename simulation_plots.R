#plotting

cases <- c(5000, 6000, 4000, 3000, 2000, 1000)
alfr <- c(0.01, 0.05, 0.10, 0.20, 0.50)

for (j in cases) {
  for (i in alfr) {
    plot(sims$Risk[which(intersect(sims$AF == i, sims$Ncases == j))], 
         sims$Power_Log[which(intersect(sims$AF == i, sims$Ncases == j))],
         xlab = "Risk", ylab = "Power")
    title(main = c("Logistic Rgression Power to risk AF = ", i), col.main = "blue")
    lines(sims$Power_Log[which(intersect(sims$AF == i, sims$Ncases == j))], 
          type="o", pch=22, lty=2, col="red")
  }  
}


