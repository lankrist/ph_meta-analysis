#plotting

alfr <- 0.2
plot(sims$Risk[which(sims$AF == alfr)], sims$Power_Log[which(sims$AF == alfr)],
     xlab = "Risk", ylab = "Power")
title(main = c("Logistic Rgression Power to risk AF = ", alfr), col.main = "blue")
lines(sims$Risk[which(sims$AF == alfr)], type="o", pch=22, lty=2, col="red")
