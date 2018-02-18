barplot(t(abs(wSparse)), col="blue", ylim=c(-1,1))
par(new=TRUE)
barplot(-t(abs(wRidge)), col="orange", ylim=c(-1,1))
legend("topright", c('Lasso', 'Ridge'), col=c("blue", "orange"), pch=15)