




k <- 200
t <- seq(min(min(f),min(f0))-0.1, 
         max(max(f),max(f0))+0.1,
         length=k)

par(mfrow=c(2,1))

hist_input <- hist(f0, t, plot=FALSE)
hist_input$counts <- (k/n**2)*hist_input$counts
plot(hist_input, col="blue", main="Input", xlim=c(0,0.5))

hist_synthesized <- hist(f, t, plot=FALSE)
hist_synthesized$counts <- (k/n**2)*hist_synthesized$counts
plot(hist_synthesized, col="blue", main="Input", xlim=c(0,0.5))
