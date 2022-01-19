rm(list=ls()) #reset working space
graphics.off() #closing current or all graphical windows


subst_data = read.table("marked_duplicates.error_rate_by_read_position.txt", header = TRUE, sep = "", dec = ".")



ones = subst_data[subst_data$read_number == 1 ,]

twos = subst_data[subst_data$read_number == 2 ,]


par(mfrow=c(1,1))

par(new = FALSE)


colours = seq(1,12)


jpeg("read_1.jpg",width = 800, height = 600)

for(i in ones[6:17]){
  plot(ones$position,i,type = "l", lty = 1, main = "Read 1", xlab = "Position", ylab = "Error Rate", axes = FALSE, col=sample(colours,1, replace = FALSE))
  par(new = TRUE)
  plot(ones$position,ones$error_rate,type = "l", lty = 1, lwd=5, col="red", xlab = "Position", ylab = "Error Rate")
  par(new = TRUE)
}
legend("topleft", legend="Error Rate", col="red", lty = 1, lwd = 5)

dev.off()

par(new = FALSE)

jpeg("read_2.jpg",width = 800, height = 600)

for(i in twos[6:17]){
  plot(twos$position,i,type = "l", lty = 1, main = "Read 2", xlab = "Position", ylab = "Error Rate", axes = FALSE, col=sample(colours,1, replace = FALSE))
  par(new = TRUE)
  plot(twos$position,twos$error_rate,type = "l", lty = 1, lwd=5, col="red", xlab = "Position", ylab = "Error Rate")
  par(new = TRUE)
}
legend("topleft", legend="Error Rate",
       col="red", lty = 1, lwd = 5)

dev.off()



jpeg("reads.jpg",width = 800, height = 600)


plot(twos$position,twos$error_rate,type = "l", lty = 1, lwd=3, col=1, xlab = "Position", ylab = "Error Rate", main = "Error Rates")
par(new = TRUE)
plot(ones$position,ones$error_rate,type = "l", lty = 3, lwd=3, col="red", xlab = "Position", ylab = "Error Rate")



legend("topleft", legend=c("Read 1", "Read 2"), col=c("red","black"), lty = c(3,1), lwd = 5)


dev.off()
