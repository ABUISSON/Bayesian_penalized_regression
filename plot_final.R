#We charge the results from the different models

setwd("/Users/arnaudb/Desktop/montecarlo")
sampling = read.xls ("IS.xls", sheet=1, header=FALSE)
resampling = read.xls ("resampling.xls", sheet=1, header=FALSE)
gibbs = read.xls ("Gibbssampler.xls", sheet=1, header=FALSE)
MH = c(-0.09761304,0.10753402 ,0.02588428,-0.21344828,0.29674087,0.00788370,-0.32460173,0.28767768,-0.22347610, -0.23548339, 0.09453886,-0.41330954)


#we compare the values of the beta
plot(1:12,gibbs[,3], pch=10, col='#8FBC8F', cex=0.5, ylim=c(-0.5,0.4), xlab="", ylab="")
points(c(1:12)+0.1, resampling[,1],pch=10, col='red' , cex=0.5 )
points(c(1:12)-0.1, MH, pch=10, col='#008000',cex=0.5)
points(c(1:12)+0.2,sampling[,1], pch=10, col='blue', cex=0.5)
for (k in 1:12){
  lines(c(k,k), c(gibbs[k,1],gibbs[k,2]) , "l", col='black')
}
