setwd("~/QtProjects/build-hello-Desktop_Qt_5_3_GCC_64bit-Debug")
letter="prueba"
folder="svg"

perfil <- read.table(paste("profile_",letter,".dat",sep=""))
N=as.numeric(perfil[1,1])
PROM=as.numeric(perfil[1,2])
SUPERTIME=as.numeric(perfil[1,3])
BLOCKLENGTH=as.numeric(perfil[1,4])
SAMPLE=as.numeric(perfil[1,5])
MUTTABLESIZE=as.numeric(perfil[1,6])
B=as.numeric(perfil[1,7])
BIGTIME=as.numeric(perfil[1,8])
BURNIN=as.numeric(perfil[1,9])
STRUCTURED=as.numeric(perfil[1,10])
mu=as.numeric(perfil[1,11])
rho=as.numeric(perfil[1,12])
kappa=as.numeric(perfil[1,13])
lambda=as.numeric(perfil[1,14])

R=4*N*rho*BLOCKLENGTH
C=4*N*kappa*lambda
 
numoflines=SUPERTIME
runtime=0
 
tiempo0<-0
tiempo1<-BURNIN
tiempo2<-BIGTIME-BIGTIME/2
tiempo3<-BIGTIME

dty <- read.table(paste("samplepi0[",SAMPLE,"]_",letter,".dat",sep=""))
pi0 <- colMeans(dty)
dty <- read.table(paste("samplepi1[",SAMPLE,"]_",letter,".dat",sep=""))
pi1 <-  colMeans(dty)
dty <- read.table(paste("samplepi2[",SAMPLE,"]_",letter,".dat",sep=""))
pi2 <-  colMeans(dty)
dty <- read.table(paste("samplepi3[",SAMPLE,"]_",letter,".dat",sep=""))
pi02 <-  colMeans(dty)
#dty <- read.table(paste("dupFreq_",letter,".dat",sep=""))
#dupli <- mean(dty[1:numoflines,])

D0<-pi0
D1<-pi1
D2<-pi2#
D2
D02<-pi02
#dupli <- dupli/(2*N)

timelength <- dim(dty)[2]
timelength
time <- seq(1,timelength,1)
maximus<-max(D1)
theta=4*N*mu*BLOCKLENGTH
maximus=2*theta
#maximus<-mean(D02[tiempo2:tiempo3])

svg(paste("piSample_",folder,"_",letter,".svg"))

plot(time,D0,col="darkgreen",xlab="time", ylab="Average pairwise differences", type="l",xlim=c(-1,timelength+.1*timelength),ylim=c(-.005,maximus+.1*maximus),cex=0.5,pch=20)
points(time,D1,col="red",type="l",cex=0.5,pch=20)
points(time,D2,col="purple",type="l",cex=0.5,pch=20)
points(time,D02,col="blue",type="l",cex=0.5,pch=20)
#points(time,dupli,col="brown",type="l",cex=0.5,pch=20)

theta=4*N*mu*BLOCKLENGTH
abline(h=2*theta, col ="brown")

abline(h=theta, col="brown")
abline(v=tiempo2, col="brown")
#abline(v=30, col="black")
#abline(v=31, col="black")
#abline(v=32, col="black")

points(timelength+(timelength*.02),mean(D1[tiempo2:tiempo3]),col="red",cex=1)
points(timelength+(timelength*.02),mean(D0[tiempo2:tiempo3]),col="darkgreen",cex=1)
points(timelength+(timelength*.02),mean(D2[tiempo2:tiempo3]),col="purple",cex=1)
points(timelength+(timelength*.02),mean(D02[tiempo2:tiempo3]),col="blue",cex=1)
promedio02=(mean(D0[tiempo2:tiempo3])+mean(D2[tiempo2:tiempo3]))/2

distanciaNum=.12
text(timelength+(timelength*distanciaNum),mean(D1[tiempo2:tiempo3]),mean(D1[tiempo2:tiempo3]),col="red",cex=0.6)
text(timelength+(timelength*distanciaNum),mean(D0[tiempo2:tiempo3]),mean(D0[tiempo2:tiempo3]),col="darkgreen",cex=0.5)
text(timelength+(timelength*distanciaNum),mean(D2[tiempo2:tiempo3]),mean(D2[tiempo2:tiempo3]),col="purple",cex=0.5)
text(timelength+(timelength*distanciaNum),mean(D02[tiempo2:tiempo3]),mean(D02[tiempo2:tiempo3]),col="blue",cex=0.5)
text(timelength,mean(D0[tiempo2:tiempo3])-.2*mean(D0[tiempo2:tiempo3]),promedio02,col="black",cex=1)

legend("topleft", c("Single-copy", "Original", "Duplicated","O+D"), col=c("red","darkgreen","purple","blue"), pch=20)
title(main=bquote(paste(.(folder), "_",.(letter)," N=", .(N), " Prom=", .(PROM), " Runs=", .(numoflines))), #, " Runtime=", .(runtime),"hrs")),
sub=bquote(paste(theta, "= ", .(theta), " ", R, "= ", .(R), " ", C, "= ", .(C), " ", lambda, "= ", .(lambda), " ", BlockL, "= ", .(BLOCKLENGTH)))) 

#sub=bquote(paste(mu, "= ", .(mu)*4*.(N), " ", rho, "= ", .(rho), " ", kappa, "= ", .(kappa), " ", lambda, "= ", .(lambda), " ", BlockL, "= ", .(BLOCKLENGTH)))) 

dev.off()
