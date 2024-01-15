filename <- "biocard_result"
load(sprintf("%s.RData",filename))
filename <- "biocard_result"
usePackage("splines2")
usePackage("TruncatedNormal")
library(tidyverse)

library(ggplot2)
pdf("visual.pdf",width=7)
for(biomarker in colnames(Y)){
  print(ggplot(data=df)+
          geom_line(aes_string(x="ageori",y=biomarker,group="Study.ID"),alpha=0.3)+
          geom_smooth(aes_string(x='ageori',y=biomarker),formula=y~x,
                      method='loess',alpha=0.5,color='red'))
}
dev.off()

indice <- (Burnin+1):R

# Make MCMC trace plot to check for convergence ----

pdf(sprintf("TracePlot_%s.pdf",filename),width=14)
plot(indice,sigmays[indice],type='l',xlab='index',ylab='Error Variance')
plot(indice,sigmaws[indice],type='l',xlab='index',ylab='Random Effect Variance')
matplot(indice,t(coefs[1:nX,1,indice]),type='l',xlab='index',ylab='Coefficients of X(Biomarker 1)')
matplot(indice,log(t(coefs[(nX+1):(nX+dfi-4),2,indice])),type='l',xlab='index',ylab='Log-Coefficients of Spline(Biomarker 2)')
matplot(indice,log(t(pens[,indice])),type='l',xlab='index',ylab='Log Penalty Parameters')
plot(indice,REs[1,1,indice],type='l',xlab='index',ylab='Random Intercept: Individual 1 Biomarker 1')
dev.off()

indice <- (Burnin+1):R

# Calculate Fitted Values ----
# We use the posterior mean for X-coefficients, 
# random effects, and spline coefficients

fit.coef <- apply(coefs[,,indice],c(1,2),mean)
fit.RE <- apply(REs[,,indice],c(1,2),mean)
fit.Y <- array(0, dim(Y))
for(k in 1:ncol(Y)) fit.Y[,k] <- drop(covar.list[[k]]%*%fit.coef[,k]) + fit.RE[df$ID,k]
# If an individual has no measurement of one biomarker entirely
# then the corresponding fitted value will be discarded
#for(j in 1:K){
#  empty <- which(long_ss[,j]==0)
#  fit.Y[df$ID %in% empty,j] <- NA
#}
fit.Y[is.na(Y)] <- NA
colnames(fit.Y) <- paste(colnames(Y),"_fit",sep='')

# Draw goodness-of-fit graph

frame <- cbind(ID=df$ID,age=df$ageori,Y,fit.Y)
library(ggplot2)
pdf(sprintf("Goodness_of_fit_%s.pdf",filename),width=7)
for(biomarker in colnames(Y)){
  print(ggplot(data=frame)+
          geom_line(aes_string(x="age",y=biomarker,group="ID"),alpha=0.3)+
          geom_line(aes_string(x="age",y=paste(biomarker,"_fit",sep=''),group="ID"),
                    alpha=0.3,color='blue') +
          geom_smooth(aes_string(x='age',y=biomarker),formula=y~x,
                      method='loess',alpha=0.5,color='red'))
}
dev.off()

# Draw fitted spline curves
# fit.spline <- fit.coef[(nX+1):(nX+dfi),]
# tt <- seq(min(df$ageori),max(df$ageori),length.out=10000)
# tt.Spline <- iSpline(tt,knots = attr(B,"knots"),
#                      Boundary.knots = attr(B,"Boundary.knots"),
#                      degree = 2)
# yy <- tt.Spline %*% fit.spline
tt <- seq(0,120,length.out=10000)
yy <- array(0,c(length(tt),K))
yy_std <- yy
for(j in 1:K){
  tsp <- ibs(tt,knots=knot.list[[j]],Boundary.knots=boundary.knot,
             degree=2, intercept=TRUE)
  tsp <- tsp[,3:(ncol(tsp)-2)]
  yy[,j] <- tsp %*% fit.coef[(nX+1):(nX+dfi-4),j]
  yy_std[,j] <- yy[,j]/max(yy[,j])
}
colnames(yy) <- colnames(Y)
frame.spline <- cbind(age=tt,yy)
pdf(sprintf("Spline_%s.pdf",filename),width=14)
frame.spline <- as.data.frame(frame.spline)
# p <- ggplot(frame.spline)
# for(j in 1:K){
#   p <- p + geom_line(aes_string(x='age',y=colnames(Y)[j],
#                                 color=colnames(Y)[j]))
# }
# p <- p + scale_color_manual(values=c("MMSCORE"="darkred","logmem"="red","DSST"="orange","biec.thik"="yellow",
#                                      "Hippo_dadjust"="green","Ent_dadjust"="darkgreen","MTL1"="lightblue",
#                                      "SPARE_AD"="blue","ttau"="darkblue","ptau181"="purple","AB42AB40"="black"),
#                             labels=colnames(Y),
#                             guide="legend")
# p <- p + labs(x = "Age", y = "Biomarker Abnormality",
#               title = "Adjusted Biomarker Trajectory") +
#   #scale_y_continuous(breaks=seq(0,1,by=0.2)) +
#   theme(plot.title = element_text(hjust = 0.5))
# p <- p + geom_vline(xintercept=range(df$ageori))

colors0 <- c("darkred","red","orange","yellow","green","darkgreen","lightblue",
             "blue","darkblue","purple","black")
names(colors0) <- colnames(Y)

p <- frame.spline %>% gather("biomarker","Y",-age) %>%
  ggplot(aes(x=age,y=Y)) + geom_line(aes(color=biomarker)) +
  scale_color_manual(values=colors0)
print(p)
dev.off()

colnames(yy_std) <- colnames(Y)
frame.spline <- cbind(age=tt,yy_std)
pdf(sprintf("SplineStd_%s.pdf",filename),width=14)
frame.spline <- as.data.frame(frame.spline)
# p <- ggplot(frame.spline)
# for(j in 1:K){
#   p <- p + geom_line(aes_string(x='age',y=colnames(Y)[j],
#                                 color=as.factor(j)))
# }
# p <- p + scale_color_manual(values=c("1"="darkred","2"="red","3"="orange","4"="yellow",
#                                      "5"="green","6"="darkgreen","7"="lightblue",
#                                      "8"="blue","9"="darkblue","10"="purple","11"="black"),
#                             labels=colnames(Y),
#                             guide="legend")
# p <- p + labs(x = "Age", y = "Biomarker Abnormality",
#               title = "Adjusted Biomarker Trajectory") +
#   scale_y_continuous(breaks=seq(0,1,by=0.2)) +
#   theme(plot.title = element_text(hjust = 0.5))
# p <- p + geom_vline(xintercept=range(df$ageori))
# print(p)

p <- frame.spline %>% gather("biomarker","Y",-age) %>%
  ggplot(aes(x=age,y=Y)) + geom_line(aes(color=biomarker)) +
  scale_color_manual(values=colors0)
print(p)
dev.off()


# # Calculate age point where second derivative reaches zero
# pdf(sprintf("Inflextion Age_%s.pdf",filename))
# tt <- seq(0,120,length.out=10000)
# inflex <- array(0,c(K,length(indice)))
# layout(matrix(c(1,2,3,4),2,2))
# for(j in 1:K){
#   dsp <- dbs(tt,derivs=1,knots=knot.list[[j]],Boundary.knots=boundary.knot,
#              degree=2, intercept=TRUE)
#   dsp <- dsp[,3:(ncol(dsp)-2)]
#   darray <- dsp%*%coefs[(nX+1):(nX+dfi-4),j,indice]
#   u=apply(darray,2,function(x)(max(which(x>0))))
#   u[u==-Inf]=1
#   inflex[j,]=tt[u]
# 
#   plot(density(inflex[j,]),xlab='age',main=sprintf("%s: Mean=%.2f,Median=%.2f",colnames(Y)[j],
#                                           mean(inflex[j,]),median(inflex[j,])))
#   #yy[,j] <- tsp %*% fit.coef[(nX+1):(nX+dfi),j]
#   #yy[,j] <- yy[,j]/max(yy[,j])
# }
# dev.off()
