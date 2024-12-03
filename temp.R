puMVN <- function(mu,V,q,Ms,k=NULL,eps=1e-2,log=T,verbose=FALSE){
  stopifnot(c(length(mu)==dim(V)[1],dim(V)[1]==dim(V)[2],q<=length(mu)))
  require(mvtnorm)
  l <- length(mu)
  if(is.null(k)) k <- 1:(l-q)
  res <- rep(0,length(k))
  for(ki in 1:length(k)){
    k0 <- k[ki]
    M <- Ms[,,k0]
    new.mu <- drop(M%*%mu)
    new.V <- M%*%V%*%t(M); new.V[l+1,l+1] <- new.V[l+1,l+1]*(1+eps)
    lb <- c(rep(-Inf,q),rep(0,l-q+1))
    ub <- rep(+Inf,l+1)
    pp <- logpmvnorm(lb,ub,new.mu,new.V)
    res[ki] <- pp
  }
  if(log) return((res))
  else return(exp(res))
}
lwmean <- function(a,b,w=1/2){
  require(matrixStats)
  a0=a+log(w);b0=b+log(1-w)
  mat=rbind(a0,b0)
  return(colLogSumExps(mat))
}

logpmvnorm <- function(lb,ub,mu,Sigma,Nmax=1e3){
  require(matrixStats)
  require(VGAM)
  precision=100
  a=lb-mu; b=ub-mu
  C=t(chol(Sigma))
  m=ncol(C)
  d=array(0,dim=c(Nmax,m))
  e=array(0,dim=c(Nmax,m))
  f=array(0,dim=c(Nmax,m))
  y=array(0,dim=c(Nmax,m-1))
  temp=array(0,dim=c(Nmax,m-1))
  d[,1]=ifelse(a[1]==-Inf,-Inf,pnorm(a[1]/C[1,1],log=T));
  e[,1]=ifelse(b[1]==+Inf,0,pnorm(b[1]/C[1,1],log=T));
  f[,1]=e[1,1]+log1mexp(e[1,1]-d[1,1])
  w=array(runif((m-1)*Nmax),dim=c(Nmax,m-1))
  for(i in 2:m){
    temp[,i-1]=lwmean(d[,i-1],e[,i-1],w[,i-1])
    temp[temp[,i-1]>=0,i-1]=0
    y[,i-1]=qnorm(temp[,i-1],log=T)
    y[y[,i-1]==+Inf,i-1]=1e3
    y[y[,i-1]==-Inf,i-1]=-1e3
    if(a[i]==-Inf)
      d[,i]=-Inf
    else
      d[,i]=pnorm((a[i]-colSums(C[i,1:(i-1)]*t(y[,1:(i-1)])))/C[i,i],log=T)
    if(b[i]==+Inf)
      e[,i]=0
    else
      e[,i]=pnorm((b[i]-colSums(C[i,1:(i-1)]*t(y[,1:(i-1)])))/C[i,i],log=T)
    diff=log1mexp(e[,i]-d[,i])
    diff=ifelse(diff!=-Inf,diff,log1mexp(.Machine$double.xmin))
    f[,i]=e[,i]+diff+f[,i-1]
  }
  res=logSumExp(f[,m],na.rm=TRUE)-log(sum(!is.na(f[,m])))
  if(is.na(res)||res==-Inf){ 
    res=(-Inf)
    View(f); View(d); View(e); View(y); View(temp)
    stop("Error")
    attr(res,"rel_error")=0
    return(res)
  }
  else{
    v=cbind(f[,m],res)
    idx=v[,1]<v[,2]
    idx[is.na(idx)]=FALSE
    v[idx,]=v[idx,c(2,1)]
    ad=v[,1]+log1mexp(v[,1]-v[,2])
    ad=ad*2
    var0=logSumExp(ad,na.rm=TRUE)-log(sum(!is.na(ad)))
    var0=var0-log(sum(!is.na(ad)))
    sd0=var0/2
    attr(res,"rel_error")=exp(sd0-res)
  }
  return(res)
}
penalty_Matrix <- function(dimension,smooth.sigma=Inf,flat.sigma=Inf,weight=NULL){
  if(is.null(weight)){
    weight <- (0:(dimension-1))*((dimension-1):0)
  }
  weight <- weight + 0.001
  weight <- weight/max(weight)
  prec <- matrix(0,dimension,dimension)
  for(j in 1:(dimension-1)){
    u <- rep(0,dimension); u[j] <- 1; u[j+1] <- -1;
    prec <- prec + u%*%t(u)/(smooth.sigma/dimension/dimension)
  }
  prec <- prec + diag(1/flat.sigma/weight)
  V <- chol2inv(chol(prec))
  return(list(V=V,prec=prec))
}