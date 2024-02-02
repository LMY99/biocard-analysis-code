puMVN <- function(mu,V,q,Ms,k=NULL,eps=1e-2,log=T,verbose=FALSE){
  #require(Rmpfr)
  stopifnot(c(length(mu)==dim(V)[1],dim(V)[1]==dim(V)[2],q<=length(mu)))
  require(mvtnorm)
  l <- length(mu)
  if(is.null(k)) k <- 1:(l-q)
  res <- rep(0,length(k)); #res <- mpfr(res, 128)
  #rel_error <- numeric(0)
  for(ki in 1:length(k)){
    # M <- array(0,dim(V)+c(1,0))
    # for(i in 1:(q+1)) M[i,i] <- 1
    # if(q+2<=l){
    #   for(i in (q+2):l){
    #     M[i,i] <- 1; M[i,i-1] <- -1
    #     if(i>q+k0) M[i,] <- M[i,] * (-1)
    #   }}
    # M[l+1,l] <- 1
    k0 <- k[ki]
    M <- Ms[,,k0]
    new.mu <- drop(M%*%mu)
    new.V <- M%*%V%*%t(M); new.V[l+1,l+1] <- new.V[l+1,l+1]*(1+eps)
    #print(new.V[l+1,l+1]-eps)
    lb <- c(rep(-Inf,q),rep(0,l-q+1))
    ub <- rep(+Inf,l+1)
    #point <- c(rep(0.001/k0,k0),rep(0.001/(l+1-q-k0),l+1-q-k0))
    #if(q>0) point <- c(new.mu[1:q],point)
    # point <- ifelse(new.mu>=lb,new.mu,lb)
    #pp <- mvtnorm::dmvnorm(point,new.mu,new.V,log=TRUE) - 
    #TruncatedNormal::dtmvnorm(point,new.mu,new.V,lb,ub,log=TRUE,type='qmc')
    #  tmvtnorm::dtmvnorm(point,new.mu,new.V,lb,ub,log=TRUE) #-Inf occuring here
    #pp <- TruncatedNormal::pmvnorm(new.mu, new.V, lb, ub, type='mc')
    #pp <- mvtnorm::pmvnorm(lb,ub,new.mu,sigma=new.V)
    pp <- logpmvnorm(lb,ub,new.mu,new.V)
    res[ki] <- pp
    #rel_error <- c(rel_error,attr(pp,"rel_error"))
    #res <- c(res, log(pp))
  }
  #if(any(is.na(res))) print(res)
  #print(log(res))
  #if(verbose) print(round(rbind(res,rel_error*100),2))
  #print(res)
  if(log) return((res))
  else return(exp(res))
}
lwmean <- function(a,b,w=1/2){
  require(matrixStats)
  #require(Rmpfr)
  a0=a+log(w);b0=b+log(1-w)
  mat=rbind(a0,b0)
  return(colLogSumExps(mat))
}

logpmvnorm <- function(lb,ub,mu,Sigma,Nmax=1e3){
  require(matrixStats)
  require(VGAM)
  #require(Rmpfr)
  precision=100
  a=lb-mu; b=ub-mu
  #a=mpfr(a,precision);b=mpfr(b,precision)
  # Reordering
  # double_inf=which((a==-Inf)&(b==+Inf))
  # no_inf=which(!is.infinite(b-a))
  # no_inf=no_inf[order((b-a)[no_inf],decreasing=FALSE)]
  # a_inf=which((a==-Inf)&(b!=+Inf))
  # b_inf=which((a!=-Inf)&(b==+Inf))
  # score=rep(NA,length(a))
  # score[a_inf]=b[a_inf]; score[b_inf]=-a[b_inf]
  # one_inf=c(a_inf,b_inf)
  # one_inf=one_inf[order(score[one_inf],decreasing=FALSE)]
  # perm=c(no_inf,one_inf,double_inf)
  # a=a[perm]; b=b[perm]; Sigma=Sigma[perm,perm];
  # Monte-Carlo
  C=t(chol(Sigma))
  m=ncol(C)
  #d=rep(0,m);e=rep(0,m);f=rep(0,m);y=rep(0,m-1)
  d=array(0,dim=c(Nmax,m))
  e=array(0,dim=c(Nmax,m))
  f=array(0,dim=c(Nmax,m))
  y=array(0,dim=c(Nmax,m-1))
  temp=array(0,dim=c(Nmax,m-1))
  #d=mpfr(d,precision);e=mpfr(e,precision);f=mpfr(f,precision)
  #y=mpfr(y,precision);temp=mpfr(temp,precision)
  d[,1]=ifelse(a[1]==-Inf,-Inf,pnorm(a[1]/C[1,1],log=T));
  e[,1]=ifelse(b[1]==+Inf,0,pnorm(b[1]/C[1,1],log=T));
  f[,1]=e[1,1]+log1mexp(e[1,1]-d[1,1])
  
  # for(N in 1:Nmax){
  #   w=runif(m-1)
  #   for(i in 2:m){
  #     #print(c(e[i-1],d[i-1]))
  #     y[i-1]=qnorm(exp(d[i-1])+w[i-1]*(exp(e[i-1])-exp(d[i-1])),log=F)
  #     #if(N==1) print(c(w[i-1],y[i-1]))
  #     d[i]=ifelse(a[i]==-Inf,-Inf,
  #       pnorm((a[i]-sum(C[i,1:(i-1)]*y[1:(i-1)]))/C[i,i],log=T)          
  #     )
  #     e[i]=ifelse(b[i]==+Inf,0,
  #       pnorm((b[i]-sum(C[i,1:(i-1)]*y[1:(i-1)]))/C[i,i],log=T)        
  #     )
  #     #if(N==1) print(c(d[i],e[i]))
  #     f[i]=e[i]+log1mexp(e[i]-d[i])+f[i-1]
  #     #if(N==1) print(f[i])
  #   }
  #   #print(y)
  #   #if(N==1) print("***")
  #   #delta=f[m]+log1p(-exp(Intsum-f[m]))-log(N)
  #   #Intsum=logSumExp(c(Intsum,delta))
  #   ff[N]=f[m]
  # }
  #cat(c(sum(is.na(d[,1])),sum(is.na(e[,1])),sum(is.na(f[,1]))))
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
    #cat(c(sum(is.na(d[,i])),sum(is.na(e[,i])),sum(is.na(f[,i]))))
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
  #cat(sum(is.na(ad))," ",sum(is.na(f[,m])), "\n")
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