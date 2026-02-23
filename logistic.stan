data {
  int<lower=1> Nobs;
  int<lower=1> Npred;
  int<lower=1> Nout;
  int<lower=1> Nind;
  matrix[Nobs,Npred] x;
  vector[Nobs] t;
  real y[Nobs,Nout];
  array[Nobs,Nout] int<lower=0,upper=1> y_obs;
  array[Nobs] int<lower=1,upper=Nind> jj;
}
parameters {
  matrix[Nout,Npred] beta;
  real<lower=0> lscale[Nout];
  real<lower=0> lpos[Nout];
  real<lower=0> lamp[Nout];
  real<lower=0> sigmarandom;
  real<lower=0> sigmaerror;
  matrix[Nout,Nind] randomint;
}
model {
  real mu[Nobs,Nout];
  for(n in 1:Nout){
    for(k in 1:Npred){
      beta[n,k] ~ normal(0, sqrt(10000));
    }
  }
  sigmarandom ~ inv_gamma(3,0.5);
  sigmaerror ~ inv_gamma(3,0.5);
  lscale ~ normal(5, 1);
  lpos ~ normal(70, 30);
  lamp ~ normal(2, 1);
  for(n in 1:Nout){
    for(j in 1:Nind){
      randomint[n,j] ~ normal(0,sqrt(sigmarandom));
    }
  }
  for (n in 1:Nobs){
    for (k in 1:Nout){
      mu[n,k] = randomint[k,jj[n]]+inv_logit((t[n]-lpos[k])/lscale[k])*lamp[k]+dot_product(x[n],beta[k]);
      if (y_obs[n,k]){
        y[n,k] ~ normal(mu[n,k],sqrt(sigmaerror));
      }
    }
  }
}
