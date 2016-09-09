// stan - model: simplest possible model for the
// dimerization reaction (A+A --> P), with rate constant k_on

data
{
  int<lower=1> npnts;
  real<lower=0> ts[npnts];
  real cy35_1[npnts];
  real<lower=0> mu_A0; // mean of the initial concentration
  real<lower=0> sigma_A0; // stdev of the initial concentration
}

parameters
{
  real<lower=0> sigma;
  real<lower=0> k_on;
  real<lower=0> A_0;
  real<lower=0> amp_cy35;
  real<lower=0> bg_cy35;
}

model
{
  real ideal_intens[npnts];
  A_0 ~ normal(mu_A0,sigma_A0);
  k_on ~ gamma(0.001,0.001);
  sigma ~ cauchy(0,2.5);
  amp_cy35 ~ normal(800,200);
  bg_cy35 ~ normal(5000,500);
  for(i in 1:npnts)
  {
    ideal_intens[i] <- bg_cy35 + amp_cy35*(A_0 - 1.0/(k_on*ts[i]+1.0/A_0));
    cy35_1[i] ~ normal(ideal_intens[i],sigma);
  }
}