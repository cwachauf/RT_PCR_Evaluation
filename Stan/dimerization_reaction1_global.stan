// stan - model: simplest possible model for the
// dimerization reaction (A+A --> P), with rate constant k_on

data
{
  int<lower=1> num_times; // how many time points ?
  int<lower=1> num_traces; // how many different intensity traces (differing by the initial concentration) ?
  real<lower=0> ts[num_times]; // times (should be the same for all traces...)
  real intens[num_traces,num_times]; // intensity data..
  
  real<lower=0> mu_A0s[num_traces]; // means for the prior distribution of each concentration
  real<lower=0> sigma_A0s[num_traces]; // stdevs for the prior distribution of each concentration
  
  real<lower=0> mu_amps_cy35[num_traces];
  real<lower=0> mu_bgs_cy35[num_traces];
}

parameters
{
  real<lower=0> sigma;
  real<lower=0> k_on;
  
  real <lower=0> A0s[num_traces];
  real <lower=0> amps_cy35[num_traces];
  real <lower=0> bgs_cy35[num_traces];
}

model
{
  real ideal_intens[num_traces,num_times];
  
  for(i in 1:num_traces)
  {
    A0s[i] ~ normal(mu_A0s[i],sigma_A0s[i]);
    amps_cy35[i] ~ normal(mu_amps_cy35[i],200);
    bgs_cy35[i] ~ normal(mu_bgs_cy35[i],500);
  }
  
  for(i in 1:num_traces)
  {
    for(t in 1:num_times)
    {
      ideal_intens[i,t] <- bgs_cy35[i] + amps_cy35[i]*(A0s[i] - 1.0/(k_on*ts[t]+1.0/A0s[i]));
      intens[i,t] ~ normal(ideal_intens[i,t],sigma);
    }
  }
}