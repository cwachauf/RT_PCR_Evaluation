// fit model with 2 consecutive reactions: A + D --> P
// A + DO1 --> A*  (k_on1)
// B + DO2 --> B*  (k_on1)
// A* + B* --> C   (k_on2)
// y[1] = A
// y[2] = B
// y[3] = D01
// y[4] = D02
// y[5] = A*
// y[6] = B*
// y[7] = C
functions
{
  real[] KineticSys(real t,real[] y,real[] theta,real[] x_r,int[] x_i)
  {
    
    real dydt[7];
    dydt[1] <- -theta[1]*y[1]*y[3];
    dydt[2] <- -theta[1]*y[2]*y[4];
    dydt[3] <- -theta[1]*y[1]*y[3];
    dydt[4] <- -theta[1]*y[2]*y[4];
    dydt[5] <- theta[1]*y[1]*y[3] - theta[2]*y[5]*y[6];
    dydt[6] <- theta[1]*y[2]*y[4] - theta[2]*y[5]*y[6];
    dydt[7] <- theta[2]*y[5]*y[6];
    return dydt;
  }
}

data
{
  int<lower=1> num_times; // how many time points ?
  int<lower=1> num_traces; // how many different intensity traces (differing by the initial concentration) ?
  real<lower=0> ts[num_times]; // times (should be the same for all traces...)
  real intens[num_traces,num_times]; // intensity data..
  
  real<lower=0> mu_A0s[num_traces];
  real<lower=0> sigma_A0s[num_traces];
  
  real<lower=0> mu_B0s[num_traces];
  real<lower=0> sigma_B0s[num_traces];
  
  real<lower=0> mu_D1_0s[num_traces];
  real<lower=0> sigma_D1_0s[num_traces];
  
  real<lower=0> mu_D2_0s[num_traces];
  real<lower=0> sigma_D2_0s[num_traces];
  
  real<lower=0> mu_amps_cy35[num_traces];
  real<lower=0> mu_bgs_cy35[num_traces];
}

transformed data
{
  real x_r[0];
  int x_i[0];
}

parameters
{
  real<lower=0> sigma;
  real<lower=0> k_on1;
  real<lower=0> k_on2;
  
  real <lower=0> A0s[num_traces];
  real <lower=0> B0s[num_traces];
  real <lower=0> D01s[num_traces];
  real <lower=0> D02s[num_traces];
  
  real <lower=0> amps_cy35[num_traces];
  real <lower=0> bgs_cy35[num_traces];
}

model
{
  real y_hat[num_traces,num_times,7];
  real conc_inits[num_traces,7];
  real theta[2];
  sigma ~ cauchy(0,2.5); 
 
  k_on1 ~ gamma(0.001,0.001);
  k_on2 ~ gamma(0.001,0.001);
  
  theta[1] <- k_on1;
  theta[2] <- k_on2;
  
  ## all amplitudes and background-values
  for(i in 1:num_traces)
  {
    amps_cy35[i] ~ normal(mu_amps_cy35[i],200);
    bgs_cy35[i] ~ normal(mu_bgs_cy35[i],500);
    A0s[i] ~ normal(mu_A0s[i],sigma_A0s[i]);
    B0s[i] ~ normal(mu_B0s[i],sigma_B0s[i]);
    D01s[i] ~ normal(mu_D1_0s[i],sigma_D1_0s[i]);
    D02s[i] ~ normal(mu_D2_0s[i],sigma_D2_0s[i]);
    
    conc_inits[i,1] <- A0s[i];
    conc_inits[i,2] <- B0s[i];
    conc_inits[i,3] <- D01s[i];
    conc_inits[i,4] <- D02s[i];
    conc_inits[i,5] <- 0.0;
    conc_inits[i,6] <- 0.0;
    conc_inits[i,7] <- 0.0;
    
  }
  for(i in 1:num_traces)
  {
    y_hat[i,] <- integrate_ode(KineticSys,conc_inits[i,],0.0,ts,theta,x_r,x_i);
  }
  for(i in 1:num_traces)
  {
    for(t in 1:num_times)
    {
      ##intens[i,t] <- y_hat[i,t,7]
      y_hat[i,t,7] <- y_hat[i,t,7]*amps_cy35[i] + bgs_cy35[i];
      intens[i,t] ~ normal(y_hat[i,t,7],sigma);
    }
  }
}

