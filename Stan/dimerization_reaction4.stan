functions
{
  real[] KineticSys(real t,real[] y,real[] theta,real[] x_r,int[] x_i)
  {
    
    real dydt[4];
    dydt[1] = -theta[1]*y[1]*y[2];
    dydt[2] = -theta[1]*y[1]*y[2];
    dydt[3] = theta[1]*y[1]*y[2] - theta[2]*y[3]^2;
    dydt[4] = theta[2]*y[3]^2;
    return dydt;
  }
}

// simulate some shit...

data
{
  int<lower=1> num_times; // number of data points
  int<lower=1> num_traces;
  real ts[num_times];
  real cy35_data[num_traces,num_times];
  
  real<lower=0> A0s[num_traces]; // initial concentrations for A0s (the inactivated monomers)
  real<lower=0> D0s[num_traces]; // initial concentrations for D0s (the dimerization oligos)
  real t0;
}

transformed data
{
  real x_r[0];
  int x_i[0];
}

parameters
{
  real<lower=0,upper=200> sigma; // noise, supposed to be equal in both cases...
  real<lower=0,upper=0.1> k1; // rate constants (theta[1] and theta[2])
  real<lower=0,upper=0.1> k2;
  real<lower=0,upper=10000> bgs_cy35[num_traces];
  real<lower=0,upper=3000> amps_cy35[num_traces];
}

model
{
  real y_hat[num_traces,num_times,4];
  real concs0[num_traces,4];
  real theta[2];
  
  sigma ~ cauchy(0,2.5);
  k1 ~ gamma(0.001,0.001);
  k2 ~ gamma(0.001,0.001);
  for(i in 1:num_traces)
  {
    amps_cy35[i] ~ normal(800,200);
    bgs_cy35[i] ~ normal(5000,500);
    concs0[i,1] = A0s[i];
    concs0[i,2] = D0s[i];
    concs0[i,3] = 0.0;
    concs0[i,4] = 0.0;
  }
  theta[1] = k1;
  theta[2] = k2;
  for(i in 1:num_traces)
  {
    y_hat[i,] = integrate_ode_rk45(KineticSys,concs0[i,],t0,ts,theta,x_r,x_i);
    for(t in 1:num_times)
    {
      y_hat[i,t,4] = bgs_cy35[i] + y_hat[i,t,4]*amps_cy35[i];
      cy35_data[i,t] ~ normal(y_hat[i,t,4],sigma);
    }

  }
}