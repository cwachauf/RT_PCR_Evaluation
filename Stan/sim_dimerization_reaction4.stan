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


data
{
  
  int<lower=1> num_times; // number of data points
  int<lower=1> num_traces;
  real ts[num_times];
  
  real<lower=0> k1;
  real<lower=0> k2;
  
  real bgs_cy35[num_traces];
  real amps_cy35[num_traces];
  
  //real<lower=0> A0s[num_traces]; // initial concentrations for A0s (the inactivated monomers)
  //real<lower=0> D0s[num_traces]; // initial concentrations for D0s (the dimerization oligos)
  real<lower=0> init_concs[num_traces,4];
  
  real t0;
}

transformed data
{
  real x_r[0];
  int x_i[0];
}

model
{
}

generated quantities
{
  real cy35_nf[num_traces,num_times];
  real y_hat[num_traces,num_times,4];
  real theta[2];
  
  theta[1] = k1;
  theta[2] = k2;
  
  for(i in 1:num_traces)
  {
    y_hat[i,] = integrate_ode_rk45(KineticSys,init_concs[i,],t0,ts,theta,x_r,x_i);
    for(t in 1:num_times)
    {
      cy35_nf[i,t] = bgs_cy35[i] + amps_cy35[i]*y_hat[i,t,4];
    }
  }
}