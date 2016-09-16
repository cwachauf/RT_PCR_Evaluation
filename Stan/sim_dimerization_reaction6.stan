// model system with 
// activation step and back-reaction
// C + C <--> D (with rate constant k_on,2 and k_off)
// leading to the following system of equations:
  // 
  
  functions
{
  real[] KineticSys(real t,real[] y,real[] theta,real[] x_r,int[] x_i)
  {
    real dydt[2];
    dydt[1] = -theta[1]*y[1]^2 + theta[2]*y[2];
    dydt[2] = theta[1]*y[1]^2 - theta[2]*y[2];
    return dydt;
  } 
} 

data
{
  int<lower=1> num_times; // number of data points
  int<lower=1> num_traces;
  
  real bgs_cy35[num_traces];
  real amps_cy35[num_traces];
  
  real k_on;
  real k_off;
  
  real ts[num_times]; // times at which model data is given
  real<lower=0> init_concs[num_traces,2];
  real t0;
}

transformed data
{
  real x_r[0];
  int x_i[0];
  real theta[2];
  theta[1] = k_on;
  theta[2] = k_off;
}

model
{
}

generated quantities
{
  real cy35_nf[num_traces,num_times];
  real y_hat[num_traces,num_times,2];
  
  for(i in 1:num_traces)
  {
    y_hat[i,] = integrate_ode_rk45(KineticSys,init_concs[i,],t0,ts,theta,x_r,x_i);
    for(t in 1:num_times)
    {
      cy35_nf[i,t] = bgs_cy35[i] + amps_cy35[i]*y_hat[i,t,2];
    }
  }
}
