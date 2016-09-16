// model system with 
// activation step and back-reaction
// A + B --> C   (with rate constant k_on,1)
// C + C <--> D (with rate constant k_on,2 and k_off)
// leading to the following system of equations:
  // dA/dt = - k_on,1 * A*B
  // dB/dt = -k_on,1*A*B
  // dC/dt = k_on,1*A*B - k_on,2*C^2 + k_off,2*D
  // dD/dt = k_on,2*C^2 - k_off,2*D
  
functions
{
  real[] KineticSys(real t,real[] y,real[] theta,real[] x_r,int[] x_i)
  {
    real dydt[4];
    dydt[1] = -theta[1]*y[1]*y[2];
    dydt[2] = -theta[1]*y[1]*y[2];
    dydt[3] = theta[1]*y[1]*y[2] - theta[2]*y[3]^2 + theta[3]*y[4];
    dydt[4] = theta[2]*y[3]^2 - theta[3]*y[4];
    return dydt;
  } 
} 
  
data
{
  int<lower=1> num_times; // number of data points
  int<lower=1> num_traces;
  
  real ts[num_times]; // times at which model data is given
  real cy35_data[num_traces,num_times];
    
  real<lower=0> init_concs[num_traces,4];
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
  real<lower=0,upper=0.1> k_on1; // rate constants (theta[1] and theta[2])
  real<lower=0,upper=0.008> k_on2;
  real<lower=0,upper=0.008> k_off2;
  
  real<lower=0,upper=12000> bgs_cy35[num_traces];
  real<lower=0,upper=2000> amps_cy35[num_traces];
}
  
model
{
  real y_hat[num_traces,num_times,4];
  real theta[3];
  // specify the priors:
  sigma ~ normal(40,5);
  k_on1 ~ gamma(0.001,0.001);
  k_on2 ~ gamma(0.001,0.001);
  k_off2 ~ gamma(0.001,0.001);
  
  theta[1] = k_on1;
  theta[2] = k_on2;
  theta[3] = k_off2;
  for(i in 1:num_traces)
  {
    bgs_cy35[i] ~ normal(6000,1000);
    amps_cy35[i] ~ normal(800,200);
  }
    
  for(i in 1:num_traces)
  {
    y_hat[i,] = integrate_ode_rk45(KineticSys,init_concs[i,],t0,ts,theta,x_r,x_i);
    for(t in 1:num_times)
    {
      y_hat[i,t,4] = bgs_cy35[i] + y_hat[i,t,4]*amps_cy35[i];
      cy35_data[i,t] ~ normal(y_hat[i,t,4],sigma);
    }
  }
    
}
  