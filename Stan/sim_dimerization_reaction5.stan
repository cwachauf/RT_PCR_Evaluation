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
  
  real<lower=0> k_on1;
  real<lower=0> k_on2;
  real<lower=0> k_off2;

  real<lower=0> amps_cy35[num_traces]; // amplitude (proportionality constant between intensity and concentration)
  real<lower=0> bgs_cy35[num_traces]; // background-signal
  real<lower=0> init_concs[num_traces,4];
  real t0;
}

transformed data
{
  real x_r[0];
  int x_i[0];
  real theta[3];
  theta[1] = k_on1;
  theta[2] = k_on2;
  theta[3] = k_off2;
}

model
{
}

generated quantities
{
  real cy35_nf[num_traces,num_times];
  real y_hat[num_traces,num_times,4];
  
  for(i in 1:num_traces)
  {
    y_hat[i,] = integrate_ode_rk45(KineticSys,init_concs[i,],t0,ts,theta,x_r,x_i);
    for(t in 1:num_times)
    {
      cy35_nf[i,t] = bgs_cy35[i] + amps_cy35[i]*y_hat[i,t,4];
    }
  }
}

