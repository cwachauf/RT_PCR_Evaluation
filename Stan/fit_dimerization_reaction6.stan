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
    
    real ts[num_times]; // times at which model data is given
    real cy35_data[num_traces,num_times];
    
    real<lower=0> init_concs[num_traces,2];
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
    real<lower=0,upper=0.1> k_on; // rate constants (theta[1] and theta[2])
    real<lower=0,upper=0.1> k_off;
    
    real<lower=0,upper=12000> bgs_cy35[num_traces];
    real<lower=0,upper=2000> amps_cy35;
    
  //  real<lower=0,upper=2000> amps_cy35[num_traces];
  }
  
  model
  {
    real y_hat[num_traces,num_times,2];
    real theta[2];
    // specify the priors:
    sigma ~ normal(40,5.0);
   
    k_on ~ gamma(0.01,0.01);
    k_off ~ gamma(0.01,0.01);
    
    theta[1] = k_on;
    theta[2] = k_off;
   
    amps_cy35 ~ normal(1000,100);
    
    bgs_cy35[1] ~ normal(6500,200);
    bgs_cy35[2] ~ normal(7200,200);
    bgs_cy35[3] ~ normal(8000,200);
   ## bgs_cy35[2] ~ normal(7000,100);
    ##for(i in 1:num_traces)
    #{#
    #  bgs_cy35[i] ~ normal(6000,300);
      //amps_cy35[i] ~ normal(800,200);
    #}
    
    for(i in 1:num_traces)
    {
      y_hat[i,] = integrate_ode_bdf(KineticSys,init_concs[i,],t0,ts,theta,x_r,x_i);
      for(t in 1:num_times)
      {
        //y_hat[i,t,2] = bgs_cy35[i] + y_hat[i,t,2]*amps_cy35[i];
        y_hat[i,t,2] = bgs_cy35[i] + y_hat[i,t,2]*amps_cy35;
         cy35_data[i,t] ~ normal(y_hat[i,t,2],sigma);
      }
      //y_hat[i,,2] = bgs_cy35[i] + amps_cy35*y_hat[i,,2];
      //cy35_data[i,] ~ normal(y_hat[i,,2],sigma);
    }
    
  }