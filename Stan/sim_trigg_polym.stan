functions
{
  real[] KineticSys(real t,real[] y,real[] theta,real[] x_r,int[] x_i)
  {
    real dydt[4];
    dydt[1] <- -theta[1]*y[1]*y[2];
    dydt[2] <- -theta[1]*y[1]*y[2];
    dydt[3] <- theta[1]*y[1]*y[2] - theta[2]*y[3]^2;//2*theta[2]*y[3]^2;
    dydt[4] <- theta[2]*y[3]^2;
    return dydt;
  }
}

data
{
  int<lower=1> T; // number of data points
  real y0_1[4]; // initial concentrations part2 
  real y0_2[4]; // initial concentrations part2
  real t0;    // time for the initial concentrations
  real ts[T]; // times at which model sh
  real theta[2]; // the two rate constants....
  real<lower=0> amp_cy35; // amplitude (proportionality constant between intensity and concentration)
  real<lower=0> bg_cy35; // background-signal
  real<lower=0> sigma; // noise
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
  real cy35_1_nf[T];
  real cy35_2_nf[T];
  real y_hat[T,8];
  real y_hat_1[T,4];
  real y_hat_2[T,4];
  
  y_hat_1 <- integrate_ode(KineticSys,y0_1,t0,ts,theta,x_r,x_i);
  y_hat_2 <- integrate_ode(KineticSys,y0_2,t0,ts,theta,x_r,x_i);
  for(t in 1:T)
  {
    ##print(t);
    cy35_1_nf[t] <- bg_cy35 + amp_cy35*y_hat_1[t,4];
    cy35_2_nf[t] <- bg_cy35 + amp_cy35*y_hat_2[t,4];
    
    y_hat[t,1] <- bg_cy35 + amp_cy35*y_hat_1[t,1]+normal_rng(0,sigma);
    y_hat[t,2] <- bg_cy35 + amp_cy35*y_hat_1[t,2]+normal_rng(0,sigma);
    y_hat[t,3] <- bg_cy35 + amp_cy35*y_hat_1[t,3]+normal_rng(0,sigma);
    y_hat[t,4] <- bg_cy35 + amp_cy35*y_hat_1[t,4]+normal_rng(0,sigma);
    
    y_hat[t,5] <- bg_cy35 + amp_cy35*y_hat_2[t,1]+normal_rng(0,sigma);
    y_hat[t,6] <- bg_cy35 + amp_cy35*y_hat_2[t,2]+normal_rng(0,sigma);
    y_hat[t,7] <- bg_cy35 + amp_cy35*y_hat_2[t,3]+normal_rng(0,sigma);
    y_hat[t,8] <- bg_cy35 + amp_cy35*y_hat_2[t,4]+normal_rng(0,sigma);
  }
}