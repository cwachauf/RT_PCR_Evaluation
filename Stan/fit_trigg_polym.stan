functions
{
  real[] KineticSys(real t,real[] y,real[] theta,real[] x_r,int[] x_i)
  {
    real dydt[4];
    dydt[1] <- -theta[1]*y[1]*y[2];
    dydt[2] <- -theta[1]*y[1]*y[2];
    dydt[3] <- theta[1]*y[1]*y[2] - 2*theta[2]*y[3]^2;
    dydt[4] <- theta[2]*y[3]^2;
    return dydt;
  }
}

data
{
  int<lower=1> T; // number of data points
  real cy35_1[T]; // cy35-intensity with initial concentrations part1 
  real cy35_2[T]; // cy35-intensity with initial concentrations part2 
  real ts[T]; // times at which intensities were measured...
  real a0_1;  // initial concentrations 
  real a0_2;  // initial concentrations
  real b0_1;  // initial concentrations
  real b0_2;  // initial concentrations
  real t0;    // time for the initial concentrations
}

transformed data
{
  real x_r[0];
  int x_i[0];
}

parameters
{
  real<lower=0> sigma; // noise, supposed to be equal in both cases...
  real<lower=0> k1; // rate constants (theta[1] and theta[2])
  real<lower=0> k2;
  real<lower=0> amp_cy35;
  real<lower=0> bg_cy35;
}

model
{
  real y_hat_1[T,4];
  real y_hat_2[T,4];
  real conc0_1[4];
  real conc0_2[4];
  real theta[2];
  
  sigma ~ cauchy(0,2.5);
  k1 ~ gamma(0.001,0.001);
  k2 ~ gamma(0.001,0.001);
  amp_cy35 ~ normal(800,100);
  bg_cy35 ~ normal(7300,100);
  
  theta[1] <- k1;
  theta[2] <- k2;
  
  conc0_1[1] <- a0_1;
  conc0_1[2] <- b0_1;
  conc0_1[3] <- 0;
  conc0_1[4] <- 0;
  
  conc0_2[1] <- a0_2;
  conc0_2[2] <- b0_2;
  conc0_2[3] <- 0;
  conc0_2[4] <- 0;
  
  y_hat_1 <- integrate_ode(KineticSys,conc0_1,t0,ts,theta,x_r,x_i);
  y_hat_2 <- integrate_ode(KineticSys,conc0_2,t0,ts,theta,x_r,x_i);
  
  for(t in 1:T)
  {
    y_hat_1[t,4] <- bg_cy35+amp_cy35*y_hat_1[t,4];
    y_hat_2[t,4] <- bg_cy35+amp_cy35*y_hat_2[t,4];
  }
  for(t in 1:T)
  {
    cy35_1[t] ~ normal(y_hat_1[t,4],sigma);
    cy35_2[t] ~ normal(y_hat_2[t,4],sigma);
  }
}