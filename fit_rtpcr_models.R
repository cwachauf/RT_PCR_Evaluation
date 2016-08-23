
Fit_Stan_Trigg_Polym <- function(time_seq,cy35_1,cy35_2,a0_1,b0_1,a0_2,b0_2,num_iter)
{
  require("rstan")
  data_stan_trigg_polym <- list(T=length(time_seq),cy35_1=cy35_1,cy35_2=cy35_2,
                              ts=time_seq,a0_1=a0_1,a0_2=a0_2,b0_1=b0_1,b0_2=b0_2,t0=0)

  init_stan_trigg_polym <- list(list(sigma=50,k1=1e-4,k2=5e-5,amp_cy35=1000,bg_cy35=6500))
  sf_trigg_polym <- stan(file="/Users/christianwachauf/Documents/Scripts/GithubRepo/RT_PCR_Evaluation/Stan/fit_trigg_polym.stan",data=data_stan_trigg_polym,
                       init=init_stan_trigg_polym,iter=num_iter,chain=1)
  return(sf_trigg_polym)
}

Test_RTPCR_Models <- function(df_rt_data)
{
  ## do some shit:
  times <- df_rt_data$time_seq[seq(from=1,to=length(df_rt_data$time_seq),by=3)]
  cy35_1 <- df_rt_data$`10hb_6nM_po_Cy35`[seq(from=1,to=length(df_rt_data$time_seq),by=3)]
  cy35_2 <- df_rt_data$`10hb_12nM_po_Cy35`[seq(from=1,to=length(df_rt_data$time_seq),by=3)]
  
  sf_trigg_polym <- Fit_Stan_Trigg_Polym(times,cy35_1,cy35_2,6,6,6,12,5000)
  summary(sf_trigg_polym)
  return(sf_trigg_polym)
}

Plot_Stan_Trigg_Polym_Mean_Posterior <- function(sf_trigg_polym)
{
  mat <- as.matrix(sf_trigg_polym)
  sigma_mean <- mean(mat[,1])
  k1_mean <- mean(mat[,2])
  k2_mean <- mean(mat[,3])
  amp_cy35_mean <- mean(mat[,4])
  bg_cy35_mean <- mean(mat[,5])
  times <- seq(From=1,to=30000,by=1)
  npnts <- length(times)
  #data_stan_sim_trigg_polym <- list(T=npnts,y0_1=c(6,6,0,0),y0_2=c(6,12,0,0),t0=0.0,ts=times,theta=c(k1_mean,k2_mean),
   #                                 amp_cy35=amp_cy35_mean,bg_cy35=bg_cy35_mean,sigma=sigma_mean)
  #stan_sim_trigg_polym <- stan(file="/Users/christianwachauf/Documents/Scripts/GithubRepo/RT_PCR_Evaluation/Stan/sim_trigg_polym.stan",
   #                            data=data_stan_sim_trigg_polym,iter=1,chain=1,algorithm="Fixed_param")
  out1 <- PlotKineticModel(k1_mean,k2_mean,6,12,0,0,bg_cy35_mean,amp_cy35_mean)
  return(out1)
}

PlotKineticModel <- function(k1_par,k2_par,A0,B0,C0,D0,bg_cy35,amp_cy35)
{
  ## set value of parameters
  parameters <- c(k1=k1_par,k2=k2_par)	
  ## set initial concentration values
  state <- c(A=A0,B=B0,C=C0,D=D0)
  
  ## set kinetic system...
  KineticSys<-function(t,state,parameters){
    with(as.list(c(state,parameters)),{
      dA <- -k1*A*B
      dB <- -k1*A*B
      dC <- k1*A*B - k2*C^2;#2*k2*C^2
      dD <- k2*C^2;
      list(c(dA,dB,dC,dD))
    })
  }
  
  ## time-vector:
  time <- seq(0,40000,by=1)
  
  ## solve for different initial concentrations
  library(deSolve)
  
  out <- ode(y=state,times=time,func=KineticSys,parms=parameters)
  res <- bg_cy35+out[,5]*amp_cy35;
  time <- out[,1]
  df_res <- data.frame(time=time,cy35=res)
  return(df_res)
}

#real dydt[4];
#dydt[1] <- -theta[1]*y[1]*y[2];
#dydt[2] <- -theta[1]*y[1]*y[2];
#dydt[3] <- theta[1]*y[1]*y[2] - 2*theta[2]*y[3]^2;
#dydt[4] <- theta[2]*y[3]^2;#
#data
#{
#  int<lower=1> T; // number of data points
#  real y0_1[4]; // initial concentrations part2 
#  real y0_2[4]; // initial concentrations part2
#  real t0;    // time for the initial concentrations
#  real ts[T]; // times at which model sh
#  real theta[2]; // the two rate constants....
#  real<lower=0> amp_cy35; // amplitude (proportionality constant between intensity and concentration)
#  real<lower=0> bg_cy35; // background-signal
##  real<lower=0> sigma; // noise
#}