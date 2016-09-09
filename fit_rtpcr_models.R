
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

## Stan_Fit_Simple_Dimerization_Reaction1(times,cy35_data,A0_mu,A0_sigma)
## fits a simple dimerization reaction ("A+A --> D")
Stan_Fit_Simple_Dimerization_Reaction1 <- function(times,cy35_data,A0_mu,A0_sigma,num_iters=200)
{
  require("rstan")
  data_stan_fit_dimerization1 <- list(npnts=length(times),ts=times,cy35_1=cy35_data,mu_A0=A0_mu,sigma_A0=A0_sigma)
  
  init_stan_fit_dimerization1 <- list(list(sigma=40,k_on=1e-4,A_0=4,amp_cy35=750,bg_cy35=5000))
  
  sf_dimerization_reaction1 <- stan(file="/Users/christianwachauf/Documents/Scripts/GithubRepo/RT_PCR_Evaluation/Stan/dimerization_reaction1.stan",data=data_stan_fit_dimerization1,
                                    init=init_stan_fit_dimerization1,iter=num_iters,chain=1)
  return(sf_dimerization_reaction1)
}


## obtain mean values from the fit and
## plot the function with these values:
Plot_Mean_Posterior_Simple_Dimerization_Reaction1 <- function(sf_dimerization_reaction1,times)
{
  mat <- as.matrix(sf_dimerization_reaction1)
  k_on_mean <- mean(mat[,2])
  A_0_mean <- mean(mat[,3])
  amp_cy35_mean <- mean(mat[,4])
  bg_cy35_mean <- mean(mat[,5])
  int_cy35 <- array(0,dim=c(length(times)))
  print(k_on_mean)
  print(A_0_mean)
  for(i in 1:length(times))
  {
    int_cy35[i] <- Dimer_Concentration_Simple_Dimerization1(A_0_mean,k_on_mean,times[i]) 
    int_cy35[i] <- int_cy35[i]*amp_cy35_mean + bg_cy35_mean
  }
  df_int <- data.frame(times,int_cy35)
  return(df_int)
}

## Dimer_Concentration_Simple_Dimerization1(A_0,k_on,t)
## return the dimer concentration for the simple
## reaction: A + A --> D, with rate constant k_on and
## initial concentration A_0, at time t
Dimer_Concentration_Simple_Dimerization1 <- function(A_0,k_on,t)
{
  result <- A_0 - 1.0/(k_on*t+1.0/A_0)
  return(result)
}

## little helper function that neeeds to be adjusted for the individual plots
Make_Plot_With_Data_And_Fit <- function(data_rtpcr,sf_dimerization_reaction1)
{
   time_data <- data_rtpcr$time_seq
   cy35_data <- data_rtpcr$`24hb_4nM_40nM_pos_Cy35`
   plot(time_data,cy35_data,xlim=c(0,10000),xlab="time [s]",ylab="intensity [a.u.]")
   time_sequ <- seq(from=0,to=10000,by=5)
   df_plot_data <- Plot_Mean_Posterior_Simple_Dimerization_Reaction1(sf_dimerization_reaction1,time_sequ)
   points(df_plot_data$times,df_plot_data$int_cy35,type="l",col="red")
   
}

## Return_Prediction_Intervals_For_Posterior_Dimerization_Reaction1(sf_dimerization_reaction,times)
Return_Prediction_Intervals_For_Posterior_Dimerization_Reaction1 <- function(sf_dimerization_reaction,times)
{
  require("coda")
  mat <- as.matrix(sf_dimerization_reaction)
  num_samples <- nrow(mat)
  
  predictions_intensity <- array(0,dim=c(num_samples,length(times)))
  
  int_mean <- array(0,dim=c(length(times)))
  int_lp <- array(0,dim=c(length(times)))
  int_up <- array(0,dim=c(length(times)))
  
  ## loop through the rows:
  for(i in 1:num_samples)
  {
    ## obtain parameters:
    sigma_curr <- mat[i,1]
    k_on_curr <- mat[i,2]
    A_0_curr <- mat[i,3]
    amp_cy35_curr <- mat[i,4]
    bg_cy35_curr <- mat[i,5]
    
    int_cy35 <- array(0,dim=c(length(times)))
    if(i%%1000==1)
    {
      print("iteration: ")
      print(i)
    }
    for(j in 1:length(times))
    {
      predictions_intensity[i,j] <- Dimer_Concentration_Simple_Dimerization1(A_0_curr,k_on_curr,times[j])
      predictions_intensity[i,j] <- predictions_intensity[i,j]*amp_cy35_curr + bg_cy35_curr
      ## add noise (i)
      predictions_intensity[i,j] <- predictions_intensity[i,j] + rnorm(n=1,mean=0,sd=sigma_curr)
    }  
  }
  ## for each time point:
  for(j in 1:length(times))
  {
    int_mean[j] <- mean(predictions_intensity[,j])
    hdi_values <- HPDinterval(as.mcmc(predictions_intensity[,j]))
    int_lp[j] <- hdi_values[1];
    int_up[j] <- hdi_values[2];
  }
  df_preds <- data.frame(times,int_mean,int_lp,int_up)
  return(df_preds)
}

Plot_Data_Dimerization_Reaction1_With_Preds <- function(times,cy35_data,df_preds)
{
  plot(times,cy35_data,xlab="time [s]",ylab="intensity [a.u.]")
  ## add prediction bands
  points(df_preds$times,df_preds$int_mean,type="l",lwd=1,col="red")
  points(df_preds$times,df_preds$int_lp,type="l",lwd=1,col="blue")
  points(df_preds$times,df_preds$int_up,type="l",lwd=1,col="blue")
}

Get_Gamma_Distribution_Parameters_From_Moments <- function(mean,stdev)
{
  ## shape parameter:
  k <- (mean^2)/(stdev^2)
  ## scale parameter:
  theta <- (stdev^2)/mean
  df_gamma_dist_params <- data.frame(k,theta)
  return(df_gamma_dist_params)
}

Make_Plot_Of_Rate_Parameter <- function(rate_samples)
{
  require("coda")
  mw <- mean(rate_samples)
  sdev <- sd(rate_samples)
  df_gamma_dist_params <- Get_Gamma_Distribution_Parameters_From_Moments(mw,sdev)
  
  ks <- seq(from=4e-4,to=12e-4,by=2e-6)
  probs <- dgamma(ks,shape=df_gamma_dist_params$k,scale=df_gamma_dist_params$theta)
  plot(ks,probs,type="l",lty=2,lwd=2,xlab="k_on [1/(nM s)]",ylab="probability density")
  hist(mat[,2],breaks=20,add=TRUE,freq=FALSE)
  
  hdi_values <- HPDinterval(as.mcmc(mat[,2]))
  lp <- hdi_values[1];
  up <- hdi_values[2];
  xvs_1 <- c(lp,lp)
  xvs_2 <- c(up,up)
  points(xvs_1,c(0,5000),type="l",lty=2,col="red")
  points(xvs_2,c(0,5000),type="l",lty=2,col="red")
  points(c(mw,mw),c(0,5000),type="l",lty=2,col="blue")
  df_cred_interv <- data.frame(mw,lp,up)
  return(df_cred_interv)
}
  

