
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
  
  ##draw initial concentration from the prior 
  A0_init <- rnorm(n=1,mean=A0_mu,sd=A0_sigma)
  init_stan_fit_dimerization1 <- list(list(sigma=40,k_on=1e-4,A_0=A0_init,amp_cy35=750,bg_cy35=5000))
  
  sf_dimerization_reaction1 <- stan(file="/Users/christianwachauf/Documents/Scripts/GithubRepo/RT_PCR_Evaluation/Stan/dimerization_reaction1.stan",data=data_stan_fit_dimerization1,
                                    init=init_stan_fit_dimerization1,iter=num_iters,chain=1)
  return(sf_dimerization_reaction1)
}

Stan_Fit_Simple_Dimerization_Reaction1_Global <- function(times,cy35_curves,A0_mus,A0_sigmas,amps_mus,bgs_mus,num_iters=200)
{
  require("rstan")
  num_traces <- length(A0_mus)
  data_stan_fit_dimerization1_gf <- list(num_times=length(times),num_traces=num_traces,ts=times,intens=cy35_curves,mu_A0s=A0_mus,sigma_A0s=A0_sigmas,
                                      mu_amps_cy35=amps_mus,mu_bgs_cy35=bgs_mus)
  
  A0s_init <- rnorm(n=num_traces,mean=A0_mus,sd=A0_sigmas)
  amps_cy35_init <- rnorm(n=num_traces,mean=amps_mus,sd=200)
  bgs_cy35_init <- rnorm(n=num_traces,mean=bgs_mus,sd=500)
  
  init_stan_fit_dimerization1_gf <- list(list(sigma=40,k_on=1e-4,A0s=A0s_init,amps_cy35=amps_cy35_init,bgs_cy35 = bgs_cy35_init))
  sf_dimerization_reaction1_gf <- stan(file="/Users/christianwachauf/Documents/Scripts/GithubRepo/RT_PCR_Evaluation/Stan/dimerization_reaction1_global.stan",
                                       data=data_stan_fit_dimerization1_gf,init=init_stan_fit_dimerization1_gf,iter=num_iters,chain=1)
  return(sf_dimerization_reaction1_gf)
}


Test_Stan_Fit_DimReac1Glob <- function(data_rtpcr,num_iters=200)
{
  times <- data_rtpcr$time_seq[1:720]
  cy35_2nM <- data_rtpcr$`24hb_2nM_40nM_dos_Cy35`[1:720]##24hb_2nM_40nM_dos_Cy35[1:720]
  cy35_4nM <- data_rtpcr$`24hb_4nM_40nM_dos_Cy35`[1:720]##24hb_4nM_40nM_dos_Cy35[1:720]
  cy35_data <- array(0,dim=c(2,length(times)))
  
  cy35_data[1,] <- cy35_2nM
  cy35_data[2,] <- cy35_4nM
  
  A0_mus <- c(2,4)
  A0_sigmas <- c(0.25,0.5)
  amps_mus <- c(800,800)
  bgs_mus <- c(5000,5000)
  sf_dimerization_reaction1_gf <- Stan_Fit_Simple_Dimerization_Reaction1_Global(times,cy35_data,A0_mus=A0_mus,
                                                                              A0_sigmas = A0_sigmas,amps_mus=amps_mus,bgs_mus=bgs_mus,num_iters=num_iters)
  return(sf_dimerization_reaction1_gf)
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
   cy35_data <- data_rtpcr$`24hb_2nM_40nM_pos_Cy35`
   plot(time_data,cy35_data,xlim=c(0,10000),xlab="time [s]",ylab="intensity [a.u.]")
   time_sequ <- seq(from=0,to=10000,by=5)
   df_plot_data <- Plot_Mean_Posterior_Simple_Dimerization_Reaction1(sf_dimerization_reaction1,time_sequ)
   points(df_plot_data$times,df_plot_data$int_cy35,type="l",col="red")
   
}

Return_Prediction_Intervals_For_Posterior_Dimerization_Reaction1_Global <- function(sf_dimerization1_gf,times)
{
  require("coda")
  mat_gf <- as.matrix(sf_dimerization1_gf)
  

  num_samples <- nrow(mat_gf)
  ##num_samples <- 100 ## temporary !!
  
  num_cols <- ncol(mat_gf)
  num_times <- length(times)
  print(num_cols)
  num_traces <- (num_cols-3)/3
  print("number of traces: ")
  print(num_traces)
  
  
  int_mean <- array(0,dim=c(num_traces,num_times))
  int_lp <- array(0,dim=c(num_traces,num_times))
  int_up <- array(0,dim=c(num_traces,num_times))
  
  A0s_curr <- array(0,dim=c(num_traces))
  amps_cy35_curr <- array(0,dim=c(num_traces))
  bgs_cy35_curr <- array(0,dim=c(num_traces))
  predictions_intensity <- array(0,dim=c(num_traces,num_samples,length(times)))
  for(i in 1:num_samples)
  {
      sigma_curr <- mat_gf[i,1]
      k_on_curr <- mat_gf[i,2]
      for(j in 1:num_traces)
      {
        A0s_curr[j] <- mat_gf[i,2+j]
        amps_cy35_curr[j] <- mat_gf[i,2+num_traces+j]
        bgs_cy35_curr[j] <- mat_gf[i,2+2*num_traces+j]
        for(t in 1:length(times))
        {
          predictions_intensity[j,i,t] <- Dimer_Concentration_Simple_Dimerization1(A0s_curr[j],k_on_curr,times[t])
          predictions_intensity[j,i,t] <- predictions_intensity[j,i,t]* amps_cy35_curr[j] + bgs_cy35_curr[j]
        ## add noise (i)
          predictions_intensity[j,i,t] <- predictions_intensity[j,i,t] + rnorm(n=1,mean=0,sd=sigma_curr)
        }
      }
  }
  for(j in 1:num_traces)
  {
    for(t in 1:length(times))
    {
      ##t_mean[j,t] <- mean(predictions_intensity[,j])
      int_mean[j,t] <- mean(predictions_intensity[j,,t])
      hdi_values <- HPDinterval(as.mcmc(predictions_intensity[j,,t]))
      int_lp[j,t] <- hdi_values[1];
      int_up[j,t] <- hdi_values[2];
    }
    ##df_preds <- data.frame(times,int_mean,int_lp,int_up)
    ##return(df_preds)
  }
  df_preds <- data.frame(int_mean[1,],int_lp[1,],int_up[1,],int_mean[2,],int_lp[2,],int_up[2,])
  return(df_preds)
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

Plot_Data_Dimerization_Reaction1_With_Preds_Global <- function(times,cy35_data,df_preds)
{
  max_val <- max(cy35_data)
  min_val <- min(cy35_data)
  
  d_intens <- (max_val - min_val)
  
  max_val <- max_val + 0.05*d_intens;
  min_val <- min_val - 0.05*d_intens;
  
  
  plot(times,cy35_data[1,],xlab="time [s]",ylab="intensity [a.u.]",ylim=c(min_val,max_val))
  points(times,cy35_data[2,])
  
  points(times,df_preds[,1],type="l",lwd=2,col="red")
  points(times,df_preds[,4],type="l",lwd=2,col="red")
  
  points(times,df_preds[,2],type="l",lwd=1,col="blue")
  points(times,df_preds[,3],type="l",lwd=1,col="blue")
  points(times,df_preds[,5],type="l",lwd=1,col="blue")
  points(times,df_preds[,6],type="l",lwd=1,col="blue")
}
  
Plot_Data_Dimerization_Reaction1_With_Preds <- function(times,cy35_data,df_preds)
{
  plot(times,cy35_data,xlab="time [s]",ylab="intensity [a.u.]")
  ## add prediction bands
  points(df_preds$times,df_preds$int_mean,type="l",lwd=2,col="red")
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
  
  ks <- seq(from=8e-4,to=12e-4,by=2e-6)
  probs <- dgamma(ks,shape=df_gamma_dist_params$k,scale=df_gamma_dist_params$theta)
  plot(ks,probs,type="l",lty=2,lwd=2,xlab="k_on [1/(nM s)]",ylab="probability density")
  hist(rate_samples,breaks=10,add=TRUE,freq=FALSE)
  
  hdi_values <- HPDinterval(as.mcmc(rate_samples))
  lp <- hdi_values[1];
  up <- hdi_values[2];
  xvs_1 <- c(lp,lp)
  xvs_2 <- c(up,up)
  points(xvs_1,c(0,90000),type="l",lty=2,col="red")
  points(xvs_2,c(0,90000),type="l",lty=2,col="red")
  points(c(mw,mw),c(0,90000),type="l",lty=2,col="blue")
  df_cred_interv <- data.frame(mw,lp,up)
  return(df_cred_interv)
}

Test_Plot_Model2 <- function(data_rtpcr,num_iter=200)
{
  require("rstan")
  ts <- data_rtpcr$time_seq[seq(from=1,to=720,by=3)];
  num_traces <- 2;
  intens <- array(0,dim=c(num_traces,length(ts)))
  intens[1,] <- data_rtpcr$`24hb_4nM_12nM_dos_Cy35`[seq(from=1,to=720,by=3)]
  intens[2,] <- data_rtpcr$`24hb_4nM_40nM_dos_Cy35`[seq(from=1,to=720,by=3)]
  
  mu_A0s <- c(4,4);
  mu_B0s <- c(4,4);
  
  sigma_A0s <- c(0.4,0.4);
  sigma_B0s <- c(0.4,0.4);
  
  mu_D1_0s <- c(12,40);
  sigma_D1_0s <- c(1.2,4.0);
  
  mu_D2_0s <- c(12,40);
  sigma_D2_0s <- c(1.2,4.0);
  
  init_A0s <- rnorm(n=num_traces,mean=mu_A0s,sd=sigma_A0s)
  init_B0s <- rnorm(n=num_traces,mean=mu_B0s,sd=sigma_B0s)
  init_D1_0s <- rnorm(n=num_traces,mean=mu_D1_0s,sd=sigma_D1_0s)
  init_D2_0s <- rnorm(n=num_traces,mean=mu_D2_0s,sd=sigma_D2_0s)
  
  mu_amps_cy35 <- c(800,800)
  mu_bgs_cy35 <- c(5000,5000)
  
  data_stan_dim_reaction2 <- list(num_times = length(ts),num_traces=num_traces,ts=ts,
                                intens=intens,mu_A0s=mu_A0s,sigma_A0s=sigma_A0s,mu_B0s=mu_B0s,sigma_B0s=sigma_B0s,
                                mu_D1_0s=mu_D1_0s,sigma_D1_0s=sigma_D1_0s,mu_D2_0s=mu_D2_0s,sigma_D2_0s=sigma_D2_0s,mu_amps_cy35=mu_amps_cy35,
                                mu_bgs_cy35=mu_bgs_cy35)
  ## init: sigma, all concentrations, k_on1,k_on2
  init_stan_dim_reaction2 <- list(list(sigma=40,k_on1=1e-4,k_on2=5e-4,A0s=init_A0s,B0s=init_B0s,D01s=init_D1_0s,D02s=init_D2_0s,amps_cy35=c(800,800),bgs_cy35=c(5000,5000)))
  
  sf_dim_reaction2 <- stan(file="/Users/christianwachauf/Documents/Scripts/GithubRepo/RT_PCR_Evaluation/Stan/dimerization_reaction2.stan",data=data_stan_dim_reaction2,
                         init=init_stan_dim_reaction2,iter=num_iter,chain=1)
  return(sf_dim_reaction2)
}

Test_Some_More_Shit <- function(data_rtpcr,num_iters=100)
{
  require("rstan")
  ts <- data_rtpcr$time_seq[seq(from=1,to=720,by=3)];
  cy35_1 <- data_rtpcr$`24hb_4nM_12nM_dos_Cy35`[seq(from=1,to=720,by=3)]
  cy35_2 <- data_rtpcr$`24hb_4nM_40nM_dos_Cy35`[seq(from=1,to=720,by=3)]
 # print(ts)
#  print(cy35_1)
  a0_1 <- 4;
  d0_1 <- 12;
  a0_2 <- 4;
  d0_2 <- 40;
  
  t0 <- 0;
  
  sf_data_dimerization3 <- list(T=length(ts),cy35_1 =cy35_1,cy35_2=cy35_2,ts=ts,a0_1=a0_1,a0_2=a0_2,d0_1=d0_1,d0_2=d0_2,t0=0)
  
  sf_init_dimerization3 <- list(list(sigma=40,k1=1e-4,k2=5e-4,amp_cy35_1=800,bg_cy35_1=5000,amp_cy35_2=800,bg_cy35_1=5000))
  sf_dim_reaction3 <- stan(file="/Users/christianwachauf/Documents/Scripts/GithubRepo/RT_PCR_Evaluation/Stan/dimerization_reaction3.stan",data=sf_data_dimerization3,
                           init=sf_init_dimerization3,iter=num_iters,chain=1)
  return(sf_dim_reaction3)
}

Test_Some_More_Shit3 <- function(sf_shit)
{

  mat_shit <- as.matrix(sf_shit)
  sigma_mean <- mean(mat_shit[,1])
  k1_mean <- mean(mat_shit[,2])
  k2_mean <- mean(mat_shit[,3])
  amp_cy35_1_mean <- mean(mat_shit[,4])
  bg_cy35_1_mean <- mean(mat_shit[,5])
  amp_cy35_2_mean <- mean(mat_shit[,6])
  bg_cy35_2_mean <- mean(mat_shit[,7])
  
  ts <- seq(from=30,to=7200,by=30)
  stan_data_fit_dim_react3 <- list(T=length(ts),y0_1=c(4,12,0,0),y0_2=c(4,40,0,0),t0=0,ts=ts,theta=c(k1_mean,k2_mean),amp_cy35_1 = amp_cy35_1_mean,
                                   bg_cy35_1 = bg_cy35_1_mean,amp_cy35_2 = amp_cy35_2_mean,bg_cy35_2=bg_cy35_2_mean,sigma=sigma_mean)
  fit_stan_dim_react3 <- stan(file="/Users/christianwachauf/Documents/Scripts/GithubRepo/RT_PCR_Evaluation/Stan/fit_dimerization_reaction3.stan",data=stan_data_fit_dim_react3,
                              iter=1,chain=1,algorithm="Fixed_param")
  return(fit_stan_dim_react3)
}


Plot_Shit_Dim_React3 <- function(data_rtpcr,sf_fit_shit)
{
  array_res <- as.array(sf_fit_shit)
  plot(data_rtpcr$time_seq[1:720],data_rtpcr$`24hb_4nM_12nM_dos_Cy35`[1:720],xlab="time [s]",ylab="intensity [a.u.]",ylim=c(4800,9000))
  points(data_rtpcr$time_seq[1:720],data_rtpcr$`24hb_4nM_40nM_dos_Cy35`[1:720])
  points(ts,array_res[1:240],type="l",col="red",lwd=2)
  points(ts,array_res[241:480],type="l",col="red",lwd=2)
  
}


Stan_Fit_Dimerization_Reaction4 <- function(data_rt_pcr,num_iters)
{
  num_traces <- 4;
  indices <- seq(from=1,to=720,by=3)
  num_times <- length(indices)
  times <- data_rtpcr$time_seq[indices]
  cy35_data <- array(0,dim=c(num_traces,num_times))
  cy35_data[1,] <- data_rtpcr$`24hb_2nM_12nM_dos_Cy35`[indices];
  cy35_data[2,] <- data_rtpcr$`24hb_2nM_40nM_dos_Cy35`[indices];
  cy35_data[3,] <- data_rtpcr$`24hb_4nM_12nM_dos_Cy35`[indices];
  cy35_data[4,] <- data_rtpcr$`24hb_4nM_40nM_dos_Cy35`[indices];
  
  A0s <- c(2,2,4,4);
  D0s <- c(12,40,12,40);
  
  bgs_cy35_init <- c(5000,5000,5000,5000)
  amps_cy35_init <- c(800,800,800,800)
  
  data_sf_dim_react4 <- list(num_times=num_times,num_traces=num_traces,ts=times,cy35_data=cy35_data,A0s=A0s,D0s=D0s,t0=0)
  init_sf_dim_react4 <- list(list(sigma=40,k1=1e-4,k2=5e-4,bgs_cy35=bgs_cy35_init,amps_cy35=amps_cy35_init))
  sf_dim_react4 <- stan(file="/Users/christianwachauf/Documents/Scripts/GithubRepo/RT_PCR_Evaluation/Stan/dimerization_reaction4.stan",
                        data=data_sf_dim_react4,init=init_sf_dim_react4,chain=1,iter=num_iters)
  return(sf_dim_react4)
}

Plot_Dimerization_Reaction4 <- function(data_rtpcr,sf_dim_react4)
{
 mat_dr4 <- as.matrix(sf_dim_react4)
 num_traces <- 4
 ts <- seq(from=30,to=7200,by=30)
 print(head(mat_dr4))
 k1_mean <- mean(mat_dr4[,2])
 k2_mean <- mean(mat_dr4[,3])
 bgs_cy35_means <- array(0,dim=c(num_traces))
 amps_cy35_means <- array(0,dim=c(num_traces))
 for(i in 1:num_traces)
 {
   bgs_cy35_means[i] <- mean(mat_dr4[,3+i])
   amps_cy35_means[i] <- mean(mat_dr4[,3+num_traces+i])
 }
 print(bgs_cy35_means)
 print(amps_cy35_means)
 #plot(bgs_cy35_means,amps_cy35_means)
 init_concs <- array(0,dim=c(num_traces,4))
 
 init_concs[1,] <- c(2,12,0,0)
 init_concs[2,] <- c(2,40,0,0)
 init_concs[3,] <- c(4,12,0,0)
 init_concs[4,] <- c(4,40,0,0)
 
 data_stan_sim_dimreact4 <- list(num_times=length(ts),num_traces=4,ts=ts,k1=k1_mean,k2=k2_mean,
                                  bgs_cy35 = bgs_cy35_means,amps_cy35=amps_cy35_means,init_concs=init_concs,t0=0)
 
 stan_sim_dimreact4 <- stan(file="/Users/christianwachauf/Documents/Scripts/GithubRepo/RT_PCR_Evaluation/Stan/sim_dimerization_reaction4.stan",
                            data=data_stan_sim_dimreact4,iter=1,chain=1,algorithm="Fixed_param") 
 ar_dr4 <- as.array(stan_sim_dimreact4)
 
 ## start plotting the data:
 data_indices <- seq(from=1,to=720,by=3)
 
 plot(data_rtpcr$time_seq[data_indices],data_rtpcr$`24hb_4nM_12nM_dos_Cy35`[data_indices],ylim=c(4800,9000),xlab="time [s]",ylab="intensity [a.u.]")
 points(data_rtpcr$time_seq[data_indices],data_rtpcr$`24hb_4nM_40nM_dos_Cy35`[data_indices])
 points(data_rtpcr$time_seq[data_indices],data_rtpcr$`24hb_2nM_12nM_dos_Cy35`[data_indices])
 points(data_rtpcr$time_seq[data_indices],data_rtpcr$`24hb_2nM_40nM_dos_Cy35`[data_indices])
 points(ts,ar_dr4[seq(from=1,to=960,by=4)],type="l",col="red",lwd=2)
 points(ts,ar_dr4[seq(from=2,to=960,by=4)],type="l",col="red",lwd=2)
 points(ts,ar_dr4[seq(from=3,to=960,by=4)],type="l",col="red",lwd=2)
 points(ts,ar_dr4[seq(from=4,to=960,by=4)],type="l",col="red",lwd=2)
 
 plot(data_rtpcr$time_seq[data_indices],data_rtpcr$`24hb_2nM_12nM_dos_Cy35`[data_indices])
 points(ts,ar_dr4[seq(from=1,to=960,by=4)],type="l",col="red",lwd=2)
 
 plot(data_rtpcr$time_seq[data_indices],data_rtpcr$`24hb_2nM_40nM_dos_Cy35`[data_indices])
 points(ts,ar_dr4[seq(from=2,to=960,by=4)],type="l",col="red",lwd=2)
 
 plot(data_rtpcr$time_seq[data_indices],data_rtpcr$`24hb_4nM_12nM_dos_Cy35`[data_indices])
 points(ts,ar_dr4[seq(from=3,to=960,by=4)],type="l",col="red",lwd=2)
 
 plot(data_rtpcr$time_seq[data_indices],data_rtpcr$`24hb_4nM_40nM_dos_Cy35`[data_indices])
 points(ts,ar_dr4[seq(from=4,to=960,by=4)],type="l",col="red",lwd=2)
 
 return(stan_sim_dimreact4)
}

## use Stan to simulate some dimerization reactions
## that are following the kinetic scheme:
## A + B --> C   (k_on1)
## C + C <--> D  (k_on2, k_off2)
Stan_Sim_Dimerization_Reaction5 <- function()
{
  require("rstan")
  times <- seq(from=10,to=30000,by=10)
  num_times <- length(times)
  num_traces <- 2
  k_on1 <- 2e-4  ## 2e-4 1/(nM s)
  k_on2 <- 10e-4 ## 10e-4 1/(nM s)
  k_off2 <- 1e-5 ## 1e-5 1/s
  amps_cy35 <- c(700,700)
  bgs_cy35 <- c(5000,5000)
  
  init_concs <- array(0,dim=c(2,4))
  init_concs[1,] <- c(4,40,0,0)
  init_concs[2,] <- c(2,40,0,0)
  t0<-0
  
  data_sim_stan_dimerization5 <- list(num_times=num_times,num_traces=num_traces,ts=times,k_on1=k_on1,k_on2=k_on2,k_off2=k_off2,amps_cy35=amps_cy35,
                                      bgs_cy35=bgs_cy35,init_concs=init_concs,t0=t0)
  
  stan_sim_dimerization5 <- stan(file="/Users/christianwachauf/Documents/Scripts/GithubRepo/RT_PCR_Evaluation/Stan/sim_dimerization_reaction5.stan",data = data_sim_stan_dimerization5,
                                 iter=1,chain=1,algorithm="Fixed_param")
  return(stan_sim_dimerization5)
}

## Stan_Fit_Dimerization_Reaction5(cy35_data,time)
Stan_Fit_Dimerization_Reaction5 <- function(data_rtpcr,num_iters=100)
{
  ##num_traces <- nrow(cy35_data)
  ##print("number of different traces: ")
  ##print(num_traces)
  num_traces <- 1
  indices <- seq(from=1,to=720,by=3)
  num_times <- length(indices)
  times <- data_rtpcr$time_seq[indices]
  cy35_data <- array(0,dim=c(1,num_times))
  cy35_data[1,] <- data_rtpcr$`24hb_6dim_oligos_Cy35`[indices]
  
  init_concs <- array(0,dim=c(1,4))
  init_concs[1,] <- c(4,40,0,0);
  t0 <- 0;
  data_sf_dim_react5 <- list(num_times=num_times,num_traces=num_traces,ts=times,cy35_data=cy35_data,init_concs=init_concs,t0=t0)
  
  bgs_cy35_init <- array(0,dim=c(1))
  amps_cy35_init <- array(0,dim=c(1))
  bgs_cy35_init[1] <- 5000
  amps_cy35_init[1] <- 700
  init_sf_dim_react5 <- list(list(sigma=40,k_on1=2e-4,k_on2=10e-4,k_off2=1e-5,bgs_cy35=bgs_cy35_init,amps_cy35=amps_cy35_init))
  sf_dim_react5 <- stan(file="/Users/christianwachauf/Documents/Scripts/GithubRepo/RT_PCR_Evaluation/Stan/fit_dimerization_reaction5.stan",data=data_sf_dim_react5,
                        init=init_sf_dim_react5,iter=num_iters,chain=1)
  return(sf_dim_react5)
}

#parameters
#{
#  real<lower=0,upper=1000> sigma; // noise, supposed to be equal in both cases...
#  real<lower=0,upper=0.1> k_on1; // rate constants (theta[1] and theta[2])
#  real<lower=0,upper=0.1> k_on2;
#  real<lower=0,upper=0.1> k_off2;
  
#  real<lower=0,upper=10000> bgs_cy35[num_traces];
#  real<lower=0,upper=3000> amps_cy35[num_traces];
#}#

#data
#{
#  int<lower=1> num_times; // number of data points
#  int<lower=1> num_traces;
#  
#  real ts[num_times]; // times at which model data is given#
#  real cy35_data[num_traces,num_times];
#  
#  real<lower=0> init_concs[num_traces,4];
#  real t0;
#}

## Get_Equilibrium_Concentrations(c_A0,c_B0,K_D)
## returns for a bimolecular reaction with initial concentrations
## "c_A0","c_B0" und "K_D"
Get_Equilibrium_Concentrations <- function(c_A0,c_B0,K_D)
{
  term1 <- K_D + c_B0 - c_A0
  c_A_eq <- -0.5*(term1) + sqrt(0.25*term1^2+K_D*c_A0)
  c_B_eq <- c_A_eq + c_B0 - c_A0
  c_AB_eq <- c_A_eq*c_B_eq/K_D
  
  df_eq_concs <- data.frame(c_A_eq,c_B_eq,c_AB_eq)
  return(df_eq_concs)
}
