## evaluation protocol for the "15.09.2016" - data of dimerization experiments
## experiments with 4 dimerization oligos vs 24 dimerization oligos..

source("fit_rtpcr_models.R")
source("parsertdata.r")
Load_Raw_Data <- function()
{
  filename <- "measurement_24hb_dimerization_4dim_oligos_vs_24dim_oligos_15_09_2016.xls"
  filename_mp <-"rtpcr_protocols/rtpcr_protocol1.xls"
  list_ommp <- Get_List_Ommp_From_File(filename_mp)
  num_wells <- 5
  well_names <- list("24hb_2nM_40nM_4dos","24hb_4nM_40nM_4dos","24hb_6nM_40nM_4dos","24hb_4nM_40nM_24dos","buffer")
  mp <- InitMeasurementProtocol(list_ommp,num_wells,well_names)
  rtpcr_data <- LoadRTPCRData(filename,mp)
}

Evaluate_Data <- function(rtpcr_data,num_iters=50)
{
  require("rstan")
  indices <- seq(from=1,to=1440,by=3)
  plot(rtpcr_data$time_seq[indices],rtpcr_data$`24hb_2nM_40nM_4dos_Cy35`[indices],ylim=c(6000,12200))
  points(rtpcr_data$time_seq[indices],rtpcr_data$`24hb_4nM_40nM_4dos_Cy35`[indices])
  points(rtpcr_data$time_seq[indices],rtpcr_data$`24hb_6nM_40nM_4dos_Cy35`[indices])
  
  num_traces <- 3
  num_times <- length(indices)
  times <- rtpcr_data$time_seq[indices]
  cy35_data <- array(0,dim=c(num_traces,num_times))
  
  cy35_data[1,] <- rtpcr_data$`24hb_2nM_40nM_4dos_Cy35`[indices]
  cy35_data[2,] <- rtpcr_data$`24hb_4nM_40nM_4dos_Cy35`[indices]
  cy35_data[3,] <- rtpcr_data$`24hb_6nM_40nM_4dos_Cy35`[indices]
  
  init_concs <- array(0,dim=c(num_traces,4))
  init_concs[1,] <- c(2,40,0,0);
  init_concs[2,] <- c(4,40,0,0);
  init_concs[3,] <- c(6,40,0,0);
  
  t0 <- 0
  
  data_sf_dim_react5 <- list(num_times=num_times,num_traces=num_traces,ts=times,cy35_data=cy35_data,init_concs=init_concs,t0=t0)
  
  bgs_cy35_init <- c(6000,6000,6000)
  amps_cy35_init <- c(800,800,800)
  init_sf_dim_react5 <- list(list(sigma=40,k_on1=2e-4,k_on2=10e-4,k_off2=1e-5,bgs_cy35=bgs_cy35_init,amps_cy35=amps_cy35_init))
  sf_dim_react5 <- stan(file="/Users/christianwachauf/Documents/Scripts/GithubRepo/RT_PCR_Evaluation/Stan/fit_dimerization_reaction5.stan",data=data_sf_dim_react5,
                        init=init_sf_dim_react5,iter=num_iters,chain=1)
  return(sf_dim_react5)
}

Evaluate_Data_Dim_React6 <- function(rtpcr_data,num_iters=50)
{
  require("rstan")
  indices <- seq(from=1,to=720,by=3)
  plot(rtpcr_data$time_seq[indices],rtpcr_data$`24hb_2nM_40nM_4dos_Cy35`[indices],ylim=c(6000,12200))
  points(rtpcr_data$time_seq[indices],rtpcr_data$`24hb_4nM_40nM_4dos_Cy35`[indices])
  points(rtpcr_data$time_seq[indices],rtpcr_data$`24hb_6nM_40nM_4dos_Cy35`[indices])
  
  num_traces <- 3
  num_times <- length(indices)
  times <- rtpcr_data$time_seq[indices]
  cy35_data <- array(0,dim=c(num_traces,num_times))
  
  cy35_data[1,] <- rtpcr_data$`24hb_2nM_40nM_4dos_Cy35`[indices]
  cy35_data[2,] <- rtpcr_data$`24hb_4nM_40nM_4dos_Cy35`[indices]
  cy35_data[3,] <- rtpcr_data$`24hb_6nM_40nM_4dos_Cy35`[indices]
  
  init_concs <- array(0,dim=c(num_traces,2))
 ## init_concs[1,] <- c(2,40);
  
  init_concs[1,] <- c(2,0);
  init_concs[2,] <- c(4,0);
  init_concs[3,] <- c(6,0);
  
  t0 <- 0
  
  data_sf_dim_react6 <- list(num_times=num_times,num_traces=num_traces,ts=times,cy35_data=cy35_data,init_concs=init_concs,t0=t0)
  
  bgs_cy35_init <- c(6000,7200,8000)
  #amps_cy35_init <- c(800,800,800)
  amps_cy35_init <- 1000
  init_sf_dim_react6 <- list(list(sigma=40,k_on=10e-4,k_off=1e-4,bgs_cy35=bgs_cy35_init,amps_cy35=amps_cy35_init))
  sf_dim_react6 <- stan(file="/Users/christianwachauf/Documents/Scripts/GithubRepo/RT_PCR_Evaluation/Stan/fit_dimerization_reaction6.stan",data=data_sf_dim_react6,
                        init=init_sf_dim_react6,iter=num_iters,chain=1)
  return(sf_dim_react6)
}

Simulate_Dim_React6 <- function()
{
  times <- seq(from=30,to=7200,by=30)
  num_times <- length(times)
  num_traces <- 2
  init_concs <- array(0,dim=c(num_traces,2))
  init_concs[1,] <- c(4,0)
  init_concs[2,] <- c(6,0)
  
  bgs_cy35 <- c(7200,7800)
  amps_cy35 <- c(900,900)
  t0 <- 0
  k_on <- 4e-4;
  k_off <- 2e-5;
  data_stan_sim_dim_react6 <- list(num_times=num_times,num_traces=num_traces,bgs_cy35=bgs_cy35,amps_cy35=amps_cy35,k_on=k_on,k_off=k_off,ts=times,init_concs=init_concs,t0=t0)
  stan_sim_dim_react6 <- stan(file="/Users/christianwachauf/Documents/Scripts/GithubRepo/RT_PCR_Evaluation/Stan/sim_dimerization_reaction6.stan",data=data_stan_sim_dim_react6,
                              iter=1,chain=1,algorithm="Fixed_param")
  return(stan_sim_dim_react6)

}

Evaluate_Test_Data <- function(df_data,num_iters=50)
{
  require("rstan")
 ## indices <- seq(from=1,to=720,by=3)
  ##plot(rtpcr_data$time_seq[indices],rtpcr_data$`24hb_2nM_40nM_4dos_Cy35`[indices],ylim=c(6000,12200))
##  points(rtpcr_data$time_seq[indices],rtpcr_data$`24hb_4nM_40nM_4dos_Cy35`[indices])
 ## points(rtpcr_data$time_seq[indices],rtpcr_data$`24hb_6nM_40nM_4dos_Cy35`[indices])
  
  num_traces <- 2
  
  times <- seq(from=30,to=7200,by=30)
  num_times <- length(times)
  cy35_data <- array(0,dim=c(num_traces,num_times))
  
  cy35_data[1,] <- df_data$data1
  cy35_data[2,] <- df_data$data2
  ## cy35_data[3,] <- rtpcr_data$`24hb_6nM_40nM_4dos_Cy35`[indices]
  
  init_concs <- array(0,dim=c(num_traces,2))
  ## init_concs[1,] <- c(2,40);
  init_concs[1,] <- c(4,0);
  init_concs[2,] <- c(6,0);
  
  t0 <- 0
  
  data_sf_dim_react6 <- list(num_times=num_times,num_traces=num_traces,ts=times,cy35_data=cy35_data,init_concs=init_concs,t0=t0)
  
  bgs_cy35_init <- c(7200,7800)
  amps_cy35_init <- 920
  
  init_sf_dim_react6 <- list(list(sigma=40,k_on=5e-4,k_off=3e-5,bgs_cy35=bgs_cy35_init,amps_cy35=amps_cy35_init))
  sf_dim_react6 <- stan(file="/Users/christianwachauf/Documents/Scripts/GithubRepo/RT_PCR_Evaluation/Stan/fit_dimerization_reaction6.stan",data=data_sf_dim_react6,
                        init=init_sf_dim_react6,iter=num_iters,chain=1)
  return(sf_dim_react6)
}

Plot_From_Mean_Posterior <- function(sf_dim_react6)
{
  mat <- as.matrix(sf_dim_react6)
  k_on_mean <- mean(mat[,2])
  k_off_mean <- mean(mat[,3])
  
  
}

