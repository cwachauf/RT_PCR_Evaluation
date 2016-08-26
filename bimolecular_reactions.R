## Fitting and Plotting functions for bimolecular reactions:
Calculate_Bimolecular_conc_A_eq <- function(A0,B0,K_D)
{
  cA_T <- A0
  cB_T <- B0
  term1 <- K_D + cB_T - cA_T
  conc_A_eq <- -0.5*(term1)+sqrt(0.25*term1^2+K_D*cA_T)
  return(conc_A_eq)
}

Calculate_Bimolecular_conc_AB_eq <- function(A0,B0,K_D)
{
  AB_eq <- A0 - Calculate_Bimolecular_conc_A_eq(A0,B0,K_D)
  return(AB_eq)
}

Test_Plot_Conc_A <- function()
{
  KD <- 0.1
  B0 <- 2
  conc_A0 <- seq(from=0,to=20,by=0.05)
  conc_AB_eq <- array(0,dim=c(length(conc_A0)))
  for(i in 1:length(conc_A0))
  {
    conc_AB_eq[i] <- Calculate_Bimolecular_conc_AB_eq(conc_A0[i],B0,KD)
  }
  plot(conc_A0,conc_AB_eq)
}