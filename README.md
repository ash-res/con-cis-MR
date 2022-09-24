"Conditional inference in cis-Mendelian randomization using weak genetic factors"  
by Ashish Patel, Dipender Gill, Paul Newcombe, and Stephen Burgess

load R code to use the F-LIML and S-LIML methods: source(sliml.R)  
load R code to use the F-AR, F-LM, and F-CLR methods: source(fclr.R) 

Required data: two-sample summary data and an LD (genetic variant correlation) matrix

dx = p-vector of estimated delta_X (genetic variant--exposure covariances)  
dy = p-vector of estimated delta_Y (genetic variant--outcome covariances)  
vz = p-vector of variance estimates for genetic variants  
ld = p x p genetic variant correlation matrix (LD matrix)  
vx = variance of exposure  
vy = variance of outcome  
nx = sample size of variant--exposure association study  
ny = sample size of variant--outcome association study
  
(the LD matrix can be from based on the same sample as one of the genetic association studies, or can be from a separate reference panel).

To recover dx (likewise dy) and vz from univariable linear regression summary data, use sumdat function

Inputs:  
beta = p-vector of estimated variant--exposure effect coefficients from univariable linear regressions  
se = p-vector of standard errors of estimated variant--exposure effect coefficients from univariable linear regressions  
n = sample size of variant--exposure association study  
trait.var = sample variance of exposure  

Outputs:   
sumdat(beta,se,n,trait.var)$del = dx (p-vector of estimated delta_X; genetic variant--exposure covariances)  
sumdat(beta,se,n,trait.var)$vz = vz (p-vector of variance estimates for genetic variants)
