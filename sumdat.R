# Recover genetic variant--trait covariances and variance of genotypes from univariable linear regression summary data

## Inputs
# beta = p-vector of estimated variant--trait effect coefficients from univariable linear regressions
# se = p-vector of standard errors of estimated variant--trait effect coefficients from univariable linear regressions
# n = (scalar) sample size of genetic association study
# trait.var = (scalar) sample variance of trait

## Outputs
# del = p-vector of estimated covariances between variants and trait
# vz = p-vector of variance of genotypes

sumdat <- function(beta,se,n,trait.var){
  vz <- trait.var/((n*(se^2))+(beta^2))
  del <- vz*beta
  res.list <- list("delta"=del, "vz"=vz)
  return(res.list)
}