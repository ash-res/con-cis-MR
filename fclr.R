# F-CLR test of H0:theta=tet.0

## Inputs
# dx = p-vector of estimated delta_X (genetic variant--exposure covariances)
# dy = p-vector of estimated delta_Y (genetic variant--outcome covariances)
# vz = p-vector of variance estimates for genetic variant
# ld = p x p matrix of genetic correlation estimates (LD matrix)
# tet.0 = tested value of the causal association (H0:theta=tet.0)
# vx = variance of exposure
# vy = variance of outcome
# nx = sample size of variant--exposure association study
# ny = sample size of variant--outcome association study
# r = number of factors (default: factors explain 99 percent of genetic variation)

## Outputs
# far = p-value of the F-AR test of the null of no causal association of X and Y
# flm = p-value of the F-LM test of the null of no causal association of X and Y
# fclr = p-value of the F-CLR test of the null of no causal association of X and Y

FCLR <- function(dx,dy,vz,ld,tet.0,vx,vy,nx,ny,r=NULL){
  p <- nrow(ld)
  varZ0 <- ld*(sqrt(vz)%*%t(sqrt(vz)))
  if(missing(r)) {
    r=which(cumsum(prcomp(varZ0,scale=FALSE)$sdev^2/sum((prcomp(varZ0,scale=FALSE)$sdev^2)))>0.99)[1]
  } else {
    r=r
  }
  
  # estimate factor loadings
  lambda <- sqrt(p)*(eigen(varZ0)$vectors[,1:r])
  evec <- eigen((t(lambda)%*%lambda))$vectors
  eval <- eigen((t(lambda)%*%lambda))$values
  lambda <- lambda%*%(solve(evec%*%diag(sqrt(eval))%*%t(evec)))
  dim(lambda) <- c(p,r)
  
  # estimate relevant derivative and variance quantities
  sv0 <- vx-as.numeric(t((t(lambda)%*%dx))%*%solve(t(lambda)%*%varZ0%*%lambda)%*%t(lambda)%*%dx)
  se0 <- vy-as.numeric(t((t(lambda)%*%dy))%*%solve(t(lambda)%*%varZ0%*%lambda)%*%t(lambda)%*%dy)
  Om <- function(tet){t(lambda)%*%(varZ0*((se0/ny)+((tet^2)*(sv0/nx))))%*%lambda}
  G <- -t(lambda)%*%dx
  g <- function(tet){t(lambda)%*%(dy-dx*tet)}
  
  # ID-robust methods
  delG <- (t(lambda)%*%(varZ0*(sv0/nx))%*%lambda)*tet.0
  delGG <- t(lambda)%*%(varZ0*(sv0/nx))%*%lambda
  g0 <- g(tet.0)
  Om <-  t(lambda)%*%(varZ0*((se0/ny)+((tet.0^2)*(sv0/nx))))%*%lambda
  Om_inv <- solve(Om)
  tilG <- G-(t(delG)%*%Om_inv%*%g0)
  var.tilG <- delGG-(t(delG)%*%Om_inv%*%delG)
  
  evec <- eigen(var.tilG)$vectors
  eval <- eigen(var.tilG)$values
  inv.se.tilG <- solve(evec%*%diag(sqrt(eval))%*%t(evec))
  evec <- eigen(Om_inv)$vectors
  eval <- eigen(Om_inv)$values
  sqrt.Om_inv0 <- evec%*%diag(sqrt(eval))%*%t(evec)
  
  S <- sqrt.Om_inv0%*%g0
  T <- inv.se.tilG%*%tilG
  Qs <- as.numeric(t(S)%*%S)
  Qt <- as.numeric(t(T)%*%T)
  Qst <- as.numeric(t(S)%*%T)
  
  AR <- Qs
  LM <- (Qst^2)/Qt
  CLR <- (Qs-Qt+sqrt(((Qs+Qt)^2)-4*(Qs*Qt-(Qst)^2)))/2
  K4 <- gamma(r/2)/(sqrt(pi)*gamma((r-1)/2))
  pval <- function(m,qt)
  {integrand <- function(m,qt,s2){pchisq(((qt+m)/(1+(qt*(s2^2)/m))),r)*((1-s2^2)^((r-3)/2))}
  pval0 <- 1-2*K4*integrate(integrand,lower=0,upper=1,m=m,qt=qt,subdivisions=5000,stop.on.error=FALSE)$value
  return(pval0)}
  
  AR.pval <- 1-pchisq(AR,r)
  LM.pval <- 1-pchisq(LM,1)
  CLR.pval <- pval(CLR,Qt)
    
  res.list <- list("far"=AR.pval, "flm"=LM.pval, "fclr"=CLR.pval)
  return(res.list)
}