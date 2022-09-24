# F-LIML and S-LIML tests of H0:theta=tet.0

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
# k0 = (1-upsilon) screening threshold for selecting relevant factors (default: screen at 0.01 level tests)

## Outputs
# fliml.est = F-LIML estimate of the causal association of X and Y
# fliml.se = standard error of F-LIML estimate
# fliml = p-value of the F-LIML test of the null of no causal association
# r = selected number of factors after screening at threshold k0
# sliml.est = S-LIML estimate of the causal association of X and Y
# N.B. cil and ciu are not the confidence intervals for tet0, but in fact... 
# cil: lower 2.5% tail of the distribution of sliml.est under the null H0: true theta is tet.0 
# ciu: upper 2.5% tail of the distribution of sliml.est under the null H0: true theta is tet.0
# result: if sliml.est does not lie between cil and ciu then we reject H0


SLIML <- function(dx,dy,vz,ld,tet.0,vx,vy,nx,ny,r=NULL,k0=NULL){
  p <- nrow(ld)
  varZ0 <- ld*(sqrt(vz)%*%t(sqrt(vz)))
  if(missing(r)) {
    r=which(cumsum(prcomp(varZ0,scale=FALSE)$sdev^2/sum((prcomp(varZ0,scale=FALSE)$sdev^2)))>0.99)[1]
  } else {
    r=r
  }
  if(missing(k0)) {
    k0=1-0.01
  } else {
    k0=k0
  }
  
  # estimate factor loadings
  lambda <- sqrt(p)*(eigen(varZ0)$vectors[,1:r])
  evec <- eigen((t(lambda)%*%lambda))$vectors
  eval <- eigen((t(lambda)%*%lambda))$values
  lambda <- lambda%*%(solve(evec%*%diag(sqrt(eval))%*%t(evec)))
  dim(lambda) <- c(p,r)
  
  # F-LIML inference
  sv0 <- vx-as.numeric(t((t(lambda)%*%dx))%*%solve(t(lambda)%*%varZ0%*%lambda)%*%t(lambda)%*%dx)
  se0 <- vy-as.numeric(t((t(lambda)%*%dy))%*%solve(t(lambda)%*%varZ0%*%lambda)%*%t(lambda)%*%dy)
  Om <- function(tet){t(lambda)%*%(varZ0*((se0/ny)+((tet^2)*(sv0/nx))))%*%lambda}
  G <- -t(lambda)%*%dx
  g <- function(tet){t(lambda)%*%(dy-dx*tet)}
  Q <- function(tet){as.numeric(t(g(tet))%*%solve(Om(tet))%*%g(tet))}
  init.val <- seq(-1,1,0.2)
  Q.init <- vector(,length=length(init.val))
  for(l in 1:length(init.val)){
    Q.init[l]<-optim(init.val[l], Q, method="Brent",lower=-1e2,upper=1e2)$value
  }
  tet_est <- optim(init.val[which.min(Q.init)[[1]]], Q, method="Brent",lower=-1e2,upper=1e2)$par
  tet_se <- sqrt(as.numeric(solve(t(G)%*%solve(Om(tet_est))%*%G)))
  fliml.pval <- (1-pnorm(abs(tet_est/tet_se),0,1))*2

  # S-LIML
  c0 <- rep(qnorm(1-((1-k0)/2)),r)
  delGG <- t(lambda)%*%(varZ0*(sv0/nx))%*%lambda
  D <- diag(diag(delGG))
  t0 <- (1/sqrt(diag(D)))*as.vector(G)
  sel <- (abs(t0)>c0)*1
  r.sel <- sum(sel==1)
  
  if(r.sel>0){
    lambda.sel <- lambda[,which(sel==1)]
    G.sel <- -t(lambda.sel)%*%dx
    sv0.sel <- 1-as.numeric(t(G.sel)%*%solve(t(lambda.sel)%*%varZ0%*%lambda.sel)%*%G.sel)
    Om.sel <- function(tet){t(lambda.sel)%*%(varZ0*((se0/ny)+((tet^2)*(sv0.sel/nx))))%*%lambda.sel}
    g.sel <- function(tet){t(lambda.sel)%*%(dy-dx*tet)}
    Q.sel <- function(tet){as.numeric(t(g.sel(tet))%*%solve(Om.sel(tet))%*%g.sel(tet))}
    init.val <- seq(-1,1,0.2)
    Q.init <- vector(,length=length(init.val))
    for(l in 1:length(init.val)){
      Q.init[l]<-optim(init.val[l], Q.sel, method="Brent",lower=-1e2,upper=1e2)$value
    }
    tet_sel <- optim(init.val[which.min(Q.init)[[1]]], Q.sel, method="Brent",lower=-1e2,upper=1e2)$par
    Om_sel <- Om.sel(tet_sel)
    V.sel <- 1/as.numeric(t(G.sel)%*%solve(Om_sel)%*%G.sel)
    C <- as.vector(-sqrt(solve(D))%*%(t(lambda)%*%(varZ0*(sv0/nx))%*%lambda.sel)%*%solve(Om_sel)%*%G.sel%*%V.sel*tet.0)
    u <- t0-(C*(1/V.sel)*tet_sel)
    A <- function(Z){
      v <- as.vector(u+(C*(1/V.sel)*(tet.0+(sqrt(V.sel)*Z))))
      a1 <- prod((abs(v[which(sel==1)])>c0[which(sel==1)])*1)
      a2 <- prod((abs(v[which(sel==0)])<=c0[which(sel==0)])*1)
      A.est <- a1*a2
      return(A.est)
    }
    z <- rnorm(5e3,0,1)
    A <- sapply(z,A)
    Pn <- function(Y){mean(ifelse(z <= ((1/sqrt(V.sel))*(Y-tet.0)),1,0)*A)/mean(A)}
    y <- seq(-3,3,0.001)
    P <- sapply(y,Pn)
    cil <- as.numeric(y[max(which(P<=0.025))])
    ciu <- as.numeric(y[min(which(P>=0.975))])
    if(tet_sel<=ciu & tet_sel>=cil){result <- (paste0("H0: S-LIML cannot reject a true value of ",tet.0, " for a 0.05 level test"))} else  {result <- (paste0("H1: S-LIML rejects a true value of ",tet.0, " for a 0.05 level test"))} 
    res.list <- list("fliml.est"=tet_est, "fliml.se"=tet_se, "fliml"=fliml.pval,"r"=r.sel, "sliml.est"=tet_sel, "cil"=cil, "ciu"=ciu, "result"=result)
  }
  
  if(r.sel==0){
    res.list <- list("fliml.est"=tet_est, "fliml.se"=tet_se, "fliml"=fliml.pval, "r"=r.sel, "sliml.est"=NA, "sliml.se"=NA, "sliml"=NA, "cil"=NA, "ciu"=NA, "result"=NA)
  }
  return(res.list)
}