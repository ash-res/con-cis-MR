# F-LIML, F-CLR, F-LM, F-AR tests

# INPUT
# bx: p-vector of genetic variant-exposure associations 
# by: p-vector of genetic variant-outcome associations
# sx: p-vector of standard errors of genetic variant-exposure associations 
# sy: p-vector of standard errors of genetic variant-exposure associations
# ld: (p x p) genetic variant correlation matrix
# tet0: value to test the null hypothesis, H0: true theta = tet0
# r (optional): the number of factors. Default: chooses the number of factors that explain 99% of total variation in the p genetic variants

# OUTPUT
# tet_est: value for the F-LIML estimator which uses all r factors as instruments
# tet_se: standard error for the F-LIML estimate
# J: value of Hansen (1982, ecta.) overidentifying restrictions test statistic (or Cochran's Q statistic)
# AR.p: p-value for the factor-based AR test of H0: true theta = tet0
# LM.p: p-value for the factor-based LM test of H0: true theta = tet0
# CLR.p: p-value for the factor-based CLR test of H0: true theta = tet0

F.CLR <- function(bx,by,sx,sy,ld,tet0,r=NULL){
  p <- nrow(ld)
  if(missing(r)) {
    r=which(cumsum(prcomp(ld,scale=FALSE)$sdev^2/sum((prcomp(ld,scale=FALSE)$sdev^2)))>0.99)[1]
  } else {
    r=r
  }
lambda <- sqrt(p)*(eigen(ld)$vectors[,1:r])
evec <- eigen((t(lambda)%*%lambda))$vectors
eval <- eigen((t(lambda)%*%lambda))$values
lambda <- lambda%*%(solve(evec%*%diag(sqrt(eval))%*%t(evec)))
dim(lambda) <- c(p,r)

SigX <- diag(sx)%*%ld%*%diag(sx)
SigY <- diag(sy)%*%ld%*%diag(sy)

Om <- function(tet){t(lambda)%*%(SigY + (tet^2)*SigX)%*%lambda}
G <- -t(lambda)%*%bx
g <- function(tet){t(lambda)%*%(by-bx*tet)}
Q <- function(tet){as.numeric(t(g(tet))%*%solve(Om(tet))%*%g(tet))}
init.val <- seq(-1,1,0.2)
Q.init <- vector(,length=length(init.val))
for(l in 1:length(init.val)){
  Q.init[l]<-optim(init.val[l], Q, method="Brent",lower=-1e2,upper=1e2)$value
}
tet_est <- optim(init.val[which.min(Q.init)[[1]]], Q, method="Brent",lower=-1e2,upper=1e2)$par
tet_se <- sqrt(as.numeric(solve(t(G)%*%solve(Om(tet_est))%*%G)))
J <- as.numeric(Q(tet_est))

tet0 <- 0
delG <- (t(lambda)%*%SigX%*%lambda)*tet0
delGG <- t(lambda)%*%SigX%*%lambda
g0 <- g(tet0)
Om <-  t(lambda)%*%(SigY + (tet0^2)*SigX)%*%lambda
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

res.list <- list("tet_est"=tet_est, "tet_se"=tet_se, "J"=J, "AR.p"=AR.pval, "LM.p"=LM.pval, "CLR.p"=CLR.pval)
return(res.list)
}
