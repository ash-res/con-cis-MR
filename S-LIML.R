# Selected LIML tests

# INPUT
# bx: p-vector of genetic variant-exposure associations 
# by: p-vector of genetic variant-outcome associations
# sx: p-vector of standard errors of genetic variant-exposure associations 
# sy: p-vector of standard errors of genetic variant-exposure associations
# ld: (p x p) genetic variant correlation matrix
# tet0: value to test the null hypothesis, H0: true theta = tet0
# d0: value (0,1) to conduct d0-level pre-tests; e.g. d0 = 0.05 corresponds to a 95% confidence level
# r (optional): the number of factors. Default: chooses the number of factors that explain 99% of total variation in the p genetic variants

# OUTPUT
# tet_est: value for the Selected LIML estimator
# factors: the number of factors r* selected by the pre-tests
# thres: equals 0 if at least 1 factor passed the pre-testing threshold of relevance
# N.B. cil and ciu are not the confidence intervals for tet0, but in fact... 
# cil: lower 2.5% tail of the distribution of tet_est under the null H0: true theta is tet0 
# ciu: upper 2.5% tail of the distribution of tet_est under the null H0: true theta is tet0
# result: if tet_est does not lie between cil and ciu then we reject H0


S.LIML <- function(bx,by,sx,sy,ld,tet0,d0,r=NULL){
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

D <- vector(,r)
for (k in 1:r){D[k] <- as.numeric(t(lambda[,k])%*%SigX%*%lambda[,k])}
D <- diag(D)
G <- -t(lambda)%*%bx

T <- vector(,r)
for (k in 1:r){T[k] <- as.numeric((1/sqrt(D[k,k]))*G[k])}
c0 <- qnorm(1-(d0/2))

S <- vector(,r)
for (k in 1:r){
  S[k] <- ifelse(abs(as.numeric(T[k])) > c0,1,0)
}
r0 <- sum(S==1)


if(r0>0){
R <- which(S==1)
gam <- diag(r)[,R]
Om.S <- function(tet){t(gam)%*%t(lambda)%*%(SigY + ((tet^2)*SigX))%*%lambda%*%gam}
G.S <- t(gam)%*%G
g.S <- function(tet){t(gam)%*%t(lambda)%*%(by-(bx*tet))}
Q.S <- function(tet){as.numeric(t(g.S(tet))%*%solve(Om.S(tet))%*%g.S(tet))}
init.val <- seq(-2,2,0.2)
Q.init <- vector(,length=length(init.val))
for(l in 1:length(init.val)){
  Q.init[l]<-optim(init.val[l], Q.S, method="Brent",lower=-1e2,upper=1e2)$value
}
tet_est <- optim(init.val[which.min(Q.init)[[1]]], Q.S, method="Brent",lower=-1e2,upper=1e2)$par
Om.S <- Om.S(tet_est)
V.S <- 1/(as.numeric(t(G.S)%*%solve(Om.S)%*%G.S))

C.G <- - as.vector(sqrt(solve(D))%*%(t(lambda)%*%SigX%*%lambda)%*%gam%*%solve(Om.S)%*%(G.S)%*%V.S*tet_est)
u <- as.vector(T-(C.G*(1/V.S)*(tet_est-tet0)))

ub <- function(Z){as.vector(u+(C.G*(1/V.S)*sqrt(V.S)*Z))}
K <- rnorm(2e3,0,1)
L <- vector(,length=length(K))
for (k in 1:length(K)){L[k] <- (prod((abs(ub(K[k]))[R]>c0)*1)*prod((abs(ub(K[k]))[-R]<=c0)*1))}

w <- seq(-20,20,0.001)
P <- vector(,length=length(w))
for (l in 1:length(w)){P[l] <- mean(ifelse(K<=(w[l]/sqrt(V.S)),1,0)*L)/mean(L)}

cil <- tet0+as.numeric(w[max(which(P<=0.025))])
ciu <- tet0+as.numeric(w[min(which(P>=0.975))])

if(tet_est<=ciu & tet_est>=cil){result <- (paste0("H0: cannot reject a true value of ",tet0, " for a 0.05 level test"))} else  {result <- (paste0("H1: reject a true value of ",tet0, " for a 0.05 level test"))} 

res.list <- list("tet_est"=tet_est, "factors"=r0, "thres"=ifelse(sum(S)==0,1,0), "cil"=cil, "ciu"=ciu, "result" = result)
}

if(r0==0){
  res.list <- list("tet_est"=NA, "factors"=NA, "thres"=NA, "cil"=NA, "ciu"=NA, "result" = NA)
}
return(res.list)
}
