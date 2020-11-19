######################################################################
##############
# the spatial panel data model:
#   Y - T*Y = phis*S*Y + phist*S*T*Y + phits*T*S*Y + (I-T)*X*beta1 + S*(I-T)*X*beta2 + eps ; i=1,...,N ; t=2,...T
##############
adgLasso_SpatioTemporal_cp<-function( x, y, pt )
{
  # adaptive group lasso from gglasso with BIC stopping rule 
  
  # Arguments
  # y:  the a dinamic penal data response with a list of length n  .
  # x:  a list of length n, which is the number of observations,each list
  #     contains a pxd  matrix. p is the number of series, d is the number 
  #     of covariate.
  # pt: number of observations in each time t
  # penalty:	the method to be used by group Lasso,group mcp or group scad.
  # Suggest lambda.min=0.001 for p > d, it could be others depends on application.
  
  # Value
  # adgLasso_SpatioTemporal_cp returns an object of class "adgLasso_SpatioTemporal_cp" .
  # An object of class "adgLasso_SpatioTemporal_cp" is a list containing the following components: 
  #
  # object:	the adaptive group lasso object based on the R function gglasso
  # coeff: a list of the estimated coefficients of the panel model 
  # mse: MSE of the fitted model
  # bic: bic of different lambda
  # inv.Gamma: inverse covariance matrix
  # cov.mat: estimated covariance matrix
  # x.ind: index of nonzero coeff
  # ccp.local:	location of the common change points. 
  # p_value: the p_value of the estimated cofficients
  # t_value: t test of the estimation
  
  require(gglasso)
  require(Matrix)
  #penalty <- match.arg(penalty)
  
  d <- ncol(x)
  n <- length(pt)
  
  #   if (q<=p) 
  #     warning(paste("the the design matrix is undefined, ncol should not larger than nrow"))
  
  X_tmp <- NULL
  X_tmp[[1]] <- x
  tmp <- x
  ind_tmp <- c(0, Reduce("+", pt, accumulate = TRUE) )
  for (l in 1:(n-1) ){
    tmp[c(1:ind_tmp[l+1]),] <- 0
    X_tmp[[l+1]] <- tmp
  }   
  X <- do.call("cbind",X_tmp)  ## X is a design matrix after transformation
  # sX <- Matrix(X, sparse = TRUE)  ## store X as a sparse Matrix
  
  #get the initial estimate of beta by least square method
  theta.ls<-matrix(0,n*d,1)
  theta.ls[1:d,] <- lsfit(x[(ind_tmp[1]+1):ind_tmp[2],], y[(ind_tmp[1]+1):ind_tmp[2],1], intercept = F)$coef
  for(i in 2:n ) {
    theta.ls[((i-1)*d+1):(i*d),] <- lsfit(x[(ind_tmp[i]+1):ind_tmp[i+1],], y[(ind_tmp[i]+1):ind_tmp[i+1],1], intercept = F)$coef - 
      lsfit(x[(ind_tmp[i-1]+1):ind_tmp[i],], y[(ind_tmp[i-1]+1):ind_tmp[i],], intercept = F)$coef
  } 
  
  ## get the weight w
  #for nu=1  and k=1
  w<-rep( sqrt( (theta.ls[1:d,] %*% theta.ls[1:d,])/d ), d) 
  #for k=2,....,n
  for(i in 2:n ) {
    s <- NULL
    s[i] <- sqrt( (theta.ls[((i - 1) * d + 1):(i * d)] %*% theta.ls[((i - 1) * d + 1):(i * d)])/d )
    w <- c(w,rep(s[i], d))
    rm(s)
  }
  
  xs <- scale(X,center=FALSE,scale=1/w)  #  scales by the weights
  
  #do group selection by the R function gglasso()
  index<- rep(1:n,each=d)
  Y <- y
  object <- gglasso(xs, Y, group=index, loss="ls" , intercept=FALSE) 
  
  # get min BIC
  # bic=log(n)*object$df+n*log(RSS/n)   # rss/n version
  bic <- c()
  for(i in 1:length(object$lambda)) {
    bic[i] <- log(dim(xs)[1])*object$df[i] + dim(xs)[1]*log( t(Y-xs %*% object$beta[,i])%*%(Y-xs %*% object$beta[,i]) /length(Y) )
  }
  
  step.bic=which.min(bic)            # step with min BIC
  lambda.1bic=object$lambda[step.bic]
  grlasso.coef <- coef(object, s =lambda.1bic)[-1]
  coeff <- w*grlasso.coef                # get back in right scale
  names(coeff) <- object$group
  coeff_not0 <- coeff[which(coeff!=0)]
  st <- sum(coeff !=0)         # number nonzero
  
  group.ind <- as.numeric(unique(names(coeff_not0)))  
  ccp.local<- group.ind[-1] #omit the 1st group
  
  fit<-predict(object, xs, s=lambda.1bic, type="link")
  names(fit) <- rep(1:n, times = pt)
  names(Y) <- rep(1:n, times = pt)
  mse=sum((Y-fit)^2)/(sum(pt)-st)   ## model fit MSE
  
  theta.hat <- NULL
  theta.hat[[1]] <- coeff[1:d]
  for(i in 2:n){
    theta.hat[[i]] <- coeff[((i-1)*d+1):(i*d)] + theta.hat[[i-1]]
  }
  
  xlist <- NULL
  for(i in 1:n){
    xlist[[i]] <- x[(ind_tmp[i]+1):ind_tmp[i+1],]
  }
  
  Gamma <- NULL;sigma2 <- NULL
  for (i in 1:(length(group.ind)-1)) {
    tmp <- sapply(xlist[group.ind[i]:(group.ind[i+1]-1)],
                  function(a,b){(t(a)%*%a)/(nrow(a)*b)}, b=group.ind[i+1]-group.ind[i], simplify=F)
    Gamma[[i]] <- Reduce("+", tmp)
    sigma2[i] <- sum((Y[which(names(Y)==i)]-fit[which(names(fit)==i)])^2)/(length(Y[which(names(Y)==i)])-d)
  }
  tmp <- sapply(xlist[group.ind[length(group.ind)]:n],
                function(a,b){(t(a)%*%a)/(nrow(a)*b)}, b=n+1-group.ind[length(group.ind)], simplify=F)
  Gamma[[length(group.ind)]] <- Reduce("+", tmp)
  sigma2[length(group.ind)] <- sum((Y[which(names(Y)==length(group.ind))]-
                                      fit[which(names(fit)==length(group.ind))])^2)/(length(Y[which(names(Y)==length(group.ind))])-d)
  rm(tmp)
  ## get the covariance matrix of theta
  inv.Gamma <- sapply(Gamma, solve, simplify = F)  
  cov.mat <- mapply(function(a,b){a*b}, sigma2, inv.Gamma, SIMPLIFY = F)
  
  ## get the t-value for significient test
  t_value <- NULL
  p_value <- NULL
  for(i in 1:(length(group.ind)-1) ){
    t_value[[i]] <- theta.hat[[i]]/sqrt(diag(cov.mat[[i]])/sum(pt[group.ind[i]:(group.ind[i+1]-1)]))
    p_value[[i]] <- 2*stats::pt(- abs(t_value[[i]]), df= sum(pt[group.ind[i]:(group.ind[i+1]-1)]) - 
                                  (group.ind[i+1]-group.ind[i])*d)
  } 
  t_value[[length(group.ind)]] <- theta.hat[[length(group.ind)]]/sqrt(diag(cov.mat[[length(group.ind)]])/
                                                                        sum(pt[group.ind[length(group.ind)]:n])) 
  p_value[[length(group.ind)]] <- 2*stats::pt(-abs(t_value[[length(group.ind)]]), 
                                              df= sum(pt[group.ind[length(group.ind)]:n]) - (group.ind[i+1]-group.ind[i])*d )
  
  # this next line just finds the variable id of coeff. not equal 0
  if(st>0) x.ind<-as.vector(which(coeff !=0)) else x.ind <- 0 
  return(list(object=object, coeff=theta.hat[group.ind], cov.mat=cov.mat, t_value=t_value, p_value=p_value,
              sigma2=sigma2, mse=mse, cov.mat=cov.mat, x.ind=x.ind, ccp.local= ccp.local,bic=bic))
} 

########## simulation fun #################
sim.fun<-function(N, d, S, ST, TS, lagIT, lagSIT, invmat, b1, b2, sd=0.1 ) {  
  
  x <- matrix(rnorm(N*d, 0,1), N, d)
  lagTXb <- matrix(0,N,1)
  lagSTXb <- matrix(0,N,1)
  a1 <- as.matrix(lagIT)
  a2 <- as.matrix(lagSIT)
  for(k in 1:n)
  {
    lagTXb[((k-1)*p+1):(k*p),] <- (a1%*%x)[((k-1)*p+1):(k*p),]%*%b1[,k]
    lagSTXb[((k-1)*p+1):(k*p),] <- (a2%*%x)[((k-1)*p+1):(k*p),]%*%b2[,k]
  }
  rm(a1,a2)
  lagTXb <- Matrix(lagTXb, sparse = T)
  lagSTXb <- Matrix(lagSTXb, sparse = T)
  y <- invmat%*%(lagTXb + lagSTXb + rnorm(N,0,sd))
  yt <- y- Tmat %*% y
  yt <- as.matrix(yt)
  z <- cbind2(cbind2(lagIT%*%x,lagSIT%*%x), cbind2(cbind2(S%*%y, ST%*%y),TS%*%y))
  z <- as.matrix(z)
  pt <- rep(p,n)
  
  #select by BIC
  sim.object <- adgLasso_SpatioTemporal_cp(x=z, y=yt, pt=pt)
  sim.local <- sim.object$ccp.local
  sim.coeff <- sim.object$coeff
  sim.cov <- sim.object$cov.mat
  sim.t_value <- sim.object$t_value
  sim.p_value <- sim.object$p_value
  #   sim.t_value <- sim.object$t_value
  #   sim.p_value <- sim.object$p_value
  
  cat(paste("Common Change-point Location:",sim.local,"\n"),"\n")
  return(list(sim.local=sim.local, sim.cov=sim.cov, sim.coeff=sim.coeff,
              sim.t_value= sim.t_value, sim.p_value= sim.p_value))
  
}

ptm <- proc.time()
sim <-replicate( 10, sim.fun(p, n, d, I, Tmat, S, ST, TS, phis, phist,
                                           phits, b1, b2, sd=.1) ) # sd=.1
proc.time() - ptm

library(foreach)
library(doMC)
registerDoMC(20)
ptm <- proc.time()
r <- foreach(1:500) %dopar% {
   sim.fun(N, d, S, ST, TS, lagIT, lagSIT, invmat, b1, b2, sd = 0.1)
}
proc.time() - ptm


###############  get spatiotemporal matrix ########################
p = 100  
n = 20   
d = 4  
N = n*p

library(Matrix)
set.seed(111)
D <- matrix(0,N,N)
xcoord <- rnorm(N)
ycoord <- rnorm(N)
for (j in 1:(N-1)) {
  for(i in (j+1):(N)){
    D[i,j] <- (xcoord[i]-xcoord[j])^2+(ycoord[i]-ycoord[j])^2
  }
} 
# D <- Matrix(D, sparse = T)

# ms=5
# rho <- 0.8
# lagcoef <- rho^(1:ms)
# # %standardizing these to sum to 1
# lagcoef <- lagcoef/sum(lagcoef)
# lagcoef <- as.list(lagcoef)


# %weighting and summing to form spatial average
ms <- 5
S <- NULL
for(i in 1:ms){
  S[[i]] <- matrix(0,N,N)
  for(j in (i+1):(N)){
    S[[i]][j, which(order(D[j,1:(j-1)])==i)] <- 1
  }
}
S <- Reduce("+", S)
for(i in 2:(N)){
  S[i,] <- S[i,]/sum(S[i,])
}
rm(D)
S <- Matrix(S, sparse = T)

mt=60
Tmat <- matrix(0,N,N)
for(i in 2:(mt+1)) {
  Tmat[i,1:(i-1)] <- rep(1/(i-1),i-1)
}
for(i in (mt+2):N) {
  Tmat[i,(i-mt):(i-1)] <- rep(1/(mt),mt)
}
Tmat <- Matrix(Tmat, sparse = T)

# tri <- 200
# for(i in (tri+1):N) {
#   S[i:N, 1:(i-tri)] <- 0
# }

ST <- S%*%Tmat
TS <- Tmat%*%S
I <- Diagonal(N)
lagIT <- (I- Tmat)
lagSIT <- S%*%(I- Tmat)

# #### bands matrix ######
# tri <- p
# for(i in (tri+1):N) {
#   S[i:N, 1:(i-tri)] <- 0
# }
# S <- Matrix(S, sparse = T)
# 
# ST <- S%*%Tmat
# TS <- Tmat%*%S
# I <- diag(N)
# lagIT <- (I- Tmat)
# lagSIT <- S%*%(I- Tmat)
# ST <- Matrix(ST, sparse = T)
# TS <- Matrix(ST, sparse = T)



###########################   no common break  ################################
phis <- rep(0.9, n)
phist <- rep(-0.5, n)
phits <- rep(-0.4, n)

set.seed(111)
b1 <-matrix(1,d,n)

b2 <-matrix(1.5,d,n)

a <- as.matrix(S)
lagS <- Matrix(0,N,N)
for(k in 1:n)
{
  lagS[((k-1)*p+1):(k*p),] <- phis[k]*a[((k-1)*p+1):(k*p),]
}
rm(a)
lagS <- Matrix(lagS, sparse = T)

a <- as.matrix(ST)
lagST <- matrix(0,N,N)
for(k in 1:n)
{
  lagST[((k-1)*p+1):(k*p),] <- phist[k]* a[((k-1)*p+1):(k*p),]
}
rm(a)
lagST <- Matrix(lagST, sparse = T)

a <- as.matrix(TS)
lagTS <- matrix(0,N,N)
for(k in 1:n)
{
  lagTS[((k-1)*p+1):(k*p),] <- phits[k]*a[((k-1)*p+1):(k*p),]
}
rm(a)
lagTS <- Matrix(lagTS, sparse = T)

invmat <- Matrix::solve(I- Tmat- lagS- lagST- lagTS)


###########################   single common break  ################################
ccp_tr <- n/2 # location of change point
# N1 <- (ccp_tr-1)*p
# N2 <- (n-ccp_tr+1)*p

phis <- c(rep(0.9, ccp_tr-1), rep(0.7, n-ccp_tr+1) )
phist <- c(rep(-0.5, ccp_tr-1), rep(-0.4, n-ccp_tr+1) )
phits <- c(rep(-0.4, ccp_tr-1), rep(-0.3, n-ccp_tr+1) )

set.seed(111)
b1 <-matrix(0,d,n)
b1[,1:(ccp_tr-1)] <- matrix(1,d,ccp_tr-1)
b1[,ccp_tr:n] <- b1[,ccp_tr-1] + runif(d,-2,2)

b2 <-matrix(0,d,n)
b2[,1:(ccp_tr-1)] <- matrix(1,d,ccp_tr-1)
b2[,ccp_tr:n] <- b2[,ccp_tr-1] + runif(d,-2,2)

a <- as.matrix(S)
lagS <- Matrix(0,N,N)
for(k in 1:n)
{
  lagS[((k-1)*p+1):(k*p),] <- phis[k]*a[((k-1)*p+1):(k*p),]
}
rm(a)
lagS <- Matrix(lagS, sparse = T)

a <- as.matrix(ST)
lagST <- matrix(0,N,N)
for(k in 1:n)
{
  lagST[((k-1)*p+1):(k*p),] <- phist[k]* a[((k-1)*p+1):(k*p),]
}
rm(a)
lagST <- Matrix(lagST, sparse = T)

a <- as.matrix(TS)
lagTS <- matrix(0,N,N)
for(k in 1:n)
{
  lagTS[((k-1)*p+1):(k*p),] <- phits[k]*a[((k-1)*p+1):(k*p),]
}
rm(a)
lagTS <- Matrix(lagTS, sparse = T)

invmat <- Matrix::solve(I- Tmat- lagS- lagST- lagTS)

######################### 5 common breaks   ##########################
ccp_tr <- c( floor(n/6), floor(n/3), floor(n/2), floor(2*n/3), floor(5*n/6) )   ## location of true change point

phis <- c(rep(0.8, ccp_tr[1]-1), rep(0.6, ccp_tr[2]- ccp_tr[1]), rep(0.9, ccp_tr[3]- ccp_tr[2]),
  rep(0.7, ccp_tr[4]- ccp_tr[3]),rep(0.9, ccp_tr[5]- ccp_tr[4]),rep(0.65, n- ccp_tr[5]+1))
phist <- c(rep(-0.4, ccp_tr[1]-1), rep(-0.3, ccp_tr[2]- ccp_tr[1]), rep(-0.45, ccp_tr[3]- ccp_tr[2]),
          rep(-0.35, ccp_tr[4]- ccp_tr[3]),rep(-0.45, ccp_tr[5]- ccp_tr[4]),rep(-0.35, n- ccp_tr[5]+1))
phits <- c(rep(-0.4, ccp_tr[1]-1), rep(-0.3, ccp_tr[2]- ccp_tr[1]), rep(-0.45, ccp_tr[3]- ccp_tr[2]),
           rep(-0.3, ccp_tr[4]- ccp_tr[3]),rep(-0.45, ccp_tr[5]- ccp_tr[4]),rep(-0.3, n- ccp_tr[5]+1))

set.seed(111)
b1 <-matrix(0,d,n)
b1[,1:(ccp_tr[1]-1)] <- runif(d,-1.5,1.5)
for(i in 1:(length(ccp_tr)-1) ) {
  b1[,(ccp_tr[i]):(ccp_tr[i+1]-1)] <- b1[,ccp_tr[i]-1] +  runif(d,-0.5,0.5)
}
b1[,ccp_tr[length(ccp_tr)]:n] <- b1[,ccp_tr[length(ccp_tr)]-1] + runif(d,-0.5,0.5)

b2 <-matrix(0,d,n)
b2[,1:(ccp_tr[1]-1)] <- runif(d,-1.5,1.5)
for(i in 1:(length(ccp_tr)-1) ) {
  b2[,(ccp_tr[i]):(ccp_tr[i+1]-1)] <- b2[,ccp_tr[i]-1] + runif(d,-0.5,0.5)
}
b2[,ccp_tr[length(ccp_tr)]:n] <- b2[,ccp_tr[length(ccp_tr)]-1] + runif(d,-0.5,0.5)

a <- as.matrix(S)
lagS <- matrix(0,N,N)
for(k in 1:n)
{
  lagS[((k-1)*p+1):(k*p),] <- phis[k]*a[((k-1)*p+1):(k*p),]
}
rm(a)
lagS <- Matrix(lagS, sparse = T)

a <- as.matrix(ST)
lagST <- matrix(0,N,N)
for(k in 1:n)
{
  lagST[((k-1)*p+1):(k*p),] <- phist[k]* a[((k-1)*p+1):(k*p),]
}
rm(a)
lagST <- Matrix(lagST, sparse = T)

a <- as.matrix(TS)
lagTS <- matrix(0,N,N)
for(k in 1:n)
{
  lagTS[((k-1)*p+1):(k*p),] <- phits[k]*a[((k-1)*p+1):(k*p),]
}
rm(a)
lagTS <- Matrix(lagTS, sparse = T)

invmat <- Matrix::solve(I- Tmat- lagS- lagST- lagTS)

