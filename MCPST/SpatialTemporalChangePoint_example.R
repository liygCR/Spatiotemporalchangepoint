library(R.matlab)

##  locations.mat, yst.mat, xst.mat are extract from the Matlab Spatial Statistics Toolbox
p_t <- readMat("locations.mat")
p_t <- p_t$tdums
y <- readMat("yst.mat")
y <- y$yst
x <- readMat("xst.mat")
x <- x$xst
xorig <- x
pt <- c(476,colSums(p_t)) 

x[,11] <- x[,11]-x[,10]
x[,12] <- x[,12]-x[,10]
x <- x[,-10]

##############
# the spatial panel data model:
#   Y - T*Y = phis*S*Y + phist*S*T*Y + phits*T*S*Y + (I-T)*X*beta1 + S*(I-T)*X*beta2 + eps
# with restriction phis=-(phist+phits)
##############
adgLasso_SpatioTemporal_example1<-function( xorig, x, y, pt )
{
  # adaptive group lasso from gglasso with BIC stopping rule 
  
  # Arguments
  # y:  the a dinamic penal data response with a list of length n  .
  # x:  a list of length n, which is the number of observations,each list
  #     contains a pxd  matrix. p is the number of series, d is the number 
  #     of covariate, which includs the intercept.
  # pt:	the number of series in each time
  # Suggest lambda.min=0.001 for p > d, it could be others depends on application.
  
  # Value
  # adgLasso_DP_panel_cp returns an object of class "adgLasso_DinamicPanel_cp" .
  # An object of class "adgLasso_DinamicPanel_cp" is a list containing the following components: 
  #
  # object:	the adaptive group lasso object based on the R function gglasso
  # coeff: a list of the estimated coefficients of the panel model 
  # mse: MSE of the fitted model
  # bic: bic of different lambda
  # lambda.1bic: value of lambda that gives minimum bic 
  # x.ind: index of nonzero coeff
  # ccp.local:	location of the common change points. 
  
  require(gglasso)
  require(Matrix)
  #penalty <- match.arg(penalty)
  
  d <- ncol(x)
  n <- length(pt)
  xorig <- xorig
  
  #   if (q<=p) 
  #     warning(paste("the the design matrix is undefined, ncol should not larger than nrow"))
  
  X_tmp <- NULL
  X_tmp[[1]] <- x
  tmp <- x
  ind_tmp <- c(0, Reduce("+", pt, accumulate = TRUE) )
  for (l in 1:(n-1) ){
    tmp[c(1:ind_tmp[l+1]) ,] <- 0
    X_tmp[[l+1]] <- tmp
  }   
  X <- do.call("cbind",X_tmp)  ## X is a design matrix after transformation
  # sX <- Matrix(X, sparse = TRUE)  ## store X as a sparse Matrix
  
  #get the initial estimate of beta by least square method
  theta.ls<-matrix(0,n*d,1)
  theta.ls[1:d,] <- lsfit(x[(ind_tmp[1]+1):ind_tmp[2],-1], y[(ind_tmp[1]+1):ind_tmp[2],1])$coef
  for(i in 2:n ) {
    theta.ls[((i-1)*d+1):(i*d),] <- lsfit(x[(ind_tmp[i]+1):ind_tmp[i+1],-1], y[(ind_tmp[i]+1):ind_tmp[i+1],1])$coef - 
      lsfit(x[(ind_tmp[i-1]+1):ind_tmp[i],-1], y[(ind_tmp[i-1]+1):ind_tmp[i],])$coef
  } 
  
  ## get the weight w
  #for nu=1  and k=1
  w<-rep( sqrt( (theta.ls[1:d,] %*% theta.ls[1:d,])/d ), d) 
  #for k=2,....,n
  for(i in 2:n ) {
    s <- NULL
    s[i] <- ( (theta.ls[((i - 1) * d + 1):(i * d)] %*% theta.ls[((i - 1) * d + 1):(i * d)])/d )^(1/2)
    w <- c(w,rep(s[i], d))
    rm(s)
  }
  
  
  xs <- scale(X,center=FALSE,scale=1/w)  #  scales by the weights
  
  
  #do group selection by the R function gglasso()
  index<- rep(1:n,each=d)
  Y <- y
  object <- gglasso(xs, Y, group=index, loss="ls" , intercept=FALSE) 
  
  # get min BIC
  beta.tmp <- coef(object)[-1,]
  theta.hat.tmp <- array(dim = c(d,n,100))
  df.tmp <- matrix(0,8,100)
  for(k in 1:100) {
    theta.hat.tmp[,1,k] <- beta.tmp[1:d,k]
    for(i in 2:n){
      theta.hat.tmp[,i,k] <- beta.tmp[((i-1)*d+1):(i*d),k]
    }
    df.tmp[,k] <- apply(theta.hat.tmp[, ,k], 2, function(x){norm(x,type = "2")} )
  }
  df <- apply(df.tmp, 2, function(t){t%*%(rep(d,n)-1)})
  bic <- c()
  for(i in 1:length(object$lambda)) {
    bic[i] <- log(dim(xs)[1])*( sum(df.tmp[,i]>0)+ df[i] )/dim(xs)[1] + log( t(Y-xs %*% object$beta[,i])%*%(Y-xs %*% object$beta[,i]) /length(Y) )
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
  
  numb <- cumsum(c(sum(pt[1:ccp.local[1]]), sum(pt[(1+ccp.local[1]):ccp.local[2]]), sum(pt[(1+ccp.local[2]):n]) ))
  sumresid2 <- sum((lsfit( x[1:numb[1],-1],y[1:numb[1]])$resid)^2)+
    sum((lsfit( x[(1+numb[1]):numb[2],-1],y[(1+numb[1]):numb[2]])$resid)^2)+      
    sum((lsfit( x[(1+numb[2]):numb[3],-1],y[(1+numb[2]):numb[3]])$resid)^2)
  R2 <- 1 - sumresid2/( t(Y-mean(Y))%*%(Y-mean(Y)) )
  
  library(lrmest)
  theta.hat <- NULL
  theta.hat[[1]] <- rls(y[1:numb[1]]~xorig[1:numb[1],-1],r=0,R=c(rep(0,9),1,1,1), 
                        delt = 0)$`*****Restricted Least Square Estimator*****`[,1]
  theta.hat[[2]] <- rls(y[(1+numb[1]):numb[2]]~xorig[(1+numb[1]):numb[2],-1],r=0,R=c(rep(0,9),1,1,1), 
                        delt = 0)$`*****Restricted Least Square Estimator*****`[,1]
  theta.hat[[3]] <- rls(y[(1+numb[2]):numb[3]]~xorig[(1+numb[2]):numb[3],-1],r=0,R=c(rep(0,9),1,1,1), 
                        delt = 0)$`*****Restricted Least Square Estimator*****`[,1]
    
  t_value <- NULL
  t_value[[1]] <- theta.hat[[1]]/rls(y[1:numb[1]]~xorig[1:numb[1],-1],r=0,R=c(rep(0,9),1,1,1), 
                                     delt = 0)$`*****Restricted Least Square Estimator*****`[,2]
  t_value[[2]] <- theta.hat[[2]]/rls(y[(1+numb[1]):numb[2]]~xorig[(1+numb[1]):numb[2],-1],r=0,R=c(rep(0,9),1,1,1), 
                                     delt = 0)$`*****Restricted Least Square Estimator*****`[,2]
  t_value[[3]] <- theta.hat[[3]]/rls(y[(1+numb[2]):numb[3]]~xorig[(1+numb[2]):numb[3],-1],r=0,R=c(rep(0,9),1,1,1), 
                                     delt = 0)$`*****Restricted Least Square Estimator*****`[,2]
  
  p_value <- NULL
  for(i in 1:(length(group.ind)-1) ){
    p_value[[i]] <- 2*stats::pt(- abs(t_value[[i]]), df= sum(pt[group.ind[i]:(group.ind[i+1]-1)]) - 
                                  (group.ind[i+1]-group.ind[i])*d)
  } 
  p_value[[length(group.ind)]] <- 2*stats::pt(-abs(t_value[[length(group.ind)]]), 
                                              df= sum(pt[group.ind[length(group.ind)]:n]) - (group.ind[i+1]-group.ind[i])*d )
  
  # this next line just finds the variable id of coeff. not equal 0
  if(st>0) x.ind<-as.vector(which(coeff !=0)) else x.ind <- 0 
  return(list(object=object, coeff=theta.hat[group.ind], t_value=t_value, p_value=p_value,
              ccp.local= ccp.local, R2 = R2, SSE=sumresid2 ))
} 

object <- adgLasso_SpatioTemporal_example1(xorig, x, y, pt)
ccp.local <- object$ccp.local
theta.hat <- object$coeff
t_value <- object$t_value
p_value <- object$p_value
R2 <- object$R2
SSE <- object$SSE

plot(y,type = "l", xaxt="n", xlab = "Time", ylab = "", main = "The Baton Rouge real estate price data" )
axis(1, at = c(1,cumsum(pt)), c(seq(1985,1993)))
change <- cumsum( c(sum(pt[1:ccp.local[1]]), sum(pt[(1+ccp.local[1]):ccp.local[2]]) ) )
abline(v = change, lty = 2, lwd=2 )

############ model5 in Pace al et.(2000) #####################
library(lrmest)
beta.hat <- rls(y ~ xorig[,-1],r=0,R=c(rep(0,9),1,1,1), 
                  delt = 0)$`*****Restricted Least Square Estimator*****`[,1]
t_value1 <- beta.hat/rls(y ~ xorig[,-1],r=0,R=c(rep(0,9),1,1,1), 
                delt = 0)$`*****Restricted Least Square Estimator*****`[,2]
SSE1 <- sum((y-xorig %*% beta.hat)^2)
R2 <- 1 - sum((y-xorig %*% beta.hat)^2)/( t(y-mean(y))%*%(y-mean(y)) )

####### prediction compared ############
### prediction by pooled model
xnew <- xorig[5229:5243,]
ynew <- y[5229:5243,]

beta.hat.pool <- rls(y[1:5228,] ~ xorig[1:5228,-1],r=0,R=c(rep(0,9),1,1,1), 
                delt = 0)$`*****Restricted Least Square Estimator*****`[,1]

y.pred.pool <- xnew %*% beta.hat.pool
MASE.pool <- mean( abs(ynew - y.pred.pool)/mean(abs(ynew- mean(ynew)) ) )

MAPE.pool <- mean( abs( (ynew - y.pred.pool)/ynew)  )
wMAPE.pool <- sum(abs(ynew - y.pred.pool))/ sum(abs(ynew))



### based on the last period
beta.hat.cp <- rls(y[4474:5228,] ~ xorig[4474:5228,-1],r=0,R=c(rep(0,9),1,1,1), 
                delt = 0)$`*****Restricted Least Square Estimator*****`[,1]

y.pred.cp <- xnew %*% beta.hat.cp
MASE.cp <- mean( abs(ynew - y.pred.cp)/mean(abs(ynew- mean(ynew)) ) )

MAPE.cp <- mean( abs( (ynew - y.pred.cp)/ynew)  )
wMAPE.cp <- sum(abs(ynew - y.pred.cp))/ sum(abs(ynew))


############################################################################
# the spatial panel data model:
#   Y = phit*T*Y + phis*S*Y + phist*S*T*Y + phits*T*S*Y + X*beta1 + eps 
############################################################################
adgLasso_SpatioTemporal_example2<-function( x, y, pt )
{
  # adaptive group lasso from gglasso with BIC stopping rule 
  
  # Arguments
  # y:  the a dinamic penal data response with a list of length n  .
  # x:  a list of length n, which is the number of observations,each list
  #     contains a pxd  matrix. p is the number of series, d is the number 
  #     of covariate, which includs the intercept.
  # pt:	the number of series in each time
  # Suggest lambda.min=0.001 for p > d, it could be others depends on application.
  
  # Value
  # adgLasso_DP_panel_cp returns an object of class "adgLasso_DinamicPanel_cp" .
  # An object of class "adgLasso_DinamicPanel_cp" is a list containing the following components: 
  #
  # object:	the adaptive group lasso object based on the R function gglasso
  # coeff: a list of the estimated coefficients of the panel model 
  # mse: MSE of the fitted model
  # bic: bic of different lambda
  # lambda.1bic: value of lambda that gives minimum bic 
  # x.ind: index of nonzero coeff
  # ccp.local:	location of the common change points. 
  
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
    tmp[c(1:ind_tmp[l+1]) ,] <- 0
    X_tmp[[l+1]] <- tmp
  }   
  X <- do.call("cbind",X_tmp)  ## X is a design matrix after transformation
  # sX <- Matrix(X, sparse = TRUE)  ## store X as a sparse Matrix
  
  #get the initial estimate of beta by least square method
  theta.ls<-matrix(0,n*d,1)
  theta.ls[1:d,] <- lsfit(x[(ind_tmp[1]+1):ind_tmp[2],-1], y[(ind_tmp[1]+1):ind_tmp[2],1])$coef
  for(i in 2:n ) {
    theta.ls[((i-1)*d+1):(i*d),] <- lsfit(x[(ind_tmp[i]+1):ind_tmp[i+1],-1], y[(ind_tmp[i]+1):ind_tmp[i+1],1])$coef - 
      lsfit(x[(ind_tmp[i-1]+1):ind_tmp[i],-1], y[(ind_tmp[i-1]+1):ind_tmp[i],])$coef
  } 
  
  ## get the weight w
  #for nu=1  and k=1
  w<-rep( sqrt( (theta.ls[1:d,] %*% theta.ls[1:d,])/d ), d) 
  #for k=2,....,n
  for(i in 2:n ) {
    s <- NULL
    s[i] <- ( (theta.ls[((i - 1) * d + 1):(i * d)] %*% theta.ls[((i - 1) * d + 1):(i * d)])/d )^(1/2)
    w <- c(w,rep(s[i], d))
    rm(s)
  }
  
  
  xs <- scale(X,center=FALSE,scale=1/w)  #  scales by the weights
  
  
  #do group selection by the R function gglasso()
  index<- rep(1:n,each=d)
  Y <- y
  object <- gglasso(xs, Y, group=index, loss="ls" , intercept=FALSE) 
  
  # get min BIC
  #   beta.tmp <- coef(object)[-1,]
  #   theta.hat.tmp <- array(dim = c(d,n,100))
  #   df.tmp <- matrix(0,8,100)
  #   for(k in 1:100) {
  #     theta.hat.tmp[,1,k] <- beta.tmp[1:d,k]
  #     for(i in 2:n){
  #       theta.hat.tmp[,i,k] <- beta.tmp[((i-1)*d+1):(i*d),k]
  #     }
  #     df.tmp[,k] <- apply(theta.hat.tmp[, ,k], 2, function(x){norm(x,type = "2")} )
  #   }
  #   df <- apply(df.tmp, 2, function(t){t%*%(rep(d,n)-1)})
  #   bic <- c()
  #   for(i in 1:length(object$lambda)) {
  #     bic[i] <- log(dim(xs)[1])*( sum(df.tmp[,i]>0)+ df[i] )/dim(xs)[1] + log( t(Y-xs %*% object$beta[,i])%*%(Y-xs %*% object$beta[,i]) /length(Y) )
  #   }
  
  # get min BIC
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
  
  numb <- cumsum(c(sum(pt[1:ccp.local[1]]), sum(pt[(1+ccp.local[1]):ccp.local[2]]), sum(pt[(1+ccp.local[2]):n]) ))
  sumresid2 <- sum((lsfit( x[1:numb[1],-1],y[1:numb[1]])$resid)^2)+
    sum((lsfit( x[(1+numb[1]):numb[2],-1],y[(1+numb[1]):numb[2]])$resid)^2)+      
    sum((lsfit( x[(1+numb[2]):numb[3],-1],y[(1+numb[2]):numb[3]])$resid)^2)
  R2 <- 1 - sumresid2/( t(Y-mean(Y))%*%(Y-mean(Y)) )
  
  theta.hat <- NULL
  theta.hat[[1]] <- summary(lm( Y[1:numb[1]]~x[1:numb[1],-1]))$coefficients[,1]
  theta.hat[[2]] <- summary(lm( Y[(1+numb[1]):numb[2]]~x[(1+numb[1]):numb[2],-1]))$coefficients[,1]
  theta.hat[[3]] <- summary(lm( Y[(1+numb[2]):numb[3]]~x[(1+numb[2]):numb[3],-1]))$coefficients[,1]
  
  SE <- NULL
  SE[[1]] <- summary(lm(Y[1:numb[1]] ~ x[1:numb[1],-1]))$coefficients[, 2]
  SE[[2]] <- summary(lm(Y[(1 + numb[1]):numb[2]] ~ x[(1 + numb[1]):numb[2],-1]))$coefficients[, 2]
  SE[[3]] <- summary(lm(Y[(1 + numb[2]):numb[3]] ~ x[(1 + numb[2]):numb[3],-1]))$coefficients[, 2]
  
  t_value <- NULL
  t_value[[1]] <- summary(lm( Y[1:numb[1]] ~ x[1:numb[1],-1]))$coefficients[,3]
  t_value[[2]] <- summary(lm( Y[(1+numb[1]):numb[2]] ~ x[(1+numb[1]):numb[2],-1]))$coefficients[,3]
  t_value[[3]] <- summary(lm( Y[(1+numb[2]):numb[3]] ~ x[(1+numb[2]):numb[3],-1]))$coefficients[,3]
  
  p_value <- NULL
  p_value[[1]] <- summary(lm( Y[1:numb[1]] ~ x[1:numb[1],-1]))$coefficients[,4]
  p_value[[2]] <- summary(lm( Y[(1+numb[1]):numb[2]] ~ x[(1+numb[1]):numb[2],-1]))$coefficients[,4]
  p_value[[3]] <- summary(lm( Y[(1+numb[2]):numb[3]] ~ x[(1+numb[2]):numb[3],-1]))$coefficients[,4]

  
  # this next line just finds the variable id of coeff. not equal 0
  if(st>0) x.ind<-as.vector(which(coeff !=0)) else x.ind <- 0 
  return(list(object=object, coeff=theta.hat[group.ind], t_value=t_value, p_value=p_value, SE=SE,
              ccp.local= ccp.local, R2 = R2, SSE=sumresid2 ))
} 


p_t <- readMat("locations.mat")
p_t <- p_t$tdums
y <- readMat("y2.mat")
y <- y$yst2
x <- readMat("xst2.mat")
x <- x$xst2

pt <- c(476,colSums(p_t)) 

object <- adgLasso_SpatioTemporal_example2(x, y, pt)
ccp.local <- object$ccp.local
theta.hat <- object$coeff
t_value <- object$t_value
p_value <- object$p_value
R2 <- object$R2
SSE <- object$SSE

plot(y,type = "l", xaxt="n", xlab = "Time", ylab = "", main = "House Data in Baton Rouge" )
axis(1, at = c(1,cumsum(pt)), c(seq(1985,1993)))
change <- cumsum(c(sum(pt[1:ccp.local[1]]), sum(pt[(1+ccp.local[1]):ccp.local[2]]) ))
abline(v = change, lty = 1, lwd = 2 )

######### OLS without common breaks ###############
summary(lm( y ~ x[,-1]))

##SSE
sum( lm( y ~ x[,-1])$res^2 )

################### prediction compared ########################
### prediction by pooled model
xnew2 <- x[5234:5243,]
ynew2 <- y[5234:5243,]

fit2 <- lm(y[1:5233,] ~ x[1:5233,-1] )
summary(fit2)

y.pred.pool2 <- xnew2 %*% fit2$coeff
MASE.pool2 <- mean( abs(ynew2 - y.pred.pool2)/mean(abs(ynew2 - mean(ynew2)) ) )

MAPE.pool2 <- mean( abs( (ynew2 - y.pred.pool2)/ynew2 )  )
wMAPE.pool2 <- sum(abs(ynew2 - y.pred.pool2))/ sum(abs(ynew2))


### based on the last period
fit.cp2 <- lm(y[3374:5233,] ~ x[3374:5233,-1] )
summary(fit.cp2)

y.pred.cp2 <- xnew2 %*% fit.cp2$coeff
MASE.cp2 <- mean( abs(ynew2 - y.pred.cp2)/mean(abs(ynew2 - mean(ynew2)) ) )

MAPE.cp2 <- mean( abs( (ynew2 - y.pred.cp2)/ynew2)  )
wMAPE.cp2 <- sum(abs(ynew2 - y.pred.cp2))/ sum(abs(ynew2))
