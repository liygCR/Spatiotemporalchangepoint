library(gdata)
library(IDPmisc)
load("HKdist_mat.RData")
load("HKtranpric")

## 1). 2006-01-01 to 2016-01-01
datsub <- subset(dat, Trans.Date > as.Date("2008-01-01") & 
                   Trans.Date < as.Date("2017-01-01"))
datsub$`Price(m)` <- log(datsub$`Price(m)`)
datsub <- NaRV.omit(datsub)

n = nrow(datsub)
Tlag = 90 ### time lag is 90(120, 180, 30) days 

ptm <- proc.time()
W <- matrix(0, n, n)
### weight matrix
for (i in 2:n) {
  for (j in 1:(i-1)) {
    W[i, j] <-
      ifelse(
        datsub$`building name`[i] != datsub$`building name`[j] &
          datsub$Trans.Date[i] - datsub$Trans.Date[j] <= Tlag,
        1 / dist_mat[datsub$`building name`[i], datsub$`building name`[j]],
        ifelse(
          datsub$`building name`[i] == datsub$`building name`[j] &
            datsub$Trans.Date[i] - datsub$Trans.Date[j] <= Tlag,
          1 / min(lowerTriangle(dist_mat)), 0
          )
        )
  }
  print(i)
}
proc.time() - ptm
W <- W/1000

Vnum <- function(dd) { 
  ### Tlag = 90

  v1 <- sapply(dd, function(x, y) {
    sum(as.numeric(y - x) <= 30 &
          as.numeric(y - x) >= 1)
  }, y = datsub$Trans.Date, simplify = T)
  
  v2 <- sapply(dd, function(x, y) {
    sum(as.numeric(y - x) <= 60 &
          as.numeric(y - x) >= 31)
  }, y = datsub$Trans.Date, simplify = T)
  
  v3 <- sapply(dd, function(x, y) {
    sum(as.numeric(y - x) <= 90 &
          as.numeric(y - x) >= 61)
  }, y = datsub$Trans.Date, simplify = T)
  
  V <- v1 + v2/2 + v3/3
  return(V)
  
}

V <-  Vnum(datsub$Trans.Date)

Y <- as.matrix(datsub$`Price(m)`)
X <- cbind( datsub$`building age`, datsub$floor, datsub$`Area (Gross)`, 
            datsub$`building age`^2, datsub$floor^2, datsub$`Area (Gross)`^2)

pt <- apply( as.matrix(2008:2015),1, function(x){nrow(subset(datsub, format.Date(Trans.Date, "%Y")== as.character(x) ))})

x <- apply(X, 2, function(a) {scale(a)})
xall <- cbind(1, W%*%Y, V*(W%*%Y), x)

# y <- as.matrix(scale(Y, scale = FALSE))
y <- as.matrix(Y[- grep("2016", datsub$Trans.Date),])
xx <- xall[- grep("2016", datsub$Trans.Date),]


########
adgLasso_SpatioTemporal_example<-function( x, y, pt )
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
  
  # # get min BIC
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
  
  cl <- ccp.local-1
  numb <- cumsum(c(sum(pt[1:cl[1]]), sum(pt[(1 + cl[1]):cl[2]]),
             sum(pt[(1 + cl[2]):cl[3]]), sum(pt[(1 + cl[3]):n])))
  sumresid2 <- sum((lsfit(x[1:numb[1], -1], y[1:numb[1]])$resid) ^ 2) +
    sum((lsfit(x[(1 + numb[1]):numb[2], -1], y[(1 + numb[1]):numb[2]])$resid) ^
          2) +
    sum((lsfit(x[(1 + numb[2]):numb[3], -1], y[(1 + numb[2]):numb[3]])$resid) ^
          2) +
    sum((lsfit(x[(1 + numb[3]):numb[4], -1], y[(1 + numb[3]):numb[4]])$resid) ^
          2)
  R2 <- 1 - sumresid2 / (t(Y - mean(Y)) %*% (Y - mean(Y)))
  
  theta.hat <- NULL
  theta.hat[[1]] <- summary(lm( Y[1:numb[1]]~x[1:numb[1],-1]))$coefficients[,1]
  theta.hat[[2]] <- summary(lm( Y[(1+numb[1]):numb[2]]~x[(1+numb[1]):numb[2],-1]))$coefficients[,1]
  theta.hat[[3]] <- summary(lm( Y[(1+numb[2]):numb[3]]~x[(1+numb[2]):numb[3],-1]))$coefficients[,1]
  theta.hat[[4]] <- summary(lm( Y[(1+numb[3]):numb[4]]~x[(1+numb[3]):numb[4],-1]))$coefficients[,1]
  
  SE <- NULL
  SE[[1]] <- summary(lm(Y[1:numb[1]] ~ x[1:numb[1],-1]))$coefficients[, 2]
  SE[[2]] <- summary(lm(Y[(1 + numb[1]):numb[2]] ~ x[(1 + numb[1]):numb[2],-1]))$coefficients[, 2]
  SE[[3]] <- summary(lm(Y[(1 + numb[2]):numb[3]] ~ x[(1 + numb[2]):numb[3],-1]))$coefficients[, 2]
  SE[[4]] <- summary(lm(Y[(1 + numb[3]):numb[4]] ~ x[(1 + numb[3]):numb[4],-1]))$coefficients[, 2]
  
  t_value <- NULL
  t_value[[1]] <- summary(lm( Y[1:numb[1]] ~ x[1:numb[1],-1]))$coefficients[,3]
  t_value[[2]] <- summary(lm( Y[(1+numb[1]):numb[2]] ~ x[(1+numb[1]):numb[2],-1]))$coefficients[,3]
  t_value[[3]] <- summary(lm( Y[(1+numb[2]):numb[3]] ~ x[(1+numb[2]):numb[3],-1]))$coefficients[,3]
  t_value[[4]] <- summary(lm( Y[(1+numb[3]):numb[4]] ~ x[(1+numb[3]):numb[4],-1]))$coefficients[,3]
  
  p_value <- NULL
  p_value[[1]] <- summary(lm( Y[1:numb[1]] ~ x[1:numb[1],-1]))$coefficients[,4]
  p_value[[2]] <- summary(lm( Y[(1+numb[1]):numb[2]] ~ x[(1+numb[1]):numb[2],-1]))$coefficients[,4]
  p_value[[3]] <- summary(lm( Y[(1+numb[2]):numb[3]] ~ x[(1+numb[2]):numb[3],-1]))$coefficients[,4]
  p_value[[4]] <- summary(lm( Y[(1+numb[3]):numb[4]] ~ x[(1+numb[3]):numb[4],-1]))$coefficients[,4]
  
  r.square <- c()
  r.square[1] <- summary(lm( Y[1:numb[1]]~x[1:numb[1],-1]))$r.square
  r.square[2] <- summary(lm( Y[(1+numb[1]):numb[2]]~x[(1+numb[1]):numb[2],-1]))$r.square
  r.square[3] <- summary(lm( Y[(1+numb[2]):numb[3]]~x[(1+numb[2]):numb[3],-1]))$r.square
  r.square[4] <- summary(lm( Y[(1+numb[3]):numb[4]]~x[(1+numb[3]):numb[4],-1]))$r.square
  
  # this next line just finds the variable id of coeff. not equal 0
  if(st>0) x.ind<-as.vector(which(coeff !=0)) else x.ind <- 0 
  return(list(object=object, coeff=theta.hat, t_value=t_value, p_value=p_value, SE=SE,
              ccp.local= ccp.local, R2 = R2, SSE=sumresid2, r.square = r.square ))
} 

########

object <- adgLasso_SpatioTemporal_example(xx, y, pt)

ccp.local <- object$ccp.local # 2,4,5
theta.hat <- object$coeff
t_value <- object$t_value
p_value <- object$p_value
SE <- object$SE
r.square <- object$r.square
# [1] 0.4405716 0.5588220 0.6660910 0.5553701

results.cof <- round(do.call("cbind", theta.hat),5)
rownames(results.cof) <- c("Intercept", "Spatial.Lag", "Spatial*Liquidity", "building.age",
 "floor", "Area", "building.age^2", "floor^2", "Area^2" )
results.pvalue <- round(do.call("cbind", p_value),5)
rownames(results.pvalue) <- rownames(results.cof)
results.se <- round(do.call("cbind", SE),5)
rownames(results.se) <- rownames(results.cof)

R2 <- object$R2
SSE <- object$SSE  

plot(y,type = "l", xaxt="n", xlab = "Time", ylab = "natural logarithm of housing price of transacted buildings",
 main = "House Data in Taikoo shing, Hong Kong" )
axis(1, at = c(1,cumsum(pt)), c(seq(2008,2016)))
cl <- ccp.local-1 # 1,3,4
change <- cumsum(c(sum(pt[1:cl[1]]), sum(pt[(1+cl[1]):cl[2]]),
                   sum(pt[(1+cl[2]):cl[3]]) ))
abline(v = change, lty = 1, col = 4, cex = 5 )



######### OLS without common breaks ###############
fit <- lm(y ~ xx[,-1])
summary(fit)

##SSE
sum( lm( y ~ xx[,-1])$res^2 )

R2.pool <- 1-398.7973/sum((y - mean(y) )^2)

################### prediction compared ########################
### prediction by pooled model
# da2016 <- subset(dat, format.Date(Trans.Date, "%Y")=="2016")

# Ynew <- as.matrix(log(da2016$`Price(m)`))
# Xnew <- cbind( da2016$`building age`, da2016$`Area (Gross)`, da2016$floor,
#             da2016$floor^2, da2016$`Area (Gross)`^2, da2016$`building age`^2)
# xnew <- apply(Xnew, 2, scale)


# W2016 <- matrix(0, nrow(da2016), nrow(da2016))
# ### weight matrix
# for (i in 2:nrow(da2016)) {
#   for (j in 1:(i-1)) {
#     W2016[i, j] <-
#       ifelse(
#         da2016$`building name`[i] != da2016$`building name`[j] &
#           da2016$Trans.Date[i] - da2016$Trans.Date[j] <= Tlag,
#         1 / dist_mat[da2016$`building name`[i], da2016$`building name`[j]],
#         ifelse(
#           da2016$`building name`[i] == da2016$`building name`[j] &
#             da2016$Trans.Date[i] - da2016$Trans.Date[j] <= Tlag,
#           1 / min(lowerTriangle(dist_mat)), 0
#           )
#         )
#   }
#   print(i)
# }
# W2016 <- W2016/1000

# Vnum2016 <- function(dd) {
#   v1 <- sapply(dd, function(x, y) {
#     sum(as.numeric(y - x) <= 30 &
#           as.numeric(y - x) >= 1)
#   }, y = da2016$Trans.Date, simplify = T)
  
#   v2 <- sapply(dd, function(x, y) {
#     sum(as.numeric(y - x) <= 60 &
#           as.numeric(y - x) >= 31)
#   }, y = da2016$Trans.Date, simplify = T)
  
#   v3 <- sapply(dd, function(x, y) {
#     sum(as.numeric(y - x) <= 90 &
#           as.numeric(y - x) >= 61)
#   }, y = da2016$Trans.Date, simplify = T)
  
#   V <- v1 + v2/2 + v3/3
#   return(V)
  
# }

# V2016 <-  Vnum2016(da2016$Trans.Date)

# xxnew <- cbind(1, W2016%*%Ynew, V2016*(W2016%*%Ynew), xnew)
xxnew <- xall[grep("2016", datsub$Trans.Date),]
Ynew <- Y[grep("2016", datsub$Trans.Date),]

y.pred.pool <- xxnew %*% fit$coeff
MASE.pool <- mean( abs(Ynew - y.pred.pool)/mean(abs(Ynew - mean(Ynew)) ) )

MSE.pool <- mean((Ynew - y.pred.pool)^2)

MAPE.pool <- mean( abs( (Ynew - y.pred.pool)/Ynew )  )
wMAPE.pool <- sum(abs(Ynew - y.pred.pool))/ sum(abs(Ynew))


### based on the last period
y.pred.cp <- xxnew %*% theta.hat[[4]]
MASE.cp <- mean( abs(Ynew - y.pred.cp)/mean(abs(Ynew - mean(Ynew)) ) )

MSE.cp <- mean( (Ynew - y.pred.cp)^2 )

MAPE.cp <- mean( abs( (Ynew - y.pred.cp)/Ynew)  )
wMAPE.cp <- sum(abs(Ynew - y.pred.cp))/ sum(abs(Ynew))


## ggplot with HK index
library(ggplot2)
load("hkhsi.RData")
ss <- cbind(hkid, datsub)
dd <- ss[,c(1,6)]
colnames(dd) <- colnames(hkid)
d1 <- dd[- grep("2016", datsub$Trans.Date),]
d2 <- dd[grep("2016", datsub$Trans.Date),]
d3 <- data.frame(Date = d2$Date, Index= y.pred.cp)
d4 <- data.frame(Date = d2$Date, Index= y.pred.pool)

ggplot(aes(x = Date), data = dd) + geom_line(aes(y = Index, linetype = "Transaction Price"), data = d1) + 
  geom_line(aes(y = Index, linetype = "Transaction Price in 2016"), data = d2) +
  geom_line(aes(y = Index, linetype = "Predict Price by CP"), data = d3) +
  geom_line(aes(y = Index, linetype = "Predict Price by Pooled"), data = d4) +
  geom_line(aes(y = Index, linetype = "HKHS Index"), data = hkid) + 
  scale_linetype_manual("", values=c("Transaction Price" =1, "Transaction Price in 2016"=2,
   "Predict Price by CP"=4, "Predict Price by Pooled"=9, "HKHS Index" = 3)) +
  xlab("Time") + ylab("Daily Views") + ggtitle("Natural Logarithm of House Data in Taikoo shing, Hong Kong vs. HK HSI") + 
   theme(plot.title = element_text(hjust = 0.5),axis.line = element_line(colour = "black"), panel.border = element_blank(), 
         panel.background = element_blank()) + 
  geom_vline( xintercept = as.numeric(dd$Date[change]), lty = 2 )


ggplot(aes(x = Date), data = dd) + geom_line(aes(y = Index, linetype = "Transaction prices",color="Transaction prices"), data = d1) + 
  geom_line(aes(y = Index, linetype = "Transaction prices in 2016", color="Transaction prices in 2016"), data = d2) +
  geom_line(aes(y = Index, linetype = "Predicted prices by a SLHP model with three change-points", 
                color="Predicted prices by a SLHP model with three change-points"), data = d3) +
  geom_line(aes(y = Index, linetype = "Predicted prices by a SLHP model without any change-point",
                color="Predicted prices by a SLHP model without any change-point"), data = d4) +
  geom_line(aes(y = Index, linetype = "Hong Kong HSI",color="Hong Kong HSI"), data = hkid) + 
  scale_color_manual("",values=c("Transaction prices"=1, "Transaction prices in 2016"=2,
                                 "Predicted prices by a SLHP model with three change-points"=3, 
                                 "Predicted prices by a SLHP model without any change-point"=8, "Hong Kong HSI"=4))+
  scale_linetype_manual("", values=c("Transaction prices" =1, "Transaction prices in 2016"=2,
                                     "Predicted prices by a SLHP model with three change-points"=4, 
                                     "Predicted prices by a SLHP model without any change-point"=5, "Hong Kong HSI" = 3)) +
  xlab("Time") + ylab("") + ggtitle("Natural Logarithm of House Data in Taikoo Shing, Hong Kong vs. Hong Kong HSI") + 
  guides(fill=guide_legend(ncol = 2, byrow = TRUE)) +
  theme(plot.title = element_text(hjust = 0.5),legend.position = "bottom",legend.direction = "vertical" ) + # legend.position="none"
  # theme(plot.title = element_text(hjust = 0.5),axis.line = element_line(colour = "black"), 
  #       panel.background=element_rect(fill="transparent",colour=NA),
  #       plot.background=element_rect(fill="transparent",colour=NA),
  #       legend.key = element_rect(fill = "transparent", colour = "transparent")) +
  geom_vline( xintercept = as.numeric(dd$Date[change]), lty = 2 ) 



