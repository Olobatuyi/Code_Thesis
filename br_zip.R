knitr::opts_chunk$set(echo = FALSE)
library(knitr)
require(keras)
library(tidyverse)
library(OpenMPController)
library(caret)
library(flexmix)
library(lattice)
library(MASS)
require(pscl) # alternatively can use package ZIM for zero-inflated models
library(lmtest)
library(sandwich)
library(ggplot2)
library(SIBER)
library(factoextra)
library(FactoMineR)
library(xtable)
library(clues)
library(MBCbook)
library(HDclassif)
library(dplyr)
library(SIBER)
library(pROC)
require("markdown")
require("rattle")
require("xtable")
require("stringr")
require("fBasics")
require("MASS")
require("survival")
require("STAR")
require("gamlss.dist")
require("VGAM")
require("ggplot2")
library(pscl)
library(sandwich)
library(data.table)
library(plotROC);library(gridExtra);library(grid)
cores <- parallel::detectCores()
#install.packages("optimization")
library(optimization)
#graphics.off()
library(mise)

## Contraceptive use among married women

cmd <- read.table("C:\\Users\\Olobatuyi\\Downloads\\cmcdata.txt", sep = ",")
head(cmd)

Y1 <- cmd$V10
Y <- Y1 - 1
table(Y)
w <- cbind(log(cmd$V2+0.5),log(cmd$V3+0.5), log(cmd$V4+0.5),log(cmd$V7+0.5),
            log(cmd$V8+0.5),log(cmd$V9+0.5))



xe <- prcomp(w)
summary(xe)
w <- as.matrix(xe$x[,1])
u <- cmd$V5; v <- cmd$V6
x <- cbind(w,u,v)

############### Optim SA
zipcwm <- function(Y, w, u, v, k=3, maxit = 1000, tol = 0.1, show_table = FALSE){
  
  x <- cbind(rep(1, nrow(w)), w,u,v);
  
  if(!is.data.frame(x)) x <- unname(as.matrix(x))
  eps <- sqrt(.Machine$double.eps); n <- nrow(x)
  
  W1 <- keras::to_categorical(u); W2 <- keras::to_categorical(v)
  d <- ncol(x); d1 <- ncol(w); c1 <- ncol(W1); c2 <- ncol(W2)
  L1 <- ai <- ai3 <- lb <- L <- NULL
  L2 <- lf1 <- lf2 <- lp <- k1 <- matrix(NA, nrow = n, ncol = (k-1))
  
  dm1 <- lki <- lk2 <- Linf <- al <- NULL
  pr <- matrix(NA, nrow = n, ncol = (k-1))
  lp1 <- dm <- lam <- matrix(NA, nrow = n, ncol = (k-1))
  
  b <- matrix(runif(d*(k-1)), nrow = (k-1))
  sigma <- wew <- list()
  mean <- matrix(runif(d1*(k-1)), nrow = k-1)
  
  P1i <- rep(1/k, k)
  P1 <- P1i[1]; Pi <- P1i[-1]
  
  for(l in 1:(k-1))sigma[[l]] <- diag(d1)
  
  m <- (k-1)*ncol(b) + (2 * (k-1) * d) + (k-1) + (k-1)*(d*d - d)/2
  
  py11 <- matrix(colMeans(W1), nrow = (k-1), ncol = c1, byrow = TRUE)
  py21 <- matrix(colMeans(W2), nrow = (k-1), ncol = c2, byrow = TRUE)
  
  sig1 <- sig <- array(NA, c(d1,d1,n))
  et1 <- list()
  
  count <- 2;
  
  L <- c(-16000,-15000,-14000);
  
  repeat{
    
    a <- x %*% t(b)
    lam <- exp(a)
    
    for (i in 1:n) {
      
      # Zero count
      dm1[i] <- VGAM::dzipois(Y[i], lambda = 0) * P1
      
      for(l in 1:(k-1)){
        
        # Non-zero count
        dm[i,l] <- Pi[l] * mvnfast::dmvn(w[i,], mean[l,], sigma[[l]], ncores = 8) *
          dpois(Y[i], lambda = lam[i,l]) * dmultinom(W1[i,], 1, prob = py11[l,], log = FALSE) *
          dmultinom(W2[i,], 1, prob = py21[l,], log = FALSE)
      }
      
    }
    
    po  <- dm  / rowSums(cbind(dm1, dm))
    po1 <- dm1 / rowSums(cbind(dm1, dm)) #1 - rowSums(po)
    Poc <- cbind(po1, po)
    
    ### For non-zeros
    
    for (l in 1:(k-1)) {
      
      for (i in 1:n) {
        
        lf1[i,l] <- po[i,l] * dmultinom(W1[i,], 1, prob = py11[l,], log = FALSE) 
        lf2[i,l] <- po[i,l] * dmultinom(W2[i,], 1, prob = py21[l,], log = FALSE)
        
      }
      
      lp [,l] <- po[,l ] * dpois(Y, lambda = lam[,l], log = FALSE)
      lp1[,l] <- po[,l]  * Pi[l]
      L2 [,l] <- po[,l]  * mvnfast::dmvn(w, mean[l,], sigma[[l]], ncores = 8,
                                         isChol = TRUE, log = FALSE)
      
      k1[,l] <- lf1[,l] * lf2[,l] * lp[,l] * L2[,l] * lp1[,l]
      
    }
    
    poip <- po1 * P1
    
    ad <- cbind(poip, k1)
    
    for(l in 1:(k-1)){
      
      LogLike <- function(par) {
        
        lambda <- exp(x %*% par)
        LL <- -sum(po[,l] * dpois(Y, lambda, log = TRUE))
        return(LL)
        
      }
      
      par <- b[l,]
      
      m.like <- optim_sa(fun = LogLike, start = par, 
                         lower = rep(eps,ncol(x)),
                         upper = rep(1, ncol(x)), trace = TRUE,
                         control = list(t0 = 100, nlimit = 550, t_min = 0.1,
                                        dyn_rf = F, rf = 1, r = 0.7))
      b[l,] <- m.like$par
      
    }
    
    # M-Step Updates
    P1 <- sum(po1) / n
    
    for(l in 1:(k-1)) mean[l,] <- colSums(po[,l] * w) / colSums(po)[l]
    
    for (l in 1:(k-1)) {
      
      for (i in 1:nrow(x)) {sig[,,i] <- po[i,l] * outer((w[i,] - mean[l,]),(w[i,] - mean[l,]))}
      sigma[[l]] <- apply(sig, c(1, 2), sum)/ sum(po[,l]) + diag(eps, ncol(w))
      
    }
    
    for (l in 1:(k-1)) {
      
      py11[l, ] <- (colSums(po[,l] * W1) / colSums(po)[l])
      py21[l, ] <- (colSums(po[,l] * W2) / colSums(po)[l])
      
    }
    
    #}
    
    Pi <- colSums(po) / n
    
    # Compute Log_likelihood lki, lk2,
    
    L[count] <- sum(log(rowSums(ad)))
    a_k <- (L[count+1] - L[count]) / (L[count] - L[count-1])
    L[count + 2] <- L[count] + ((1-a_k)^-1 * (L[count+1] - L[count]))
    
    if (show_table) {
      
      dif <- abs(L[count+2] - L[count+1])
      out_table = data.frame(Iteration = count, Likelihood = L[count+2], difference = dif)
      print(kable(out_table))
      
      if (count == maxit || dif < tol) break;
      
    }
    
    count <- count + 1
    
  }
  
  Z <- apply(Poc, 1, which.max)
  
  Prof <- Poc1 <- Poc
  
  for(i in 1:n){
    
    for(j in 1:k){
      
      Prof[i,j] <- ifelse(Z[i] == j, 1, 0)
      
      if(Poc[i,j] > 0){
        
        Poc1[i,j] <- log(Poc[i,j])
      }
      
    }
    
  }
  
  
  ai   <- -2*L[count + 2] - 2*m;
  ai3  <- -2*L[count + 2] - m*log(n)  
  ICL  <-  ai3 + suppressWarnings(sum(rowSums(Prof * Poc1)))
  AIcc <-  ai - 2*m*(m+1)/(n-m-1)
  AIC3 <- -2*L[count+2] - 3*m
  AICu <-  AIcc - n*log(n/(n-m-1))
  ICL  <-  ai3 + suppressWarnings(sum(rowSums(Prof * Poc1)))
  Caic <- -2*L[count+2] - m*(1+log(n))
  AWE  <- -2*L[count+2] - 2*m*(3/2 + log(n))
  
  return(list("mean" = mean, "Prob1" = P1, "prob2" = Pi, "post" = Poc, "poiwei" = b,
              "classification" = Z, "logLik" = L, "AIC" = ai, "BIC" = ai3,
              "sigma" = sigma, "ICL" = ICL, "AICc" = AIcc, "AIC3" = AIC3, 
              "AICu" = AICu, "Caic" = Caic, "AWE" = AWE, "py11" = py11, "py21" = py21))
  
}

cmw_tru_use_real3 <- function(Y, w, u, v, k, maxit = 1000, tol = 0.1, show_table = FALSE){
  
  x <- cbind(rep(1, nrow(w)), w,u,v);
  
  if(!is.data.frame(x)) x <- unname(as.matrix(x))
  eps <- sqrt(.Machine$double.eps); n <- nrow(x)
  
  W1 <- keras::to_categorical(u); W2 <- keras::to_categorical(v)
  d <- ncol(x); d1 <- ncol(w); c1 <- ncol(W1); c2 <- ncol(W2)
  L1 <- ai <- ai3 <- lb <- L <- NULL
  L2 <- lf1 <- lf2 <- lp <- k1 <- matrix(NA, nrow = n, ncol = (k-1))
  
  dm1 <- lki <- lk2 <- Linf <- al <- NULL
  pr <- matrix(NA, nrow = n, ncol = (k-1))
  lp1 <- dm <- lam <- matrix(NA, nrow = n, ncol = (k-1))
  
  b <- matrix(runif(d*(k-1)), nrow = (k-1))
  sigma <- wew <- list()
  mean <- matrix(runif(d1*(k-1)), nrow = k-1)
  
  P1i <- rep(1/k, k)
  P1 <- P1i[1]; Pi <- P1i[-1]
  
  for(l in 1:(k-1))sigma[[l]] <- diag(d1)
  
  m <- (k-1)*ncol(b) + (2 * (k-1) * d) + (k-1) + (k-1)*(d*d - d)/2
  
  py11 <- matrix(colMeans(W1), nrow = (k-1), ncol = c1, byrow = TRUE)
  py21 <- matrix(colMeans(W2), nrow = (k-1), ncol = c2, byrow = TRUE)
  
  sig1 <- sig <- array(NA, c(d1,d1,n))
  et1 <- list()
  
  count <- 2;
  
  L <- c(-16000,-15000,-14000);
  
  repeat{
    
    a <- x %*% t(b)
    lam <- exp(a)
    
    for (i in 1:n) {
      
      # Zero count
      dm1[i] <- VGAM::dzipois(Y[i], lambda = 0) * P1
      
      for(l in 1:(k-1)){
        
        # Non-zero count
        dm[i,l] <- Pi[l] * mvnfast::dmvn(w[i,], mean[l,], sigma[[l]], ncores = 8) *
          dpois(Y[i], lambda = lam[i,l]) * dmultinom(W1[i,], 1, prob = py11[l,], log = FALSE) *
          dmultinom(W2[i,], 1, prob = py21[l,], log = FALSE)
      }
      
    }
    
    po  <- dm  / rowSums(cbind(dm1, dm))
    po1 <- dm1 / rowSums(cbind(dm1, dm)) #1 - rowSums(po)
    Poc <- cbind(po1, po)
    
    Prof <- Poc1 <- Poc
    
    Z <- apply(Poc, 1, which.max)
    
    for (i in 1:n) {
      
      for (j in 1:k) {
        
        Prof[i,j] <- ifelse(Z[i] == j, 1, 0)
        
      }
    }
    
    po <- as.matrix(Prof[,-1])
    
    ### For non-zeros
    
    for (l in 1:(k-1)) {
      
      for (i in 1:n) {
        
        lf1[i,l] <- po[i,l] * dmultinom(W1[i,], 1, prob = py11[l,], log = FALSE) 
        lf2[i,l] <- po[i,l] * dmultinom(W2[i,], 1, prob = py21[l,], log = FALSE)
        
      }
      
      lp [,l] <- po[,l ] * dpois(Y, lambda = lam[,l], log = FALSE)
      lp1[,l] <- po[,l]  * Pi[l]
      L2 [,l] <- po[,l]  * mvnfast::dmvn(w, mean[l,], sigma[[l]], ncores = 8,
                                         isChol = TRUE, log = FALSE)
      
      k1[,l] <- lf1[,l] * lf2[,l] * lp[,l] * L2[,l] * lp1[,l]
      
    }
    
    poip <- Prof[,1] * P1
    
    ad <- cbind(poip, k1)
    
    for(l in 1:(k-1)){
      
      LogLike <- function(par) {
        
        lambda <- exp(x %*% par)
        LL <- -sum(po[,l] * dpois(Y, lambda, log = TRUE))
        return(LL)
        
      }
      
      par <- b[l,]
      
      m.like <- optim_sa(fun = LogLike, start = par, 
                         lower = rep(eps,ncol(x)),
                         upper = rep(1, ncol(x)), trace = TRUE,
                         control = list(t0 = 100, nlimit = 550, t_min = 0.1,
                                        dyn_rf = F, rf = 1, r = 0.7))
      b[l,] <- m.like$par
      
    }
    
    # M-Step Updates
    P1 <- sum(Prof[,1]) / n
    
    for(l in 1:(k-1)) mean[l,] <- colSums(po[,l] * w) / colSums(po)[l]
    
    for (l in 1:(k-1)) {
      
      for (i in 1:nrow(x)) {sig[,,i] <- po[i,l] * outer((w[i,] - mean[l,]),(w[i,] - mean[l,]))}
      sigma[[l]] <- apply(sig, c(1, 2), sum)/ sum(po[,l]) + diag(eps, ncol(w))
      
    }
    
    for (l in 1:(k-1)) {
      
      py11[l, ] <- (colSums(po[,l] * W1) / colSums(po)[l])
      py21[l, ] <- (colSums(po[,l] * W2) / colSums(po)[l])
      
    }
    
    #}
    
    Pi <- colSums(po) / n
    
    # Compute Log_likelihood lki, lk2,
    
    L[count] <- sum(log(rowSums(ad)))
    # a_k <- (L[count+1] - L[count]) / (L[count] - L[count-1])
    # L[count + 2] <- L[count] + ((1-a_k)^-1 * (L[count+1] - L[count]))
    
    if (show_table) {
      
      dif <- abs(L[count] - L[count-1])
      out_table = data.frame(Iteration = count, Likelihood = L[count], difference = dif)
      print(kable(out_table))
      
      if (count == maxit || dif < tol) break;
      
    }
    
    count <- count + 1
    
  }
  
  # Z <- apply(Poc, 1, which.max)
  
  Poc1 <- ifelse(Poc > 0, log(Poc), 0)
  
  ai   <- -2*L[count] - 2*m;
  ai3  <- -2*L[count] - m*log(n)  
  ICL  <-  ai3 + suppressWarnings(sum(rowSums(Prof * Poc1)))
  AIcc <-  ai - 2*m*(m+1)/(n-m-1)
  AIC3 <- -2*L[count] - 3*m
  AICu <-  AIcc - n*log(n/(n-m-1))
  ICL  <-  ai3 + suppressWarnings(sum(rowSums(Prof * Poc1)))
  Caic <- -2*L[count] - m*(1+log(n))
  AWE  <- -2*L[count] - 2*m*(3/2 + log(n))
  
  return(list("mean" = mean, "Prob1" = P1, "prob2" = Pi, "post" = Poc, "poiwei" = b,
              "classification" = Z, "logLik" = L, "AIC" = ai, "BIC" = ai3,
              "sigma" = sigma, "ICL" = ICL, "AICc" = AIcc, "AIC3" = AIC3, 
              "AICu" = AICu, "Caic" = Caic, "AWE" = AWE))
  
}

## Recorded results
cw <- zipcwm(Y, w, u, v, k = 3, maxit = 1000, tol = 0.05, show_table = TRUE)


FIn11 <- cw

cm <- confusionMatrix(sort(as.factor(FIn11$classification-1)), sort(as.factor(Y)))

table(FIn11$classification)

FIn11$Prob1

FIn11$BIC


multiclass.roc(ordered(Y),ordered(FIn11$classification-1))
multiclass.roc(ordered(z),ordered(FIn2$classification))
multiclass.roc(ordered(z),ordered(FIn4$classification))
multiclass.roc(ordered(z),ordered(FIn5$classification))

adjustedRand(Y,FIn11$classification-1)

library(cluster)
library(fpc)
library(tsne)
#clusplot(cmd[-10], sort(as.factor(Fe3$classification)), color = T, shade = T, labels = 2, lines = 0)

###plotly
library(parallel)
library(OpenMPController)
library(data.table)
library(tsne)
library(Rtsne)
library(factoextra)
library(FactoMineR)

tsne_out <- Rtsne(cmd[,-10], pca = F, theta=0.0, check_duplicates = FALSE)

#r <- PCA(cmd[,-10], ncp = 5, graph = FALSE)
r <- PCA(tsne_out$Y, ncp = 5, graph = FALSE)

fviz_pca_ind(r, geom = "point", axes = c(1,2), 
             habillage = (as.factor(FIn11$classification-1)),
             palette = c("red", "blue", "gold","magenta", "black", "cyan","orange", 
                         "grey", "green","brown"),
             addEllipse = T, repel = F,
             ggtheme = theme_minimal(), title = "")

adjustedRand((Y),(FIn11$classification-1))
adjustedRand((z),(FIn31$classification))
adjustedRand((z),(FIn41$classification))
adjustedRand((z),(FIn51$classification))

# Confusion Matrix 
cm <- confusionMatrix(sort(as.factor(FIn11$classification-1)), sort(as.factor(Y)))

conf <- function(x){
  
  cm_d  <- as.data.frame(x$table)
  cm_st <- data.frame(x$overall)
  cm_st$x.overall <- round(cm_st$x.overall, 2)
  cm_p <- as.data.frame(prop.table(x$table))
  cm_d$perc <- round(cm_p$Freq*100,2)
  cm_d_p <- ggplot(data = cm_d, aes(x = Prediction, y = Reference, fill = Freq)) +
    geom_tile()+
    geom_text(aes(label = paste("", Freq, ",",perc, "%")), col = "red")+
    theme_light()+
    guides(fill = FALSE)
  cm_st_p <- tableGrob(cm_st)
  
  print(grid.arrange(cm_d_p, cm_st_p, nrow = 1, ncol = 2,
                     top = textGrob("Confusion Matrix and Statistics", gp = gpar(fontsize = 25, font = 1))))
}

conf(cm)

####### FZIP
fZip <- function(Y, w, u, v, k, maxit = 1000, tol = 0.1, show_table = FALSE){
  
  x <- cbind(rep(1, nrow(w)), w,u,v);
  
  if(!is.data.frame(x)) x <- unname(as.matrix(x))
  eps <- sqrt(.Machine$double.eps); n <- nrow(x)

  d <- ncol(x);
  L1 <- ai <- ai3 <- lb <- L <- NULL
  L2 <- lf1 <- lf2 <- lp <- k1 <- matrix(NA, nrow = n, ncol = (k-1))
  
  dm1 <- lki <- lk2 <- Linf <- al <- NULL
  pr <- matrix(NA, nrow = n, ncol = (k-1))
  lp1 <- dm <- lam <- matrix(NA, nrow = n, ncol = (k-1))
  
  b <- matrix(runif(d*(k-1)), nrow = (k-1))
  
  P1i <- rep(1/k, k)
  P1 <- P1i[1]; Pi <- P1i[-1]
  
  m <- (k-1)*ncol(b) + (k-1)
  
  et1 <- list()
  
  count <- 2;
  
  L <- c(-16000,-15000,-14000);
  
  repeat{
    
    a <- x %*% t(b)
    lam <- exp(a)
    
    for (i in 1:n) {
      
      # Zero count
      dm1[i] <- VGAM::dzipois(Y[i], lambda = 0) * P1
      
      for(l in 1:(k-1)){
        
        # Non-zero count
        dm[i,l] <- Pi[l] * dpois(Y[i], lambda = lam[i,l])
      }
      
    }
    
    po  <- dm  / rowSums(cbind(dm1, dm))
    po1 <- dm1 / rowSums(cbind(dm1, dm))
    Poc <- cbind(po1, po)
    
    ### For non-zeros
    for (l in 1:(k-1)) {
      
      lp [,l] <- po[,l ] * dpois(Y, lambda = lam[,l], log = FALSE)
      lp1[,l] <- po[,l]  * Pi[l]

      k1[,l] <- lp[,l] * lp1[,l]
      
    }
    
    poip <- po1 * P1
    
    ad <- cbind(poip, k1)
    
    for(l in 1:(k-1)){
      
      LogLike <- function(par) {
        
        lambda <- exp(x %*% par)
        LL <- -sum(po[,l] * dpois(Y, lambda, log = TRUE))
        return(LL)
        
      }
      
      par <- b[l,]
      
      m.like <- optim_sa(fun = LogLike, start = par, 
                         lower = rep(eps,ncol(x)),
                         upper = rep(1, ncol(x)), trace = TRUE,
                         control = list(t0 = 100, nlimit = 550, t_min = 0.1,
                                        dyn_rf = F, rf = 1, r = 0.7))
      b[l,] <- m.like$par
      
    }
    
    # M-Step Updates
    P1 <- sum(po1) / n; Pi <- colSums(po) / n
    
    # Compute Log_likelihood lki, lk2,
    
    L[count] <- sum(log(rowSums(ad)))
    a_k <- (L[count+1] - L[count]) / (L[count] - L[count-1])
    L[count + 2] <- L[count] + ((1-a_k)^-1 * (L[count+1] - L[count]))
    
    if (show_table) {
      
      dif <- abs(L[count+2] - L[count+1])
      out_table = data.frame(Iteration = count, Likelihood = L[count+2], difference = dif)
      print(kable(out_table))
      
      if (count == maxit || dif < tol) break;
      
    }
    
    count <- count + 1
    
  }
  
  Z <- apply(Poc, 1, which.max)
  
  Prof <- Poc1 <- Poc
  
  for(i in 1:n){
    
    for(j in 1:k){
      
      Prof[i,j] <- ifelse(Z[i] == j, 1, 0)
      
      if(Poc[i,j] > 0){
        
        Poc1[i,j] <- log(Poc[i,j])
      }
      
    }
    
  }
  
  
  ai   <- -2*L[count + 2] - 2*m;
  ai3  <- -2*L[count + 2] - m*log(n)  
  ICL  <-  ai3 + suppressWarnings(sum(rowSums(Prof * Poc1)))
  AIcc <-  ai - 2*m*(m+1)/(n-m-1)
  AIC3 <- -2*L[count+2] - 3*m
  AICu <-  AIcc - n*log(n/(n-m-1))
  ICL  <-  ai3 + suppressWarnings(sum(rowSums(Prof * Poc1)))
  Caic <- -2*L[count+2] - m*(1+log(n))
  AWE  <- -2*L[count+2] - 2*m*(3/2 + log(n))
  
  return(list("Prob1" = P1, "prob2" = Pi, "post" = Poc, "poiwei" = b,
              "classification" = Z, "logLik" = L, "AIC" = ai, "BIC" = ai3,
              "ICL" = ICL, "AICc" = AIcc, "AIC3" = AIC3, 
              "AICu" = AICu, "Caic" = Caic, "AWE" = AWE))
  
}

fz <- fZip(Y, w, u, v, k = 3, maxit = 1000, tol = 0.05, show_table = TRUE)

cmfz <- confusionMatrix((as.factor(fz$classification-1)), (as.factor(Y)))

table(fz$classification)

### Poisson CWM
pcwm <- function(Y, w, u, v, k=3, maxit = 1000, tol = 0.1, show_table = FALSE){
  
  x <- cbind(rep(1, nrow(w)), w,u,v);
  
  if(!is.data.frame(x)) x <- unname(as.matrix(x))
  eps <- sqrt(.Machine$double.eps); n <- nrow(x)
  
  W1 <- keras::to_categorical(u); W2 <- keras::to_categorical(v)
  d <- ncol(x); d1 <- ncol(w); c1 <- ncol(W1); c2 <- ncol(W2)
  L1 <- ai <- ai3 <- lb <- L <- NULL
  L2 <- lf1 <- lf2 <- lp <- k1 <- matrix(NA, nrow = n, ncol = k)
  
  dm1 <- lki <- lk2 <- Linf <- al <- NULL
  pr <- matrix(NA, nrow = n, ncol = k)
  lp1 <- dm <- lam <- matrix(NA, nrow = n, ncol = k)
  
  b <- matrix(runif(d*k), nrow = k)
  sigma <- wew <- list()
  mean <- matrix(runif(d1*k), nrow = k)
  
  P1i <- rep(1/k, k)
  #P1 <- P1i[1]; Pi <- P1i[-1]
  
  for(l in 1:k)sigma[[l]] <- diag(d1)
  
  m <- k*ncol(b) + (2 * k * d) + k + k*(d*d - d)/2
  
  py11 <- matrix(colMeans(W1), nrow = (k), ncol = c1, byrow = TRUE)
  py21 <- matrix(colMeans(W2), nrow = (k), ncol = c2, byrow = TRUE)
  
  sig1 <- sig <- array(NA, c(d1,d1,n))
  et1 <- list()
  
  count <- 2;
  
  L <- c(-16000,-15000,-14000);
  
  repeat{
    
    a <- x %*% t(b)
    lam <- exp(a)
    
    for (i in 1:n) {
      
      for(l in 1:k){
        
        # Non-zero count
        dm[i,l] <- P1i[l] * mvnfast::dmvn(w[i,], mean[l,], sigma[[l]], ncores = 8) *
          dpois(Y[i], lambda = lam[i,l]) * dmultinom(W1[i,], 1, prob = py11[l,], log = FALSE) *
          dmultinom(W2[i,], 1, prob = py21[l,], log = FALSE)
      }
      
    }
    
    po  <- dm  / rowSums(dm)
    Poc <- po
    
    ### For non-zeros
    
    for (l in 1:k) {
      
      for (i in 1:n) {
        
        lf1[i,l] <- po[i,l] * dmultinom(W1[i,], 1, prob = py11[l,], log = FALSE) 
        lf2[i,l] <- po[i,l] * dmultinom(W2[i,], 1, prob = py21[l,], log = FALSE)
        
      }
      
      lp [,l] <- po[,l ] * dpois(Y, lambda = lam[,l], log = FALSE)
      lp1[,l] <- po[,l]  * P1i[l]
      L2 [,l] <- po[,l]  * mvnfast::dmvn(w, mean[l,], sigma[[l]], ncores = 8,
                                         isChol = TRUE, log = FALSE)
      
      k1[,l] <- lf1[,l] * lf2[,l] * lp[,l] * L2[,l] * lp1[,l]
      
    }
    
    ad <- k1
    
    for(l in 1:k){
      
      LogLike <- function(par) {
        
        lambda <- exp(x %*% par)
        LL <- -sum(po[,l] * dpois(Y, lambda, log = TRUE))
        return(LL)
        
      }
      
      par <- b[l,]
      
      m.like <- optim_sa(fun = LogLike, start = par, 
                         lower = rep(eps,ncol(x)),
                         upper = rep(1, ncol(x)), trace = TRUE,
                         control = list(t0 = 100, nlimit = 550, t_min = 0.1,
                                        dyn_rf = F, rf = 1, r = 0.7))
      b[l,] <- m.like$par
      
    }
    
    # M-Step Updates
    for(l in 1:k) mean[l,] <- colSums(po[,l] * w) / colSums(po)[l]
    
    for (l in 1:k) {
      
      for (i in 1:nrow(x)) {sig[,,i] <- po[i,l] * outer((w[i,] - mean[l,]),(w[i,] - mean[l,]))}
      sigma[[l]] <- apply(sig, c(1, 2), sum)/ sum(po[,l]) + diag(eps, ncol(w))
      
    }
    
    for (l in 1:k) {
      
      py11[l, ] <- (colSums(po[,l] * W1) / colSums(po)[l])
      py21[l, ] <- (colSums(po[,l] * W2) / colSums(po)[l])
      
    }
    
    #}
    
    P1i <- colSums(po) / n
    
    # Compute Log_likelihood lki, lk2,
    
    L[count] <- sum(log(rowSums(ad)))
    a_k <- (L[count+1] - L[count]) / (L[count] - L[count-1])
    L[count + 2] <- L[count] + ((1-a_k)^-1 * (L[count+1] - L[count]))
    
    if (show_table) {
      
      dif <- abs(L[count+2] - L[count+1])
      out_table = data.frame(Iteration = count, Likelihood = L[count+2], difference = dif)
      print(kable(out_table))
      
      if (count == maxit || dif < tol) break;
      
    }
    
    count <- count + 1
    
  }
  
  Z <- apply(Poc, 1, which.max)
  
  Prof <- Poc1 <- Poc
  
  for(i in 1:n){
    
    for(j in 1:k){
      
      Prof[i,j] <- ifelse(Z[i] == j, 1, 0)
      
      if(Poc[i,j] > 0){
        
        Poc1[i,j] <- log(Poc[i,j])
      }
      
    }
    
  }
  
  
  ai   <- -2*L[count + 2] - 2*m;
  ai3  <- -2*L[count + 2] - m*log(n)  
  ICL  <-  ai3 + suppressWarnings(sum(rowSums(Prof * Poc1)))
  AIcc <-  ai - 2*m*(m+1)/(n-m-1)
  AIC3 <- -2*L[count+2] - 3*m
  AICu <-  AIcc - n*log(n/(n-m-1))
  ICL  <-  ai3 + suppressWarnings(sum(rowSums(Prof * Poc1)))
  Caic <- -2*L[count+2] - m*(1+log(n))
  AWE  <- -2*L[count+2] - 2*m*(3/2 + log(n))
  
  return(list("mean" = mean, "prob" = P1i, "post" = Poc, "poiwei" = b,
              "classification" = Z, "logLik" = L, "AIC" = ai, "BIC" = ai3,
              "sigma" = sigma, "ICL" = ICL, "AICc" = AIcc, "AIC3" = AIC3, 
              "AICu" = AICu, "Caic" = Caic, "AWE" = AWE, "py11" = py11, "py21" = py21))
  
}

pc <- pcwm(Y, w, u, v, k = 3, maxit = 1000, tol = 0.05, show_table = TRUE)

pc$mean
pc$sigma
pc$poiwei
pc$AIC
FIn11$AIC

cmp <- confusionMatrix(sort(as.factor(pc$classification-1)), sort(as.factor(Y)))
table(pc$classification)
conf(cmp)

we <- glm(Y ~ x, family = poisson)


summary(we)

#### Dispersion Test

AER::dispersiontest(we)
