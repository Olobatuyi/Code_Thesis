<<<<<<< HEAD:MultinomialCWM.R
library(knitr)
library(gridExtra);library(grid);require(keras);library(tidyverse);
library(OpenMPController);library(caret);library(ggplot2);library(SIBER)
library(factoextra);library(FactoMineR);library(xtable);library(clues)
library(MBCbook);library(HDclassif);library(dplyr);library(SIBER)
require("markdown");require("rattle");require("xtable");require("stringr")
require("fBasics");require("MASS");require("survival");require("STAR")
require("gamlss.dist");require("VGAM");library(rgl);library(ellipse)
library(threejs);library(plotROC);library(ROCR);library(pROC)


# Set the number of column of X and fix G

d <- 2; k <- 2; c <- 3-1; n <- 200; ni <- 1
df <- list(); w <- sigma <- list()

Pi <- c(0.6, 0.4); mean <- matrix(c(.1,-2,2,1), nrow = k)

for(l in 1:k) {

  sigma[[l]] <- diag(d)

}

weights1 <- matrix(c(.2,.2,.2,1,1,1,1,1), nrow = d, ncol = c)
weights2 <- matrix(c(1,1,1,1,1,1,1,1), nrow = d, ncol = c)
w[[1]] <- t(weights1); w[[2]] <- t(weights2)

z <- NULL
# Generating Algorithm
x <- matrix(NA, nrow = n, ncol = d)
a <- pr <- matrix(NA, nrow = n, ncol = c)
Y <- py <- matrix(NA, nrow = n, ncol = c+1)
set.seed(20); U <- runif(n)

for (i in 1:n) {

  if(U[i] <= Pi[1]){

    z[i] <- 1
    x[i,] <- mvnfast::rmvn(1, mean[1,], sigma[[1]])
    a[i,]    <- as.matrix(x[i,] %*% t(w[[1]]))
    pr[i,]   <- exp(a[i,]) / (1 + sum(exp(a[i,])))
    pu       <- 1 / (1 + sum(exp(a[i,])))
    py[i,]   <- c(pu, pr[i,])
    Y[i,] <- t(rmultinom(1, ni, py[i,]))
    #ni <- ni + i

  }else{

    z[i] <- 2
    x[i,] <- mvnfast::rmvn(1, mean[2,], sigma[[2]])
    a[i,]    <- as.matrix(x[i,] %*% t(w[[2]]))
    pr[i,]   <- exp(a[i,]) / (1 + sum(exp(a[i,])))
    pu       <- 1 / (1 + sum(exp(a[i,])))
    py[i,]   <- c(pu, pr[i,])
    Y[i,] <- t(rmultinom(1, ni, py[i,]))
    #ni <- ni + i
  }

}


da <- data.frame(Y, x)
names(da) <- c("a", "b", "c", "x1", 'x2')
plot(da, col = z)

pos1 <- which(z == 1)
pos  <- which(z == 2)

Y1 <- Y[pos1,];Y2 <- Y[pos,]

colSums(Y1);colSums(Y2)
table(z)

# The number of a, b, c in group 1 is 10,44, 3
# The number of a, b, c in group 2 is 8,10, 25


##################### K = 2

# Set the number of column of X and fix G
d <- 2; k <- 2; ni <- 1; c <- 3-1; n <- 1000

# create an empty list
w <- sigma <- list()
z <- NULL
x <- matrix(NA, nrow = n, ncol = d)
a <- pr <- matrix(NA, nrow = n, ncol = c)
Y <- py <- matrix(NA, nrow = n, ncol = c+1)

# Z parameter
Pi <- rep(1/k,k)
#Pi <- c(5,3,2)
# X parameter
mean <- matrix(c(.1,-2,2,0,1,3), nrow = k, ncol = d)
for(l in 1:k) {sigma[[l]] <- diag(d)}

# Y parameter
# weights1 <- matrix(c(.1,.2,.3,.1), nrow = d, ncol = c);w[[1]] <- t(weights1)
# weights2 <- matrix(c(1,1,1,1), nrow = d, ncol = c);w[[2]] <- t(weights2)
# weights3 <- matrix(c(.1,.1,.1,.1), nrow = d, ncol = c);w[[3]] <- t(weights3)

weights1 <- matrix(c(5,.4,.3,0.04), nrow = d, ncol = c);    w[[1]] <- t(weights1)
weights2 <- matrix(c(.01,.02,2,1), nrow = d, ncol = c);     w[[2]] <- t(weights2)
weights3 <- matrix(c(1, 0.03,0.06,.02), nrow = d, ncol = c);w[[3]] <- t(weights3)

# Set a seed
set.seed(20);U <- runif(n)

for (i in 1:n) {
  
  if(U[i] <= Pi[1]){
    
    z[i] <- 1
    x[i,] <- mvnfast::rmvn(1, mean[1,], sigma[[1]])
    a[i,]    <- as.matrix(x[i,] %*% t(w[[1]]))
    pr[i,]   <- exp(a[i,]) / (1 + sum(exp(a[i,])))
    pu       <- 1 / (1 + sum(exp(a[i,])))
    py[i,]   <- c(pu, pr[i,])
    Y[i,] <- t(rmultinom(1, ni, py[i,]))
    
  }else{
      
      z[i]   <- 2
      x[i,]  <- mvnfast::rmvn(1, mean[2,], sigma[[2]])
      a[i,]  <- as.matrix(x[i,] %*% t(w[[2]]))
      pr[i,] <- exp(a[i,]) / (1 + sum(exp(a[i,])))
      pu     <- 1 / (1 + sum(exp(a[i,])))
      py[i,] <- c(pu, pr[i,])
      Y[i,]  <- t(rmultinom(1, ni, py[i,]))
      
    }
  }

da <- data.frame(Y, x)
names(da) <- c("a", "b", "c", "x1", 'x2')
plot(da, col = z)

pos1 <- which(z == 1)
pos2 <- which(z == 2)

Y1 <- Y[pos1,];Y2 <- Y[pos2,];#Y3 <- Y[pos3,]
colSums(Y1);colSums(Y2)#;colSums(Y3)
table(z)


##################### K = 3

# Set the number of column of X and fix G
d <- 2; k <- 3; ni <- 1; c <- 3-1; n <- 1000

# create an empty list
w <- sigma <- list()
z <- NULL
x <- matrix(NA, nrow = n, ncol = d)
a <- pr <- matrix(NA, nrow = n, ncol = c)
Y <- py <- matrix(NA, nrow = n, ncol = c+1)

# Z parameter
Pi <- rep(1/k,k)
#Pi <- c(5,3,2)
# X parameter
mean <- matrix(c(.1,-2,2,0,1,3), nrow = k, ncol = d)
for(l in 1:k) {sigma[[l]] <- diag(d)}

# Y parameter
# weights1 <- matrix(c(.1,.2,.3,.1), nrow = d, ncol = c);w[[1]] <- t(weights1)
# weights2 <- matrix(c(1,1,1,1), nrow = d, ncol = c);w[[2]] <- t(weights2)
# weights3 <- matrix(c(.1,.1,.1,.1), nrow = d, ncol = c);w[[3]] <- t(weights3)

weights1 <- matrix(c(5, .4, .3, 0.04), nrow = d, ncol = c);     w[[1]] <- t(weights1)
weights2 <- matrix(c(.01, .02, 2,  1), nrow = d, ncol = c);      w[[2]] <- t(weights2)
weights3 <- matrix(c(1,0.03,0.06,.02), nrow = d, ncol = c); w[[3]] <- t(weights3)

# Set a seed
set.seed(20);U <- runif(n)

for (i in 1:n) {

  if(U[i] <= Pi[1]){

    z[i] <- 1
    x[i,] <- mvnfast::rmvn(1, mean[1,], sigma[[1]])
    a[i,]    <- as.matrix(x[i,] %*% t(w[[1]]))
    pr[i,]   <- exp(a[i,]) / (1 + sum(exp(a[i,])))
    pu       <- 1 / (1 + sum(exp(a[i,])))
    py[i,]   <- c(pu, pr[i,])
    Y[i,] <- t(rmultinom(1, ni, py[i,]))
    
  }else

    if(U[i] > Pi[1] & U[i] < (Pi[2] + Pi[1])){

      z[i]   <- 2
      x[i,]  <- mvnfast::rmvn(1, mean[2,], sigma[[2]])
      a[i,]  <- as.matrix(x[i,] %*% t(w[[2]]))
      pr[i,] <- exp(a[i,]) / (1 + sum(exp(a[i,])))
      pu     <- 1 / (1 + sum(exp(a[i,])))
      py[i,] <- c(pu, pr[i,])
      Y[i,]  <- t(rmultinom(1, ni, py[i,]))
      
    }else{

      z[i] <- 3
      x[i,] <- mvnfast::rmvn(1, mean[3,], sigma[[3]])
      a[i,]    <- as.matrix(x[i,] %*% t(w[[3]]))
      pr[i,]   <- exp(a[i,]) / (1 + sum(exp(a[i,])))
      pu       <- 1 / (1 + sum(exp(a[i,])))
      py[i,]   <- c(pu, pr[i,])
      Y[i,] <- t(rmultinom(1, ni, py[i,]))

    }

}

da <- data.frame(Y, x)
names(da) <- c("a", "b", "c", "x1", 'x2')
plot(da, col = z)

pos1 <- which(z == 1)
pos2 <- which(z == 2)
pos3 <- which(z == 3)

Y1 <- Y[pos1,];Y2 <- Y[pos2,];Y3 <- Y[pos3,]
colSums(Y1);colSums(Y2);colSums(Y3)
table(z)



cmw_ran <- function(Y, x, k, maxit = 1000, tol = 1e-10, show_table = FALSE){
  
  if(!is.data.frame(x)) x <- unname(as.matrix(x))
  if(!is.data.frame(Y)) Y <- unname(as.matrix(Y))
  
  eps = sqrt(.Machine$double.eps)
  d <- ncol(x); c <- ncol(Y) - 1
  L <- L1 <- L2 <- ai <- ai3 <- lb <- NULL
  
  n <- nrow(x)
  a <- pr <- matrix(NA, nrow = n, ncol = c)
  py <- matrix(NA, nrow = n, ncol = c+1)
  dm <- matrix(NA, nrow = n, ncol = k)
  
  pu <- NULL
  sig <- array(NA, c(d,d,n))
  et1 <- list()
  
  w <- sigma <- list()
  m <- (k * c * (1 + d)) + (k * ((d*d - d)/2 + 2*d + 1))
  
  Pi <- rep(1/k, k)
  
  mean <- matrix(rnorm(d*k), nrow = k)
  
  for(l in 1:k){
    
    sigma[[l]] <- diag(d)
    weights <- matrix(runif(d*c), nrow = c, ncol = d)
    w[[l]] <- weights
    
  }
  
  count <- 2; 
  
  L <- c(-16000,-15000,-14000);
  
  repeat{
    
    for(l in 1:k){
      
      for (i in 1:nrow(x)) {
        
        a[i,]   <- x[i,] %*% t(w[[l]])
        pr[i,]  <- exp(a[i,]) / (1 + sum(exp(a[i,])))
        pu[i]   <- 1 / (1 + sum(exp(a[i,]))); py[i,]  <- c(pu[i], pr[i,])
        dm[i,l] <- dmultinom(Y[i,], sum(Y[i,]), prob = py[i,]) *
          mvnfast::dmvn(x[i,], mean[l,], sigma[[l]], ncores = 8, 
                        isChol = F) * Pi[l]
        
      }
      
    }
    
    po <- dm / rowSums(dm)
    
    # M-Step
    
    # Mean for the group
    
    for(l in 1:k) mean[l, ] <- colSums(po[,l] * x) / colSums(po)[l]
    
    # Update Prior
    Pi <- colSums(po) / nrow(x)
    
    for (l in 1:k) {
      
      L2[l]  <- sum(po[,l] * mvnfast::dmvn(x, mean[l,], sigma[[l]], ncores = 8, 
                                           isChol = F, log = TRUE))
    }
    
    if(count > 10){
      
      for(l in 1:k){
        
        et <- nnet::multinom(Y ~ x, weights = po[,l], trace = FALSE)
        w[[l]] <- summary(et)$coefficients[,-1]
        lb[l] <- logLik(et)
        et1[[l]] <- summary(et)
        
      }
      
    }
    
    #if(count > 20){
      
      for (l in 1:k) {
        
        for (i in 1:nrow(x)) {
          
          sig[,,i] <- po[i,l] * outer((x[i,] - mean[l,]),(x[i,] - mean[l,]))
          
        }
        
        sigma[[l]] <- apply(sig, c(1, 2), sum)/ sum(po[,l]) #+ diag(eps, ncol(x))
        
        which(is.nan(sigma[[l]]))
        
        sigma[[l]][which(is.nan(sigma[[l]]))] <- 1
      }
      
    #}
    
    # Compute Log_likelihood
    lik  <- sum(log(Pi) * po);lik2 <- sum(L2);lik3 <- sum(lb)
    
    # Compute Log_likelihood lki, lk2, 
    L[count] <- lik + lik2 + lik3
    
    a_k <- (L[count+1] - L[count]) / (L[count] - L[count-1])
    
    L[count + 2] <- L[count] + ((1-a_k)^-1 * (L[count+1] - L[count]))
    
    if (show_table) {
      
      dif <- abs(L[count+2] - L[count+1])
      
      out_table = data.frame(Iteration = count, Likelihood = L[count+2], difference = dif)
      print(kable(out_table))  
      
      if (dif < tol || count == maxit) break; 
      
    }
    
    count <- count + 1
    
  }
  
  Z <- apply(po, 1, which.max)
  
  Prof <- Poc1 <- po
  
  for(i in 1:n){
    
    for(j in 1:k){
      
      Prof[i,j] <- ifelse(po[i,j] > 0.9, 1, 0)
      
      if(po[i,j] > 0){
        
        Poc1[i,j] <- log(po[i,j]) 
        
      }
      
    }
    
  }
  
  ### Information criterion
  
  ai   <- -2*L[count + 2] + 2*m;
  ai3  <- -2*L[count + 2] + m*log(n)  
  AIcc <-  ai - 2*m*(m+1)/(n-m-1)
  AIC3 <- -2*L[count+2] - 3*m
  AICu <-  AIcc - n*log(n/(n-m-1))
  ICL  <-  ai3 + suppressWarnings(sum(rowSums(Prof * Poc1)))
  Caic <- -2*L[count+2] - m*(1+log(n))
  AWE  <- -2*L[count+2] - 2*m*(3/2 + log(n))
  
  return(list("mean" = mean, "Prob" = Pi, "post" = po,"weights" = w, "classification" = Z, 
              "logLik" = L, "AIC" = ai, "BIC" = ai3, "sigma" = sigma, "fr" = et1, "ICL" = ICL,
              "AICc" = AIcc, "AIC3" = AIC3, "AICu" = AICu, "Caic" = Caic, "AWE" = AWE))
  
}

#### 500
Fe2_500 <- cmw_ran(Y, x, k = 2, maxit = 100, tol = 0.05, show_table = T)
Fe3_500 <- cmw_ran(Y, x, k = 3, maxit = 100, tol = 0.05, show_table = T)
Fe4_500 <- cmw_ran(Y, x, k = 4, maxit = 100, tol = 0.05, show_table = T)
Fe5_500 <- cmw_ran(Y, x, k = 5, maxit = 100, tol = 0.05, show_table = T)

### 1000
Fe2_1000 <- cmw_ran(Y, x, k = 2, maxit = 100, tol = 0.05, show_table = T)
Fe3_1000 <- cmw_ran(Y, x, k = 3, maxit = 100, tol = 0.05, show_table = T)
Fe4_1000 <- cmw_ran(Y, x, k = 4, maxit = 100, tol = 0.05, show_table = T)
Fe5_1000 <- cmw_ran(Y, x, k = 5, maxit = 100, tol = 0.05, show_table = T)


# Confusion Matrix 
cm <- confusionMatrix(sort(as.factor(Fe2_500$classification)), sort(as.factor(z)))

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

## Adjusted Rand Index
clues::adjustedRand(Fe2_500$classification,z)
clues::adjustedRand(Fe3_500$classification,z)
clues::adjustedRand(Fe4_500$classification,z)
clues::adjustedRand(Fe5_500$classification,z)

# Get Estimates
Fe3_500$Prob;  Fe3_500$weights; Fe3_500$mean
Fe3_500$sigma; Fe3_500$fr

# Information criteria
Fe_500$AIC; Fe3_500$AIC; Fe4_500$AIC; Fe5_500$AIC
Fe_500$BIC; Fe3_500$BIC; Fe4_500$BIC; Fe5_500$BIC
Fe_500$ICL; Fe3_500$ICL; Fe4_500$ICL; Fe5_500$ICL
Fe_500$AWE; Fe3_500$AWE; Fe4_500$AWE; Fe5_500$AWE

Fe_500$AIC3;Fe3_500$AIC3;Fe4_500$AIC3;Fe5_500$AIC3
Fe_500$AICc;Fe3_500$AICc;Fe4_500$AICc;Fe5_500$AICc
Fe_500$AICu;Fe3_500$AICu;Fe4_500$AICu;Fe5_500$AICu
Fe_500$Caic;Fe3_500$Caic;Fe4_500$Caic;Fe5_500$Caic

# Multi-class ROC
multiclass.roc(ordered(z),ordered(Fe2_500$classification))
multiclass.roc(ordered(z),ordered(Fe3_500$classification))
multiclass.roc(ordered(z),ordered(Fe4_500$classification))
multiclass.roc(ordered(z),ordered(Fe5_500$classification))

F2 <- ifelse(Fe2_500$classification == 1, 1, 2)  
F3 <- ifelse(Fe3_500$classification == 1, 1, 2)  
F4 <- ifelse(Fe4_500$classification == 1, 1, 2) 
F5 <- ifelse(Fe5_500$classification == 1, 1, 2) 
ZA <- ifelse(z == 1, 1, 2) 

# pr2 <- prediction(F2, ZA)
# pr3 <- prediction(F3, ZA)
# pr4 <- prediction(F4, ZA)
# pr5 <- prediction(F5, ZA)

pr2 <- prediction(Fe2_500$classification, z)
pr3 <- prediction(Fe3_500$classification, z)
pr4 <- prediction(Fe4_500$classification, z)
pr5 <- prediction(Fe5_500$classification, z)

perf2 <- performance(pr2, "tpr", "fpr")
perf3 <- performance(pr3, "tpr", "fpr")
perf4 <- performance(pr4, "tpr", "fpr")
perf5 <- performance(pr5, "tpr", "fpr")

plot(perf2, col = "red", lwd = 3)
plot(perf3, add = TRUE, lwd = 3, col = "blue")
plot(perf4, add = TRUE, col = "black", lwd = 3)
plot(perf5, add = TRUE, col = "Pink", lwd = 3)
legend(0.8,0.4, cex = 1, legend = c("G = 2", "G = 3", "G = 4", "G = 5"),
       col = c("red","blue", "black", "pink"), lwd = 2)


plot(Fe3_500$logLik, type = "l", col = "blue")
plot(da, col = Fe2_11000$classification, pch = Fe2_11000$classification + 1, cex = .5)

### Contour plot
plot(density(da[,5]))
kdf <- kde2d(da[,4], da[,5], n = 50)
image(kdf,col=gray(255:0/255),axes =F);box()
contour(kdf, add = TRUE, col = Fe2_100$classification)

library(rgl)
persp3d(x = kdf, col = Fe3$classification, xlab = "")
snapshot3d(file.path("C:\\Users\\Olobatuyi\\Downloads\\ONGOING WRITEUP\\My Ph.D\\MCWM_classical\\figures",
                     "plot_100cat.png"), top = FALSE)




####################### REAL DATA

cmw_ran <- function(Y, x, k, maxit = 1000, tol = 1e-10, show_table = FALSE){
  
  if(!is.data.frame(x)) x <- unname(as.matrix(x))
  if(!is.data.frame(Y)) Y <- unname(as.matrix(Y))
  
  eps = sqrt(.Machine$double.eps)
  d <- ncol(x); c <- ncol(Y) - 1
  L <- L1 <- L2 <- ai <- ai3 <- lb <- NULL
  
  n <- nrow(x)
  a <- pr <- matrix(NA, nrow = n, ncol = c)
  py <- matrix(NA, nrow = n, ncol = c+1)
  dm <- matrix(NA, nrow = n, ncol = k)
  
  pu <- NULL
  sig <- array(NA, c(d,d,n))
  et1 <- list()
  
  w <- sigma <- list()
  m <- (k * c * (1 + d)) + (k * ((d*d - d)/2 + 2*d + 1))
  
  Pi <- rep(1/k, k)
  
  mean <- matrix(rnorm(d*k), nrow = k)
  
  for(l in 1:k){
    
    sigma[[l]] <- diag(d)
    weights <- matrix(runif(d*c), nrow = c, ncol = d)
    w[[l]] <- weights
    
  }
  
  count <- 2; 
  
  L <- c(-16000,-15000,-14000);
  
  repeat{
    
    for(l in 1:k){
      
      for (i in 1:nrow(x)) {
        
        a[i,]   <- x[i,] %*% t(w[[l]])
        pr[i,]  <- exp(a[i,]) / (1 + sum(exp(a[i,])))
        pu[i]   <- 1 / (1 + sum(exp(a[i,]))); py[i,]  <- c(pu[i], pr[i,])
        dm[i,l] <- dmultinom(Y[i,], sum(Y[i,]), prob = py[i,]) *
          mvnfast::dmvn(x[i,], mean[l,], sigma[[l]], ncores = 8, 
                        isChol = F) * Pi[l]
        
      }
      
    }
    
    po <- dm / rowSums(dm)
    
    # M-Step
    
    # Mean for the group
    
    for(l in 1:k) mean[l, ] <- colSums(po[,l] * x) / colSums(po)[l]
    
    # Update Prior
    Pi <- colSums(po) / nrow(x)
    
    for (l in 1:k) {
      
      L2[l]  <- sum(po[,l] * mvnfast::dmvn(x, mean[l,], sigma[[l]], ncores = 8, 
                                           isChol = F, log = TRUE))
    }
    
    #if(count > 10){
      
      for(l in 1:k){
        
        et <- nnet::multinom(Y ~ x, weights = po[,l], trace = FALSE)
        w[[l]] <- summary(et)$coefficients[,-1]
        lb[l] <- logLik(et)
        et1[[l]] <- summary(et)
        
      }
      
    #}
    
    #if(count > 20){
    
    for (l in 1:k) {

      for (i in 1:nrow(x)) {

        sig[,,i] <- po[i,l] * outer((x[i,] - mean[l,]),(x[i,] - mean[l,]))

      }

      sigma[[l]] <- apply(sig, c(1, 2), sum)/ sum(po[,l]) + diag(eps, ncol(x))

      # which(is.nan(sigma[[l]]))
      # 
      # sigma[[l]][which(is.nan(sigma[[l]]))] <- 1
    }

    #}
    
    # Compute Log_likelihood
    lik  <- sum(log(Pi) * po);lik2 <- sum(L2);lik3 <- sum(lb)
    
    # Compute Log_likelihood lki, lk2, 
    L[count] <- lik + lik2 + lik3
    
    a_k <- (L[count+1] - L[count]) / (L[count] - L[count-1])
    
    L[count + 2] <- L[count] + ((1-a_k)^-1 * (L[count+1] - L[count]))
    
    if (show_table) {
      
      dif <- abs(L[count+2] - L[count+1])
      
      out_table = data.frame(Iteration = count, Likelihood = L[count+2], difference = dif)
      print(kable(out_table))  
      
      if (dif < tol || count == maxit) break; 
      
    }
    
    count <- count + 1
    
  }
  
  Z <- apply(po, 1, which.max)
  
  Prof <- Poc1 <- po
  
  for(i in 1:n){
    
    for(j in 1:k){
      
      Prof[i,j] <- ifelse(po[i,j] > 0.9, 1, 0)
      
      if(po[i,j] > 0){
        
        Poc1[i,j] <- log(po[i,j]) 
        
      }
      
    }
    
  }
  
  ### Information criterion
  
  ai   <- -2*L[count + 2] + 2*m;
  ai3  <- -2*L[count + 2] + m*log(n)
  AIcc <-  ai - 2*m*(m+1)/(n-m-1)
  AIC3 <- -2*L[count+2] - 3*m
  AICu <-  AIcc - n*log(n/(n-m-1))
  ICL  <-  ai3 + suppressWarnings(sum(rowSums(Prof * Poc1)))
  Caic <- -2*L[count+2] - m*(1+log(n))
  AWE  <- -2*L[count+2] - 2*m*(3/2 + log(n))

  return(list("mean" = mean, "Prob" = Pi, "post" = po,"weights" = w, "classification" = Z, 
              "logLik" = L, "AIC" = ai, "BIC" = ai3, "sigma" = sigma, "fr" = et1, "ICL" = ICL,
              "AICc" = AIcc, "AIC3" = AIC3, "AICu" = AICu, "Caic" = Caic, "AWE" = AWE))
  
}

## Contraceptive use among married women
library(dplyr)
contraceptive <- read.table("C:\\Users\\Olobatuyi\\Downloads\\cmcdata.txt", sep = ",")

# Choose the response variable
Y <- contraceptive %>% select(V10)
# Choose the independent variable
x <- contraceptive %>% select(-V10)
# View first four rows
x %>% slice_head(n = 4)

# Perform Pca
xe <- prcomp(x)

plot(xe, type = "l")
summary(xe)

x <- xe$x[,1]

x1 <- cmd[,-10]
split <- sample(1:2, nrow(cmd), prob = c(0.6, 0.40), replace = T)
tr <- which(split==1)
length(tr)
ctr <- cmd[tr,]
dim(ctr)

levels(cmd$V9) <- 0:7
Y <- cmd$V10
Y <- Y - 1
x <- x[tr]
table(Y)

Y <- to_categorical(Y)
head(Y)
x <- cbind(log(ctr$V6+0.5), log(ctr$V7+0.5), log(ctr$V8+0.5),log(ctr$V9+0.5),
           log(ctr$V1),log(ctr$V3+0.5),log(ctr$V9+0.5), log(ctr$V5+0.5)) 

head(x)

da <- data.frame(Y, x)
plot(da)

library(flexmix)
data("dmft")
Y3 <- dmft$End
table(Y3)
levels(dmft$Ethnic) <- 0:2
levels(dmft$Treatment) <- 0:5
levels(dmft$Gender) <- 0:1
head(dmft)

xee <- prcomp(x)

x <- cbind(log(dmft$Begin+0.5), log(as.numeric(dmft$Treatment)+0.5), log(as.numeric(dmft$Ethnic)+0.5),
           log(as.numeric(dmft$Gender)+0.5))

plot(xee, type = 'l')

x1 <- xee$x[,1:2]
Y <- keras::to_categorical(Y3)


Fe2 <- cmw_ran(Y, x, k = 2, maxit = 200, tol = 0.05, show_table = T)
Fe3 <- cmw_ran(Y, x, k = 3, maxit = 200, tol = 0.05, show_table = T)
Fe4 <- cmw_ran(Y, x, k = 4, maxit = 200, tol = 0.05, show_table = T)
Fe5 <- cmw_ran(Y, x, k = 5, maxit = 200, tol = 0.05, show_table = T)


dim(Y)
z <- cmd$V10
table(Fe3$classification)
Fe3$fr
cm <- confusionMatrix(sort(as.factor(Fe3$classification)), sort(as.factor(z)))
#cm <- confusionMatrix((as.factor(Fe2$classification-1)), (as.factor(y2)))
Fe3$Prob

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

## Adjusted Rand Index
clues::adjustedRand(sort(Fe2$classification-1),sort(y2))
clues::adjustedRand(sort(Fe3$classification),sort(z))
clues::adjustedRand(Fe4_1000$classification,z)
clues::adjustedRand(Fe5_1000$classification,z)

da <- data.frame(Y, x)
plot(da, col = sort(Fe3$classification), pch = sort(Fe3$classification) + 1, cex = .5)
plot(x, col = sort(z))

X <- as.data.frame(xt)
par(mfrow = c(1,2))
V <- lda(X, clust)$scaling
x.fda <- as.matrix(X) %*% V
plot(x.fda, type = "n", xlab = "", ylab = "", pch = 1, col = cluster, cex = .2)
legend(x = -6, y = 4, cex = 0.5, legend = c(3,5,8))
text(x.fda, col = cluster, labels = c(3,5,8)[cluster], cex = .2)


#r <- PCA(cmd[,-10], ncp = 5, graph = FALSE)
r <- PCA(tsne_out$Y, ncp = 5, graph = FALSE)

fviz_pca_ind(r, geom = "point", axes = c(1,2), 
             habillage = (as.factor(Fe3$classification)),
             palette = c("red", "blue", "gold","magenta", "black", "cyan","orange", 
                         "grey", "green","brown"),
             addEllipse = T, repel = F,
             ggtheme = theme_minimal(), title = "")

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

pred <- sort(as.factor(Fe2$classification))
colors = rainbow(length(unique(pred)))
names(colors) = unique(pred)
ecb = function(x){
  plot(x, t ='p', col=colors[e6]); 
  #text(x,labels=clust, col=colors[clust+1]) 
}

tsne_iris = tsne(x, epoch_callback = ecb, k = 6, perplexity = 20, max_iter = 1000)

rm(list = ls(all = TRUE))
# compare to PCA
dev.new()
pca_iris = princomp(xe$x)$scores[,1:2]
plot(pca_iris, t='n')
text(pca_iris, labels=pred,col=colors[pred])


colors <- rainbow(length(unique(pred)))
names(colors) <- unique(pred)

ecb = function(x){
  plot(x, t ='n');
  text(x, labels = pred, col = colors[pred], cex = .5)
}

tsne_iris = tsne(x2, k = 1, epoch_callback = ecb, perplexity=5, max_iter = 300)

x2 <- cbind(x, cluster)

qw <- colors()

table(clust)
w1 <- crep(w1[1], 1006)

cl <- fit$cluster
colors = rainbow(length(unique(pred)))
names(colors) = unique(pred)
# data.frame as input
cd <- as.data.frame(cmd[,-10])
tsne_out <- Rtsne(cmd[,-10], pca = F, theta=0.0, check_duplicates = FALSE)

#tol <- read.csv("tolu.csv")
# dim(tol)
# tol <- tol[,-1]
#text(tsne_out$Y, labels = pred, col = colors[pred], asp = 0, )

r <- PCA(tsne_out$Y, ncp = 5, graph = FALSE)

fviz_pca_ind(r, geom = "point", axes = c(1,3), 
             habillage = pred,
             addEllipse = FALSE, repel = TRUE,
             ggtheme = theme_minimal(), title = "",
             xlab = "", ylab = "")

############# Heart Data

heart <- read.table("D:\\traindcm\\Data\\Heart_disease\\processed.cleveland.data", sep = ",")

dim(heart)
head(heart)

heart[,12:13] <- with(heart, c(V12 <- as.numeric(levels(V12))[V12],
                               V13 <- as.numeric(levels(V13))[V13]))

## Replace all NA with 0
heart[which(is.na(heart$V12) == TRUE), 12] = 0
heart[which(is.na(heart$V13) == TRUE), 13] = 0

head(heart)
sum(is.na(heart$V14))

x <- as.matrix(heart)

cor(heart[,-14])
x <- x[, -14]

pairs(heart[,-14])

str(heart)

hist(heart)
?PCA
RT <- prcomp(x, scale = TRUE)
x <- RT$x[,1]
plot(RT, type = "l")
summary(RT)
# library(scatterplot3d)
# library(rgl)
# library(car)
# library(plot3D)

#data(VADeaths)
#  hist3D and ribbon3D with greyish background, rotated, rescaled,...
# hist3D(z = VADeaths, scale = FALSE, expand = 0.01, bty = "g", phi = 20,
#         col = "#0072B2", border = "black", shade = 0.2, ltheta = 90,
#         space = 0.3, ticktype = "detailed", d = 2)
# 
# 
# hist3D (x = 1:5, y = 1:4, z = VADeaths,
#         bty = "g", phi = 20,  theta = -60,
#         xlab = "", ylab = "", zlab = "", main = "VADeaths",
#         col = "#0072B2", border = "black", shade = 0.8,
#         ticktype = "detailed", space = 0.15, d = 2, cex.axis = 1e-9)
# # Use text3D to label x axis
#  text3D(x = 1:5, y = rep(0.5, 5), z = rep(3, 5),
#        labels = rownames(VADeaths),
#        add = TRUE, adj = 0)
# # Use text3D to label y axis
#  text3D(x = rep(1, 4),   y = 1:4, z = rep(0, 4),
#        labels  = colnames(VADeaths),
#        add = TRUE, adj = 1)

y1 <- heart[,14]

y2 <- ifelse(y1 == 0, 0, 1)
table(y2)

# One Hot encoding
Y <- to_categorical(y)




########### FOR MNIST not working yet 
cmw_ran_batch <- function(Y, x, k, Batch = 50, maxit = 1000, tol = 1e-10, show_table = FALSE){
  
  if(!is.data.frame(x)) x <- unname(as.matrix(x))
  if(!is.data.frame(Y)) Y <- unname(as.matrix(Y))
  
  eps = sqrt(.Machine$double.eps)
  d <- ncol(x); c <- ncol(Y) - 1
  L <- L1 <- L2 <- ai <- ai3 <- lb <- NULL
  
  #n <- nrow(x)
  a <- pr <- matrix(NA, nrow = Batch, ncol = c)
  py <- matrix(NA, nrow = Batch, ncol = c+1)
  dm <- matrix(NA, nrow = Batch, ncol = k)
  
  pu <- NULL
  sig <- array(NA, c(d,d,Batch))
  et1 <- list()
  
  w <- sigma <- list()
  m <- (k * c * (1 + d)) + (k * ((d*d - d)/2 + 2*d + 1))
  
  Pi <- rep(1/k, k)
  
  mean <- matrix(rnorm(d*k), nrow = k)
  
  for(l in 1:k){
    
    sigma[[l]] <- diag(d)
    weights <- matrix(runif(d*c), nrow = c, ncol = d)
    w[[l]] <- weights
    
  }
  
  count <- 2; 
  
  L <- c(-16000,-15000,-14000);
  
  repeat{
      
      set.seed(1)
      batchIndex <- sample(1:nrow(x), size = Batch, replace = T)
      batchx <- as.matrix(x[batchIndex, ])
      batchY <- Y[batchIndex, ]
      
      for(l in 1:k){
        
        for (i in 1:nrow(batchx)) {
          
          a[i,]   <- batchx[i,] %*% t(w[[l]])
          pr[i,]  <- exp(a[i,]) / (1 + sum(exp(a[i,])))
          pu[i]   <- 1 / (1 + sum(exp(a[i,]))); py[i,]  <- c(pu[i], pr[i,])
          dm[i,l] <- dmultinom(batchY[i,], sum(batchY[i,]), prob = py[i,]) *
            mvnfast::dmvn(batchx[i,], mean[l,], sigma[[l]], ncores = 8, 
                          isChol = F) * Pi[l]
          
        }
        
      }
      
      po <- dm / rowSums(dm)
      
      # M-Step
      
      # Mean for the group
      
      for(l in 1:k) mean[l, ] <- colSums(po[,l] * batchx) / colSums(po)[l]
      
      # Update Prior
      Pi <- colSums(po) / nrow(batchx)
      
      for (l in 1:k) {
        
        L2[l]  <- sum(po[,l] * mvnfast::dmvn(batchx, mean[l,], sigma[[l]], ncores = 8, 
                                             isChol = F, log = TRUE))
      }
      
    # if(count > 10){
    #   
    #   for(l in 1:k){
    #     
    #     et <- nnet::multinom(batchY ~ batchx, weights = po[,l], trace = FALSE)
    #     w[[l]] <- summary(et)$coefficients[,-1]
    #     lb[l] <- logLik(et)
    #     et1[[l]] <- summary(et)
    #     
    #   }
    #   
    # }
      
      #if(count > 20){
      
      for (l in 1:k) {
        
        for (i in 1:nrow(batchx)) {
          
          sig[,,i] <- po[i,l] * outer((batchx[i,] - mean[l,]),(batchx[i,] - mean[l,]))
          
        }
        
        sigma[[l]] <- apply(sig, c(1, 2), sum)/ sum(po[,l]) + diag(eps, ncol(batchx))
        
        # which(is.nan(sigma[[l]]))
        # 
        # sigma[[l]][which(is.nan(sigma[[l]]))] <- 1
      }
      
      #}
      
      # Compute Log_likelihood
      lik  <- sum(log(Pi) * po);lik2 <- sum(L2);lik3 <- sum(lb)
      
      # Compute Log_likelihood lki, lk2, 
      L[count] <- lik + lik2 + lik3
      
      a_k <- (L[count+1] - L[count]) / (L[count] - L[count-1])
      
      L[count + 2] <- L[count] + ((1-a_k)^-1 * (L[count+1] - L[count]))
      
      if (show_table) {
        
        dif <- abs(L[count+2] - L[count+1])
        
        out_table = data.frame(Iteration = count, Likelihood = L[count+2], difference = dif)
        print(kable(out_table))  
        
        if (dif < tol || count == maxit) break; 
        
      }
      
      count <- count + 1
      
    }
  
  Z <- apply(po, 1, which.max)
  
  Prof <- Poc1 <- po
  
  for(i in 1:nrow(batchx)){
    
    for(j in 1:k){
      
      Prof[i,j] <- ifelse(po[i,j] > 0.9, 1, 0)
      
      if(po[i,j] > 0){
        
        Poc1[i,j] <- log(po[i,j]) 
        
      }
      
    }
    
  }
  
  ### Information criterion
  
  # ai   <- -2*L[count + 2] + 2*m;
  # ai3  <- -2*L[count + 2] + m*log(n)
  # AIcc <-  ai - 2*m*(m+1)/(n-m-1)
  # AIC3 <- -2*L[count+2] - 3*m
  # AICu <-  AIcc - n*log(n/(n-m-1))
  # ICL  <-  ai3 + suppressWarnings(sum(rowSums(Prof * Poc1)))
  # Caic <- -2*L[count+2] - m*(1+log(n))
  # AWE  <- -2*L[count+2] - 2*m*(3/2 + log(n))
  
  return(list("mean" = mean, "Prob" = Pi, "post" = po,"weights" = w, "sigma" = sigma,
              "classification" = Z))
  
}

x <- as.matrix(x)

Fe5 <- cmw_ran_batch(Y, x, k = 10, Batch = 100, maxit = 100, tol = 0.05, show_table = T)
table(trainy)
=======
library(knitr)
library(gridExtra);library(grid);require(keras);library(tidyverse);
library(OpenMPController);library(caret);library(ggplot2);library(SIBER)
library(factoextra);library(FactoMineR);library(xtable);library(clues)
library(MBCbook);library(HDclassif);library(dplyr);library(SIBER)
require("markdown");require("rattle");require("xtable");require("stringr")
require("fBasics");require("MASS");require("survival");require("STAR")
require("gamlss.dist");require("VGAM");library(rgl);library(ellipse)
library(threejs);library(plotROC);library(ROCR);library(pROC)


# Set the number of column of X and fix G

d <- 2; k <- 2; c <- 3-1; n <- 200; ni <- 1
df <- list(); w <- sigma <- list()

Pi <- c(0.6, 0.4); mean <- matrix(c(.1,-2,2,1), nrow = k)

for(l in 1:k) {

  sigma[[l]] <- diag(d)

}

weights1 <- matrix(c(.2,.2,.2,1,1,1,1,1), nrow = d, ncol = c)
weights2 <- matrix(c(1,1,1,1,1,1,1,1), nrow = d, ncol = c)
w[[1]] <- t(weights1); w[[2]] <- t(weights2)

z <- NULL
# Generating Algorithm
x <- matrix(NA, nrow = n, ncol = d)
a <- pr <- matrix(NA, nrow = n, ncol = c)
Y <- py <- matrix(NA, nrow = n, ncol = c+1)
set.seed(20); U <- runif(n)

for (i in 1:n) {

  if(U[i] <= Pi[1]){

    z[i] <- 1
    x[i,] <- mvnfast::rmvn(1, mean[1,], sigma[[1]])
    a[i,]    <- as.matrix(x[i,] %*% t(w[[1]]))
    pr[i,]   <- exp(a[i,]) / (1 + sum(exp(a[i,])))
    pu       <- 1 / (1 + sum(exp(a[i,])))
    py[i,]   <- c(pu, pr[i,])
    Y[i,] <- t(rmultinom(1, ni, py[i,]))
    #ni <- ni + i

  }else{

    z[i] <- 2
    x[i,] <- mvnfast::rmvn(1, mean[2,], sigma[[2]])
    a[i,]    <- as.matrix(x[i,] %*% t(w[[2]]))
    pr[i,]   <- exp(a[i,]) / (1 + sum(exp(a[i,])))
    pu       <- 1 / (1 + sum(exp(a[i,])))
    py[i,]   <- c(pu, pr[i,])
    Y[i,] <- t(rmultinom(1, ni, py[i,]))
    #ni <- ni + i
  }

}


da <- data.frame(Y, x)
names(da) <- c("a", "b", "c", "x1", 'x2')
plot(da, col = z)

pos1 <- which(z == 1)
pos  <- which(z == 2)

Y1 <- Y[pos1,];Y2 <- Y[pos,]

colSums(Y1);colSums(Y2)
table(z)

# The number of a, b, c in group 1 is 10,44, 3
# The number of a, b, c in group 2 is 8,10, 25


##################### K = 2

# Set the number of column of X and fix G
d <- 2; k <- 2; ni <- 1; c <- 3-1; n <- 1000

# create an empty list
w <- sigma <- list()
z <- NULL
x <- matrix(NA, nrow = n, ncol = d)
a <- pr <- matrix(NA, nrow = n, ncol = c)
Y <- py <- matrix(NA, nrow = n, ncol = c+1)

# Z parameter
Pi <- rep(1/k,k)
#Pi <- c(5,3,2)
# X parameter
mean <- matrix(c(.1,-2,2,0,1,3), nrow = k, ncol = d)
for(l in 1:k) {sigma[[l]] <- diag(d)}

# Y parameter
# weights1 <- matrix(c(.1,.2,.3,.1), nrow = d, ncol = c);w[[1]] <- t(weights1)
# weights2 <- matrix(c(1,1,1,1), nrow = d, ncol = c);w[[2]] <- t(weights2)
# weights3 <- matrix(c(.1,.1,.1,.1), nrow = d, ncol = c);w[[3]] <- t(weights3)

weights1 <- matrix(c(5,.4,.3,0.04), nrow = d, ncol = c);    w[[1]] <- t(weights1)
weights2 <- matrix(c(.01,.02,2,1), nrow = d, ncol = c);     w[[2]] <- t(weights2)
weights3 <- matrix(c(1, 0.03,0.06,.02), nrow = d, ncol = c);w[[3]] <- t(weights3)

# Set a seed
set.seed(20);U <- runif(n)

for (i in 1:n) {
  
  if(U[i] <= Pi[1]){
    
    z[i] <- 1
    x[i,] <- mvnfast::rmvn(1, mean[1,], sigma[[1]])
    a[i,]    <- as.matrix(x[i,] %*% t(w[[1]]))
    pr[i,]   <- exp(a[i,]) / (1 + sum(exp(a[i,])))
    pu       <- 1 / (1 + sum(exp(a[i,])))
    py[i,]   <- c(pu, pr[i,])
    Y[i,] <- t(rmultinom(1, ni, py[i,]))
    
  }else{
      
      z[i]   <- 2
      x[i,]  <- mvnfast::rmvn(1, mean[2,], sigma[[2]])
      a[i,]  <- as.matrix(x[i,] %*% t(w[[2]]))
      pr[i,] <- exp(a[i,]) / (1 + sum(exp(a[i,])))
      pu     <- 1 / (1 + sum(exp(a[i,])))
      py[i,] <- c(pu, pr[i,])
      Y[i,]  <- t(rmultinom(1, ni, py[i,]))
      
    }
  }

da <- data.frame(Y, x)
names(da) <- c("a", "b", "c", "x1", 'x2')
plot(da, col = z)

pos1 <- which(z == 1)
pos2 <- which(z == 2)

Y1 <- Y[pos1,];Y2 <- Y[pos2,];#Y3 <- Y[pos3,]
colSums(Y1);colSums(Y2)#;colSums(Y3)
table(z)


##################### K = 3

# Set the number of column of X and fix G
d <- 2; k <- 3; ni <- 1; c <- 3-1; n <- 1000

# create an empty list
w <- sigma <- list()
z <- NULL
x <- matrix(NA, nrow = n, ncol = d)
a <- pr <- matrix(NA, nrow = n, ncol = c)
Y <- py <- matrix(NA, nrow = n, ncol = c+1)

# Z parameter
Pi <- rep(1/k,k)
#Pi <- c(5,3,2)
# X parameter
mean <- matrix(c(.1,-2,2,0,1,3), nrow = k, ncol = d)
for(l in 1:k) {sigma[[l]] <- diag(d)}

# Y parameter
# weights1 <- matrix(c(.1,.2,.3,.1), nrow = d, ncol = c);w[[1]] <- t(weights1)
# weights2 <- matrix(c(1,1,1,1), nrow = d, ncol = c);w[[2]] <- t(weights2)
# weights3 <- matrix(c(.1,.1,.1,.1), nrow = d, ncol = c);w[[3]] <- t(weights3)

weights1 <- matrix(c(5, .4, .3, 0.04), nrow = d, ncol = c);     w[[1]] <- t(weights1)
weights2 <- matrix(c(.01, .02, 2,  1), nrow = d, ncol = c);      w[[2]] <- t(weights2)
weights3 <- matrix(c(1,0.03,0.06,.02), nrow = d, ncol = c); w[[3]] <- t(weights3)

# Set a seed
set.seed(20);U <- runif(n)

for (i in 1:n) {

  if(U[i] <= Pi[1]){

    z[i] <- 1
    x[i,] <- mvnfast::rmvn(1, mean[1,], sigma[[1]])
    a[i,]    <- as.matrix(x[i,] %*% t(w[[1]]))
    pr[i,]   <- exp(a[i,]) / (1 + sum(exp(a[i,])))
    pu       <- 1 / (1 + sum(exp(a[i,])))
    py[i,]   <- c(pu, pr[i,])
    Y[i,] <- t(rmultinom(1, ni, py[i,]))
    
  }else

    if(U[i] > Pi[1] & U[i] < (Pi[2] + Pi[1])){

      z[i]   <- 2
      x[i,]  <- mvnfast::rmvn(1, mean[2,], sigma[[2]])
      a[i,]  <- as.matrix(x[i,] %*% t(w[[2]]))
      pr[i,] <- exp(a[i,]) / (1 + sum(exp(a[i,])))
      pu     <- 1 / (1 + sum(exp(a[i,])))
      py[i,] <- c(pu, pr[i,])
      Y[i,]  <- t(rmultinom(1, ni, py[i,]))
      
    }else{

      z[i] <- 3
      x[i,] <- mvnfast::rmvn(1, mean[3,], sigma[[3]])
      a[i,]    <- as.matrix(x[i,] %*% t(w[[3]]))
      pr[i,]   <- exp(a[i,]) / (1 + sum(exp(a[i,])))
      pu       <- 1 / (1 + sum(exp(a[i,])))
      py[i,]   <- c(pu, pr[i,])
      Y[i,] <- t(rmultinom(1, ni, py[i,]))

    }

}

da <- data.frame(Y, x)
names(da) <- c("a", "b", "c", "x1", 'x2')
plot(da, col = z)

pos1 <- which(z == 1)
pos2 <- which(z == 2)
pos3 <- which(z == 3)

Y1 <- Y[pos1,];Y2 <- Y[pos2,];Y3 <- Y[pos3,]
colSums(Y1);colSums(Y2);colSums(Y3)
table(z)



cmw_ran <- function(Y, x, k, maxit = 1000, tol = 1e-10, show_table = FALSE){
  
  if(!is.data.frame(x)) x <- unname(as.matrix(x))
  if(!is.data.frame(Y)) Y <- unname(as.matrix(Y))
  
  eps = sqrt(.Machine$double.eps)
  d <- ncol(x); c <- ncol(Y) - 1
  L <- L1 <- L2 <- ai <- ai3 <- lb <- NULL
  
  n <- nrow(x)
  a <- pr <- matrix(NA, nrow = n, ncol = c)
  py <- matrix(NA, nrow = n, ncol = c+1)
  dm <- matrix(NA, nrow = n, ncol = k)
  
  pu <- NULL
  sig <- array(NA, c(d,d,n))
  et1 <- list()
  
  w <- sigma <- list()
  m <- (k * c * (1 + d)) + (k * ((d*d - d)/2 + 2*d + 1))
  
  Pi <- rep(1/k, k)
  
  mean <- matrix(rnorm(d*k), nrow = k)
  
  for(l in 1:k){
    
    sigma[[l]] <- diag(d)
    weights <- matrix(runif(d*c), nrow = c, ncol = d)
    w[[l]] <- weights
    
  }
  
  count <- 2; 
  
  L <- c(-16000,-15000,-14000);
  
  repeat{
    
    for(l in 1:k){
      
      for (i in 1:nrow(x)) {
        
        a[i,]   <- x[i,] %*% t(w[[l]])
        pr[i,]  <- exp(a[i,]) / (1 + sum(exp(a[i,])))
        pu[i]   <- 1 / (1 + sum(exp(a[i,]))); py[i,]  <- c(pu[i], pr[i,])
        dm[i,l] <- dmultinom(Y[i,], sum(Y[i,]), prob = py[i,]) *
          mvnfast::dmvn(x[i,], mean[l,], sigma[[l]], ncores = 8, 
                        isChol = F) * Pi[l]
        
      }
      
    }
    
    po <- dm / rowSums(dm)
    
    # M-Step
    
    # Mean for the group
    
    for(l in 1:k) mean[l, ] <- colSums(po[,l] * x) / colSums(po)[l]
    
    # Update Prior
    Pi <- colSums(po) / nrow(x)
    
    for (l in 1:k) {
      
      L2[l]  <- sum(po[,l] * mvnfast::dmvn(x, mean[l,], sigma[[l]], ncores = 8, 
                                           isChol = F, log = TRUE))
    }
    
    if(count > 10){
      
      for(l in 1:k){
        
        et <- nnet::multinom(Y ~ x, weights = po[,l], trace = FALSE)
        w[[l]] <- summary(et)$coefficients[,-1]
        lb[l] <- logLik(et)
        et1[[l]] <- summary(et)
        
      }
      
    }
    
    #if(count > 20){
      
      for (l in 1:k) {
        
        for (i in 1:nrow(x)) {
          
          sig[,,i] <- po[i,l] * outer((x[i,] - mean[l,]),(x[i,] - mean[l,]))
          
        }
        
        sigma[[l]] <- apply(sig, c(1, 2), sum)/ sum(po[,l]) #+ diag(eps, ncol(x))
        
        which(is.nan(sigma[[l]]))
        
        sigma[[l]][which(is.nan(sigma[[l]]))] <- 1
      }
      
    #}
    
    # Compute Log_likelihood
    lik  <- sum(log(Pi) * po);lik2 <- sum(L2);lik3 <- sum(lb)
    
    # Compute Log_likelihood lki, lk2, 
    L[count] <- lik + lik2 + lik3
    
    a_k <- (L[count+1] - L[count]) / (L[count] - L[count-1])
    
    L[count + 2] <- L[count] + ((1-a_k)^-1 * (L[count+1] - L[count]))
    
    if (show_table) {
      
      dif <- abs(L[count+2] - L[count+1])
      
      out_table = data.frame(Iteration = count, Likelihood = L[count+2], difference = dif)
      print(kable(out_table))  
      
      if (dif < tol || count == maxit) break; 
      
    }
    
    count <- count + 1
    
  }
  
  Z <- apply(po, 1, which.max)
  
  Prof <- Poc1 <- po
  
  for(i in 1:n){
    
    for(j in 1:k){
      
      Prof[i,j] <- ifelse(po[i,j] > 0.9, 1, 0)
      
      if(po[i,j] > 0){
        
        Poc1[i,j] <- log(po[i,j]) 
        
      }
      
    }
    
  }
  
  ### Information criterion
  
  ai   <- -2*L[count + 2] + 2*m;
  ai3  <- -2*L[count + 2] + m*log(n)  
  AIcc <-  ai - 2*m*(m+1)/(n-m-1)
  AIC3 <- -2*L[count+2] - 3*m
  AICu <-  AIcc - n*log(n/(n-m-1))
  ICL  <-  ai3 + suppressWarnings(sum(rowSums(Prof * Poc1)))
  Caic <- -2*L[count+2] - m*(1+log(n))
  AWE  <- -2*L[count+2] - 2*m*(3/2 + log(n))
  
  return(list("mean" = mean, "Prob" = Pi, "post" = po,"weights" = w, "classification" = Z, 
              "logLik" = L, "AIC" = ai, "BIC" = ai3, "sigma" = sigma, "fr" = et1, "ICL" = ICL,
              "AICc" = AIcc, "AIC3" = AIC3, "AICu" = AICu, "Caic" = Caic, "AWE" = AWE))
  
}

#### 500
Fe2_500 <- cmw_ran(Y, x, k = 2, maxit = 100, tol = 0.05, show_table = T)
Fe3_500 <- cmw_ran(Y, x, k = 3, maxit = 100, tol = 0.05, show_table = T)
Fe4_500 <- cmw_ran(Y, x, k = 4, maxit = 100, tol = 0.05, show_table = T)
Fe5_500 <- cmw_ran(Y, x, k = 5, maxit = 100, tol = 0.05, show_table = T)

### 1000
Fe2_1000 <- cmw_ran(Y, x, k = 2, maxit = 100, tol = 0.05, show_table = T)
Fe3_1000 <- cmw_ran(Y, x, k = 3, maxit = 100, tol = 0.05, show_table = T)
Fe4_1000 <- cmw_ran(Y, x, k = 4, maxit = 100, tol = 0.05, show_table = T)
Fe5_1000 <- cmw_ran(Y, x, k = 5, maxit = 100, tol = 0.05, show_table = T)


# Confusion Matrix 
cm <- confusionMatrix(sort(as.factor(Fe2_500$classification)), sort(as.factor(z)))

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

## Adjusted Rand Index
clues::adjustedRand(Fe2_500$classification,z)
clues::adjustedRand(Fe3_500$classification,z)
clues::adjustedRand(Fe4_500$classification,z)
clues::adjustedRand(Fe5_500$classification,z)

# Get Estimates
Fe3_500$Prob;  Fe3_500$weights; Fe3_500$mean
Fe3_500$sigma; Fe3_500$fr

# Information criteria
Fe_500$AIC; Fe3_500$AIC; Fe4_500$AIC; Fe5_500$AIC
Fe_500$BIC; Fe3_500$BIC; Fe4_500$BIC; Fe5_500$BIC
Fe_500$ICL; Fe3_500$ICL; Fe4_500$ICL; Fe5_500$ICL
Fe_500$AWE; Fe3_500$AWE; Fe4_500$AWE; Fe5_500$AWE

Fe_500$AIC3;Fe3_500$AIC3;Fe4_500$AIC3;Fe5_500$AIC3
Fe_500$AICc;Fe3_500$AICc;Fe4_500$AICc;Fe5_500$AICc
Fe_500$AICu;Fe3_500$AICu;Fe4_500$AICu;Fe5_500$AICu
Fe_500$Caic;Fe3_500$Caic;Fe4_500$Caic;Fe5_500$Caic

# Multi-class ROC
multiclass.roc(ordered(z),ordered(Fe2_500$classification))
multiclass.roc(ordered(z),ordered(Fe3_500$classification))
multiclass.roc(ordered(z),ordered(Fe4_500$classification))
multiclass.roc(ordered(z),ordered(Fe5_500$classification))

F2 <- ifelse(Fe2_500$classification == 1, 1, 2)  
F3 <- ifelse(Fe3_500$classification == 1, 1, 2)  
F4 <- ifelse(Fe4_500$classification == 1, 1, 2) 
F5 <- ifelse(Fe5_500$classification == 1, 1, 2) 
ZA <- ifelse(z == 1, 1, 2) 

# pr2 <- prediction(F2, ZA)
# pr3 <- prediction(F3, ZA)
# pr4 <- prediction(F4, ZA)
# pr5 <- prediction(F5, ZA)

pr2 <- prediction(Fe2_500$classification, z)
pr3 <- prediction(Fe3_500$classification, z)
pr4 <- prediction(Fe4_500$classification, z)
pr5 <- prediction(Fe5_500$classification, z)

perf2 <- performance(pr2, "tpr", "fpr")
perf3 <- performance(pr3, "tpr", "fpr")
perf4 <- performance(pr4, "tpr", "fpr")
perf5 <- performance(pr5, "tpr", "fpr")

plot(perf2, col = "red", lwd = 3)
plot(perf3, add = TRUE, lwd = 3, col = "blue")
plot(perf4, add = TRUE, col = "black", lwd = 3)
plot(perf5, add = TRUE, col = "Pink", lwd = 3)
legend(0.8,0.4, cex = 1, legend = c("G = 2", "G = 3", "G = 4", "G = 5"),
       col = c("red","blue", "black", "pink"), lwd = 2)


plot(Fe3_500$logLik, type = "l", col = "blue")
plot(da, col = Fe2_11000$classification, pch = Fe2_11000$classification + 1, cex = .5)

### Contour plot
plot(density(da[,5]))
kdf <- kde2d(da[,4], da[,5], n = 50)
image(kdf,col=gray(255:0/255),axes =F);box()
contour(kdf, add = TRUE, col = Fe2_100$classification)

library(rgl)
persp3d(x = kdf, col = Fe3$classification, xlab = "")
snapshot3d(file.path("C:\\Users\\Olobatuyi\\Downloads\\ONGOING WRITEUP\\My Ph.D\\MCWM_classical\\figures",
                     "plot_100cat.png"), top = FALSE)




####################### REAL DATA

cmw_ran <- function(Y, x, k, maxit = 1000, tol = 1e-10, show_table = FALSE){
  
  if(!is.data.frame(x)) x <- unname(as.matrix(x))
  if(!is.data.frame(Y)) Y <- unname(as.matrix(Y))
  
  eps = sqrt(.Machine$double.eps)
  d <- ncol(x); c <- ncol(Y) - 1
  L <- L1 <- L2 <- ai <- ai3 <- lb <- NULL
  
  n <- nrow(x)
  a <- pr <- matrix(NA, nrow = n, ncol = c)
  py <- matrix(NA, nrow = n, ncol = c+1)
  dm <- matrix(NA, nrow = n, ncol = k)
  
  pu <- NULL
  sig <- array(NA, c(d,d,n))
  et1 <- list()
  
  w <- sigma <- list()
  m <- (k * c * (1 + d)) + (k * ((d*d - d)/2 + 2*d + 1))
  
  Pi <- rep(1/k, k)
  
  mean <- matrix(rnorm(d*k), nrow = k)
  
  for(l in 1:k){
    
    sigma[[l]] <- diag(d)
    weights <- matrix(runif(d*c), nrow = c, ncol = d)
    w[[l]] <- weights
    
  }
  
  count <- 2; 
  
  L <- c(-16000,-15000,-14000);
  
  repeat{
    
    for(l in 1:k){
      
      for (i in 1:nrow(x)) {
        
        a[i,]   <- x[i,] %*% t(w[[l]])
        pr[i,]  <- exp(a[i,]) / (1 + sum(exp(a[i,])))
        pu[i]   <- 1 / (1 + sum(exp(a[i,]))); py[i,]  <- c(pu[i], pr[i,])
        dm[i,l] <- dmultinom(Y[i,], sum(Y[i,]), prob = py[i,]) *
          mvnfast::dmvn(x[i,], mean[l,], sigma[[l]], ncores = 8, 
                        isChol = F) * Pi[l]
        
      }
      
    }
    
    po <- dm / rowSums(dm)
    
    # M-Step
    
    # Mean for the group
    
    for(l in 1:k) mean[l, ] <- colSums(po[,l] * x) / colSums(po)[l]
    
    # Update Prior
    Pi <- colSums(po) / nrow(x)
    
    for (l in 1:k) {
      
      L2[l]  <- sum(po[,l] * mvnfast::dmvn(x, mean[l,], sigma[[l]], ncores = 8, 
                                           isChol = F, log = TRUE))
    }
    
    #if(count > 10){
      
      for(l in 1:k){
        
        et <- nnet::multinom(Y ~ x, weights = po[,l], trace = FALSE)
        w[[l]] <- summary(et)$coefficients[,-1]
        lb[l] <- logLik(et)
        et1[[l]] <- summary(et)
        
      }
      
    #}
    
    #if(count > 20){
    
    for (l in 1:k) {

      for (i in 1:nrow(x)) {

        sig[,,i] <- po[i,l] * outer((x[i,] - mean[l,]),(x[i,] - mean[l,]))

      }

      sigma[[l]] <- apply(sig, c(1, 2), sum)/ sum(po[,l]) + diag(eps, ncol(x))

      # which(is.nan(sigma[[l]]))
      # 
      # sigma[[l]][which(is.nan(sigma[[l]]))] <- 1
    }

    #}
    
    # Compute Log_likelihood
    lik  <- sum(log(Pi) * po);lik2 <- sum(L2);lik3 <- sum(lb)
    
    # Compute Log_likelihood lki, lk2, 
    L[count] <- lik + lik2 + lik3
    
    a_k <- (L[count+1] - L[count]) / (L[count] - L[count-1])
    
    L[count + 2] <- L[count] + ((1-a_k)^-1 * (L[count+1] - L[count]))
    
    if (show_table) {
      
      dif <- abs(L[count+2] - L[count+1])
      
      out_table = data.frame(Iteration = count, Likelihood = L[count+2], difference = dif)
      print(kable(out_table))  
      
      if (dif < tol || count == maxit) break; 
      
    }
    
    count <- count + 1
    
  }
  
  Z <- apply(po, 1, which.max)
  
  Prof <- Poc1 <- po
  
  for(i in 1:n){
    
    for(j in 1:k){
      
      Prof[i,j] <- ifelse(po[i,j] > 0.9, 1, 0)
      
      if(po[i,j] > 0){
        
        Poc1[i,j] <- log(po[i,j]) 
        
      }
      
    }
    
  }
  
  ### Information criterion
  
  ai   <- -2*L[count + 2] + 2*m;
  ai3  <- -2*L[count + 2] + m*log(n)
  AIcc <-  ai - 2*m*(m+1)/(n-m-1)
  AIC3 <- -2*L[count+2] - 3*m
  AICu <-  AIcc - n*log(n/(n-m-1))
  ICL  <-  ai3 + suppressWarnings(sum(rowSums(Prof * Poc1)))
  Caic <- -2*L[count+2] - m*(1+log(n))
  AWE  <- -2*L[count+2] - 2*m*(3/2 + log(n))

  return(list("mean" = mean, "Prob" = Pi, "post" = po,"weights" = w, "classification" = Z, 
              "logLik" = L, "AIC" = ai, "BIC" = ai3, "sigma" = sigma, "fr" = et1, "ICL" = ICL,
              "AICc" = AIcc, "AIC3" = AIC3, "AICu" = AICu, "Caic" = Caic, "AWE" = AWE))
  
}

## Contraceptive use among married women
library(dplyr)
contraceptive <- read.table("C:\\Users\\Olobatuyi\\Downloads\\cmcdata.txt", sep = ",")

# Choose the response variable
Y <- contraceptive %>% select(V10)
# Choose the independent variable
x <- contraceptive %>% select(-V10)
# View first four rows
x %>% slice_head(n = 4)

# Perform Pca
xe <- prcomp(x)

plot(xe, type = "l")
summary(xe)

x <- xe$x[,1]

x1 <- cmd[,-10]
split <- sample(1:2, nrow(cmd), prob = c(0.6, 0.40), replace = T)
tr <- which(split==1)
length(tr)
ctr <- cmd[tr,]
dim(ctr)

levels(cmd$V9) <- 0:7
Y <- cmd$V10
Y <- Y - 1
x <- x[tr]
table(Y)

Y <- to_categorical(Y)
head(Y)
x <- cbind(log(ctr$V6+0.5), log(ctr$V7+0.5), log(ctr$V8+0.5),log(ctr$V9+0.5),
           log(ctr$V1),log(ctr$V3+0.5),log(ctr$V9+0.5), log(ctr$V5+0.5)) 

head(x)

da <- data.frame(Y, x)
plot(da)

library(flexmix)
data("dmft")
Y3 <- dmft$End
table(Y3)
levels(dmft$Ethnic) <- 0:2
levels(dmft$Treatment) <- 0:5
levels(dmft$Gender) <- 0:1
head(dmft)

xee <- prcomp(x)

x <- cbind(log(dmft$Begin+0.5), log(as.numeric(dmft$Treatment)+0.5), log(as.numeric(dmft$Ethnic)+0.5),
           log(as.numeric(dmft$Gender)+0.5))

plot(xee, type = 'l')

x1 <- xee$x[,1:2]
Y <- keras::to_categorical(Y3)


Fe2 <- cmw_ran(Y, x, k = 2, maxit = 200, tol = 0.05, show_table = T)
Fe3 <- cmw_ran(Y, x, k = 3, maxit = 200, tol = 0.05, show_table = T)
Fe4 <- cmw_ran(Y, x, k = 4, maxit = 200, tol = 0.05, show_table = T)
Fe5 <- cmw_ran(Y, x, k = 5, maxit = 200, tol = 0.05, show_table = T)


dim(Y)
z <- cmd$V10
table(Fe3$classification)
Fe3$fr
cm <- confusionMatrix(sort(as.factor(Fe3$classification)), sort(as.factor(z)))
#cm <- confusionMatrix((as.factor(Fe2$classification-1)), (as.factor(y2)))
Fe3$Prob

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

## Adjusted Rand Index
clues::adjustedRand(sort(Fe2$classification-1),sort(y2))
clues::adjustedRand(sort(Fe3$classification),sort(z))
clues::adjustedRand(Fe4_1000$classification,z)
clues::adjustedRand(Fe5_1000$classification,z)

da <- data.frame(Y, x)
plot(da, col = sort(Fe3$classification), pch = sort(Fe3$classification) + 1, cex = .5)
plot(x, col = sort(z))

X <- as.data.frame(xt)
par(mfrow = c(1,2))
V <- lda(X, clust)$scaling
x.fda <- as.matrix(X) %*% V
plot(x.fda, type = "n", xlab = "", ylab = "", pch = 1, col = cluster, cex = .2)
legend(x = -6, y = 4, cex = 0.5, legend = c(3,5,8))
text(x.fda, col = cluster, labels = c(3,5,8)[cluster], cex = .2)


#r <- PCA(cmd[,-10], ncp = 5, graph = FALSE)
r <- PCA(tsne_out$Y, ncp = 5, graph = FALSE)

fviz_pca_ind(r, geom = "point", axes = c(1,2), 
             habillage = (as.factor(Fe3$classification)),
             palette = c("red", "blue", "gold","magenta", "black", "cyan","orange", 
                         "grey", "green","brown"),
             addEllipse = T, repel = F,
             ggtheme = theme_minimal(), title = "")

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

pred <- sort(as.factor(Fe2$classification))
colors = rainbow(length(unique(pred)))
names(colors) = unique(pred)
ecb = function(x){
  plot(x, t ='p', col=colors[e6]); 
  #text(x,labels=clust, col=colors[clust+1]) 
}

tsne_iris = tsne(x, epoch_callback = ecb, k = 6, perplexity = 20, max_iter = 1000)

rm(list = ls(all = TRUE))
# compare to PCA
dev.new()
pca_iris = princomp(xe$x)$scores[,1:2]
plot(pca_iris, t='n')
text(pca_iris, labels=pred,col=colors[pred])


colors <- rainbow(length(unique(pred)))
names(colors) <- unique(pred)

ecb = function(x){
  plot(x, t ='n');
  text(x, labels = pred, col = colors[pred], cex = .5)
}

tsne_iris = tsne(x2, k = 1, epoch_callback = ecb, perplexity=5, max_iter = 300)

x2 <- cbind(x, cluster)

qw <- colors()

table(clust)
w1 <- crep(w1[1], 1006)

cl <- fit$cluster
colors = rainbow(length(unique(pred)))
names(colors) = unique(pred)
# data.frame as input
cd <- as.data.frame(cmd[,-10])
tsne_out <- Rtsne(cmd[,-10], pca = F, theta=0.0, check_duplicates = FALSE)

#tol <- read.csv("tolu.csv")
# dim(tol)
# tol <- tol[,-1]
#text(tsne_out$Y, labels = pred, col = colors[pred], asp = 0, )

r <- PCA(tsne_out$Y, ncp = 5, graph = FALSE)

fviz_pca_ind(r, geom = "point", axes = c(1,3), 
             habillage = pred,
             addEllipse = FALSE, repel = TRUE,
             ggtheme = theme_minimal(), title = "",
             xlab = "", ylab = "")

############# Heart Data

heart <- read.table("D:\\traindcm\\Data\\Heart_disease\\processed.cleveland.data", sep = ",")

dim(heart)
head(heart)

heart[,12:13] <- with(heart, c(V12 <- as.numeric(levels(V12))[V12],
                               V13 <- as.numeric(levels(V13))[V13]))

## Replace all NA with 0
heart[which(is.na(heart$V12) == TRUE), 12] = 0
heart[which(is.na(heart$V13) == TRUE), 13] = 0

head(heart)
sum(is.na(heart$V14))

x <- as.matrix(heart)

cor(heart[,-14])
x <- x[, -14]

pairs(heart[,-14])

str(heart)

hist(heart)
?PCA
RT <- prcomp(x, scale = TRUE)
x <- RT$x[,1]
plot(RT, type = "l")
summary(RT)
# library(scatterplot3d)
# library(rgl)
# library(car)
# library(plot3D)

#data(VADeaths)
#  hist3D and ribbon3D with greyish background, rotated, rescaled,...
# hist3D(z = VADeaths, scale = FALSE, expand = 0.01, bty = "g", phi = 20,
#         col = "#0072B2", border = "black", shade = 0.2, ltheta = 90,
#         space = 0.3, ticktype = "detailed", d = 2)
# 
# 
# hist3D (x = 1:5, y = 1:4, z = VADeaths,
#         bty = "g", phi = 20,  theta = -60,
#         xlab = "", ylab = "", zlab = "", main = "VADeaths",
#         col = "#0072B2", border = "black", shade = 0.8,
#         ticktype = "detailed", space = 0.15, d = 2, cex.axis = 1e-9)
# # Use text3D to label x axis
#  text3D(x = 1:5, y = rep(0.5, 5), z = rep(3, 5),
#        labels = rownames(VADeaths),
#        add = TRUE, adj = 0)
# # Use text3D to label y axis
#  text3D(x = rep(1, 4),   y = 1:4, z = rep(0, 4),
#        labels  = colnames(VADeaths),
#        add = TRUE, adj = 1)

y1 <- heart[,14]

y2 <- ifelse(y1 == 0, 0, 1)
table(y2)

# One Hot encoding
Y <- to_categorical(y)




########### FOR MNIST not working yet 
cmw_ran_batch <- function(Y, x, k, Batch = 50, maxit = 1000, tol = 1e-10, show_table = FALSE){
  
  if(!is.data.frame(x)) x <- unname(as.matrix(x))
  if(!is.data.frame(Y)) Y <- unname(as.matrix(Y))
  
  eps = sqrt(.Machine$double.eps)
  d <- ncol(x); c <- ncol(Y) - 1
  L <- L1 <- L2 <- ai <- ai3 <- lb <- NULL
  
  #n <- nrow(x)
  a <- pr <- matrix(NA, nrow = Batch, ncol = c)
  py <- matrix(NA, nrow = Batch, ncol = c+1)
  dm <- matrix(NA, nrow = Batch, ncol = k)
  
  pu <- NULL
  sig <- array(NA, c(d,d,Batch))
  et1 <- list()
  
  w <- sigma <- list()
  m <- (k * c * (1 + d)) + (k * ((d*d - d)/2 + 2*d + 1))
  
  Pi <- rep(1/k, k)
  
  mean <- matrix(rnorm(d*k), nrow = k)
  
  for(l in 1:k){
    
    sigma[[l]] <- diag(d)
    weights <- matrix(runif(d*c), nrow = c, ncol = d)
    w[[l]] <- weights
    
  }
  
  count <- 2; 
  
  L <- c(-16000,-15000,-14000);
  
  repeat{
      
      set.seed(1)
      batchIndex <- sample(1:nrow(x), size = Batch, replace = T)
      batchx <- as.matrix(x[batchIndex, ])
      batchY <- Y[batchIndex, ]
      
      for(l in 1:k){
        
        for (i in 1:nrow(batchx)) {
          
          a[i,]   <- batchx[i,] %*% t(w[[l]])
          pr[i,]  <- exp(a[i,]) / (1 + sum(exp(a[i,])))
          pu[i]   <- 1 / (1 + sum(exp(a[i,]))); py[i,]  <- c(pu[i], pr[i,])
          dm[i,l] <- dmultinom(batchY[i,], sum(batchY[i,]), prob = py[i,]) *
            mvnfast::dmvn(batchx[i,], mean[l,], sigma[[l]], ncores = 8, 
                          isChol = F) * Pi[l]
          
        }
        
      }
      
      po <- dm / rowSums(dm)
      
      # M-Step
      
      # Mean for the group
      
      for(l in 1:k) mean[l, ] <- colSums(po[,l] * batchx) / colSums(po)[l]
      
      # Update Prior
      Pi <- colSums(po) / nrow(batchx)
      
      for (l in 1:k) {
        
        L2[l]  <- sum(po[,l] * mvnfast::dmvn(batchx, mean[l,], sigma[[l]], ncores = 8, 
                                             isChol = F, log = TRUE))
      }
      
    # if(count > 10){
    #   
    #   for(l in 1:k){
    #     
    #     et <- nnet::multinom(batchY ~ batchx, weights = po[,l], trace = FALSE)
    #     w[[l]] <- summary(et)$coefficients[,-1]
    #     lb[l] <- logLik(et)
    #     et1[[l]] <- summary(et)
    #     
    #   }
    #   
    # }
      
      #if(count > 20){
      
      for (l in 1:k) {
        
        for (i in 1:nrow(batchx)) {
          
          sig[,,i] <- po[i,l] * outer((batchx[i,] - mean[l,]),(batchx[i,] - mean[l,]))
          
        }
        
        sigma[[l]] <- apply(sig, c(1, 2), sum)/ sum(po[,l]) + diag(eps, ncol(batchx))
        
        # which(is.nan(sigma[[l]]))
        # 
        # sigma[[l]][which(is.nan(sigma[[l]]))] <- 1
      }
      
      #}
      
      # Compute Log_likelihood
      lik  <- sum(log(Pi) * po);lik2 <- sum(L2);lik3 <- sum(lb)
      
      # Compute Log_likelihood lki, lk2, 
      L[count] <- lik + lik2 + lik3
      
      a_k <- (L[count+1] - L[count]) / (L[count] - L[count-1])
      
      L[count + 2] <- L[count] + ((1-a_k)^-1 * (L[count+1] - L[count]))
      
      if (show_table) {
        
        dif <- abs(L[count+2] - L[count+1])
        
        out_table = data.frame(Iteration = count, Likelihood = L[count+2], difference = dif)
        print(kable(out_table))  
        
        if (dif < tol || count == maxit) break; 
        
      }
      
      count <- count + 1
      
    }
  
  Z <- apply(po, 1, which.max)
  
  Prof <- Poc1 <- po
  
  for(i in 1:nrow(batchx)){
    
    for(j in 1:k){
      
      Prof[i,j] <- ifelse(po[i,j] > 0.9, 1, 0)
      
      if(po[i,j] > 0){
        
        Poc1[i,j] <- log(po[i,j]) 
        
      }
      
    }
    
  }
  
  ### Information criterion
  
  # ai   <- -2*L[count + 2] + 2*m;
  # ai3  <- -2*L[count + 2] + m*log(n)
  # AIcc <-  ai - 2*m*(m+1)/(n-m-1)
  # AIC3 <- -2*L[count+2] - 3*m
  # AICu <-  AIcc - n*log(n/(n-m-1))
  # ICL  <-  ai3 + suppressWarnings(sum(rowSums(Prof * Poc1)))
  # Caic <- -2*L[count+2] - m*(1+log(n))
  # AWE  <- -2*L[count+2] - 2*m*(3/2 + log(n))
  
  return(list("mean" = mean, "Prob" = Pi, "post" = po,"weights" = w, "sigma" = sigma,
              "classification" = Z))
  
}

x <- as.matrix(x)

Fe5 <- cmw_ran_batch(Y, x, k = 10, Batch = 100, maxit = 100, tol = 0.05, show_table = T)
table(trainy)
>>>>>>> 181bc0ae3e45bb35c52f920c36b0bad947580001:200_cate.R
