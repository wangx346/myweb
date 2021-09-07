rm(list = ls())
library(glmnet)
library(mvtnorm) 
library(QICD)
library(MASS)
library(parallel)
library(tcltk)
library(flare)
library(rmutil)
library(caret) 

Lp.norm <- function(v,p=2){
  if(p==0){
    return(sum(v!=0))
  } else if(is.infinite(p)){
    return(max(abs(v)))
  } else{
    return(sum(abs(v)^p)^(1/p))
  }
}

est_lambda<-function(X,alpha0=0.1,c, times=1e2){
  n=dim(X)[1]
  res=NULL
  for (i in 1:times){
    epsilon_rank=sample(1:n,n)
    xi=2*epsilon_rank-(n+1)
    S=(-2/n/(n-1))*(t(X)%*%xi)
    res=c(res,max(abs(S)))
  }
  return(quantile(res,1-alpha0)*c)
}

p_diff <- function(theta, lambda, a=3.7){
  less <- theta<=lambda
  y <- pmax(a*lambda-theta,0)
  res = lambda*less+y*(1-less)/(a-1)
  return(pmax(res,1e-3))
} 

predict.slim <- function (object, newdata){
  pred.n = nrow(newdata)
  lambda.n = object$nlambda
  intcpt = matrix(rep(object$intercept, pred.n), 
                  nrow = pred.n, ncol = lambda.n, byrow = T)
  Y.pred = newdata %*% object$beta + intcpt
  
  return(Y.pred)
}

cv.sqrtlasso <- function(x,y,nfold=5,lambda_list){
  n <- dim(x)[1]
  split <- createFolds(y, k = nfold)
  Error <- NULL
  for(i in 1:nfold){
    xtrain <- x[-split[[i]],]
    ytrain <- y[-split[[i]]]
    xtest <- x[split[[i]],]
    obj <- slim(xtrain, ytrain,verbose = FALSE,lambda = lambda_list)
    yhat <- predict(obj, xtest)
    ytest <- matrix(y[split[[i]]],ncol =length(lambda_list),nrow = length(split[[i]]) )
    Error <- rbind(Error, colMeans((yhat-ytest)^2))
  }
  mse <- colMeans(Error)
  lambda_min <- lambda_list[which.min(mse)]
  se <- sd(mse)
  ind <- which(mse<=min(mse)+se)
  lambda_1se <- max(lambda_list[ind])
  return(list(lambda_min=lambda_min,lambda_1se=lambda_1se))
}

HBIC.scad <- function(x,y,beta,const = 1){
  size = dim(x)
  n = size[1]
  p_n = size[2]
  C_n = log(log(n))/n
  df = sum(abs(beta) > 0)
  return(log(mean((y - cbind(1,x)%*% beta)^2)) + log(p_n) *df * C_n/const)
}

Scad <- function(x,y){
  t1 = proc.time()[3]
  obj_lasso <- cv.glmnet(x,y)
  beta_int <- as.vector(coef(obj_lasso, s= "lambda.min"))
  lambda.max <- obj_lasso$lambda[1]
  lambda.list <- lambda.max* 0.001^(1:100/100) 
  
  beta.mat <- NULL
  hbic <- NULL
  for(i in 1:100){
    penalty <- p_diff(beta_int[-1],lambda.list[i]) 
    obj <- glmnet(x,y,lambda = lambda.list[i],penalty.factor=penalty)
    beta <- as.vector(coef(obj))
    
    beta.mat <- rbind(beta.mat, beta)
    hbic <- c(hbic, HBIC.scad(x,y, as.vector(beta)))
  }
  lambda.idx <- which.min(hbic)
  best_lambda <- lambda.list[lambda.idx]
  beta_scad <- beta.mat[lambda.idx,]
  time_scad =  proc.time()[3]-t1
  
  return(list(beta_scad = beta_scad, time_scad = time_scad, best_lambda = best_lambda))
  
}
 
HBIC <- function (y, X, beta, C_n, const=6){
  size = dim(X)
  p_n = size[2]
  df = sum(abs(beta) > 1e-06)
  return(log(sum(checkloss(y - X %*% beta))) + log(p_n) * df * C_n/const)
}

hbic.rankscad <- function(x_diff, y_diff, id, n, beta_int, lambda_list){
  hbic <- NULL
  Beta <- NULL
  C_n <- log(log(n))/n
  newx <- x_diff[id,]
  newy <- y_diff[id]
  complete <- (length(id)==length(y_diff))
  
  for(i in 1:length(lambda_list)){
    penalty <- p_diff(beta_int,lambda_list[i])
    x_update <- t(apply(newx,1,function(x)x/penalty))
    obj <- QICD(newy,x_update,tau=0.5,lambda=length(newy)/2,funname="lasso",intercept = FALSE, maxout=5)
    beta <- as.vector(obj$beta_final)/penalty
    Beta = rbind(Beta, beta)
    hbic = c(hbic, HBIC(y_diff,x_diff, as.vector(beta), C_n))
  }
  index = which.min(hbic)
  return(list(best_lambda = lambda_list[index], beta_final = Beta[index,]))
}

Rankreg <- function(x,y,lam_lasso,eta_list,c=1.01,incomplete=TRUE,C=5){
  elm_diff <- function(index,v){
    return(v[index[1],]-v[index[2],])
  }
  t1 <- proc.time()[3]
  n <- length(y)  
  index_list <- combn(1:n,2)
  x_diff <- t(apply(index_list,2,elm_diff,v=x))
  y_diff <- apply(index_list,2,elm_diff,v=y)
  n_diff <- length(y_diff)
  xbar <- colMeans(x)
  ybar <- mean(y)
  
  if(incomplete){
    id <- (1:n_diff)[as.logical(rbinom(n_diff,1,C*n/n_diff))]
    newx <- x_diff[id,]
    newy <- y_diff[id]
    obj_RankLasso <- QICD(newy,newx,tau=0.5,lambda=lam_lasso*length(newy)/2,funname="lasso",intercept = FALSE, 
                          thresh=1e-3, maxout=1)
    beta_RL <- as.vector(obj_RankLasso$beta_final)
  } else{
    id <- 1:n_diff
    obj_RankLasso <- QICD(y_diff,x_diff,tau=0.5,lambda=lam_lasso*n_diff/2,funname="lasso",intercept = FALSE,
                          thresh = 1e-3, maxout = 1)
    beta_RL <- as.vector(obj_RankLasso$beta_final)
  }
  intercpt_RL <- ybar- crossprod(beta_RL,xbar)
  
  beta_RankLasso <- c(intercpt_RL,beta_RL)
  time_RankLasso <- proc.time()[3]-t1
  
  obj_Rankscad <- hbic.rankscad(x_diff,y_diff,id,n,beta_RL,eta_list)
  beta_Rs <- as.vector(obj_Rankscad$beta_final) 
  intercpt_Rs <- ybar-crossprod(beta_Rs,xbar)
  beta_Rankscad <- c(intercpt_Rs,beta_Rs)
  time_Rankscad <- proc.time()[3]-t1
  
  return(list(beta_RankLasso = beta_RankLasso, time_RankLasso = time_RankLasso,
              beta_Rankscad = beta_Rankscad, time_Rankscad = time_Rankscad, best_eta = obj_Rankscad$best_lambda))
}

sqrtLasso <- function(x,y){
  t1 = proc.time()[3]
  obj_sqrtLasso<- slim(x,y,verbose = FALSE,nlambda = 10)
  best_lam <- cv.sqrtlasso(x,y,lambda_list = obj_sqrtLasso$lambda)$lambda_min
  ind = which(obj_sqrtLasso$lambda==best_lam)
  beta_sqrtLasso <-c(obj_sqrtLasso$intercept[ind],obj_sqrtLasso$beta[,ind]) 
  time_sqrtLasso =  proc.time()[3]-t1
  return(list(beta_sqrtLasso = beta_sqrtLasso, time_sqrtLasso = time_sqrtLasso, best_lambda = best_lam))
  
}

Lasso <- function(x,y){
  t1 = proc.time()[3]
  obj_Lasso<-cv.glmnet(x, y)
  best_lam <- obj_Lasso$lambda.min
  beta_Lasso <- as.vector(coef(obj_Lasso,s="lambda.min"))
  time_Lasso =  proc.time()[3]-t1
  return(list(beta_Lasso = beta_Lasso, time_Lasso = time_Lasso, best_lambda = best_lam))
  
}

data.generate <- function(n,p,true.beta,rho,type){
  Sigma=if(type=="AR") rho^abs(outer(1:p,1:p,"-")) else (1-rho)*diag(p)+rho*matrix(1,p,p)
  X = rmvnorm(n,mean = rep(0,p), sigma = Sigma)
  true.y = X%*%true.beta 
  return(list(X=X,true.y=true.y))
}

error <- function(err_type,n){
  if(err_type=="N_5")
    return(rnorm(n,sd=0.5)) 
  else if(err_type=="N(0,1)")
    return(rnorm(n))
  else if(err_type=="N2")
    return(rnorm(n,sd=sqrt(2)))
  else if(err_type=="Mixture Normal"){
    a<-rbinom(n,1,0.95)
    return(a*rnorm(n,sd=0.1)+(1-a)*rnorm(n,sd=10))
  }
  #return(0.9*rnorm(n)+0.1*rnorm(n,sd=5))
  else if(err_type=="Laplace")
    return(rlaplace(n))
  else if(err_type=="t")
    return(rt(n,df=4)*sqrt(2))
  else if(err_type=="Cauchy")
    return(rcauchy(n))
  else stop("Error Type doesn't exist")
}


Simulation <- function(x,true.y,true.beta,err_type,times,n,p,c=1.01,incomplete=FALSE,C=5){
  beta_list_lasso=NULL
  time_lasso=NULL
  Lam_lasso=NULL
  
  beta_list_sqrtlasso=NULL
  time_sqrtlasso=NULL
  Lam_sqrtlasso=NULL  
  
  beta_list_scad=NULL
  time_scad=NULL
  Lam_scad=NULL  
  
  beta_list_Ranklasso=NULL
  t1 = proc.time()[3]
  Lam_Ranklasso=est_lambda(x,c=c)
  est_time = proc.time()[3]-t1
  time_Ranklasso=NULL
  
  beta_list_Rankscad=NULL
  time_Rankscad=NULL
  Lam_Rankscad=NULL
   
  
  eta.max = 1.8*sqrt(log(p)/n)
  eta_list = exp(seq(log(eta.max), log(0.3*eta.max),length = 10))#(1:10*0.12)*eta_base 
  
  for(i in 1:times){
    y <- true.y+error(err_type,n)
    obj_lasso <- Lasso(x,y)
    obj_sqrtlasso <- sqrtLasso(x,y)
    obj_Rankreg <- Rankreg(x,y,Lam_Ranklasso,eta_list,incomplete=incomplete,C=C)
    obj_scad <- Scad(x,y) 
    
    
    
    beta_list_lasso <- rbind(beta_list_lasso,obj_lasso$beta_Lasso)
    time_lasso <- c(time_lasso, obj_lasso$time_Lasso)
    Lam_lasso <- c(Lam_lasso, obj_lasso$best_lambda)
    
    beta_list_sqrtlasso <- rbind(beta_list_sqrtlasso,obj_sqrtlasso$beta_sqrtLasso)
    time_sqrtlasso <- c(time_sqrtlasso, obj_sqrtlasso$time_sqrtLasso)
    Lam_sqrtlasso <- c(Lam_sqrtlasso, obj_sqrtlasso$best_lambda)
    
    beta_list_scad <- rbind(beta_list_scad,obj_scad$beta_scad)
    time_scad <- c(time_scad, obj_scad$time_scad)
    Lam_scad <- c(Lam_scad, obj_scad$best_lambda) 
    
    beta_list_Ranklasso <- rbind(beta_list_Ranklasso,obj_Rankreg$beta_RankLasso)
    time_Ranklasso <- c(time_Ranklasso, est_time+obj_Rankreg$time_RankLasso)
    
    beta_list_Rankscad <- rbind(beta_list_Rankscad,obj_Rankreg$beta_Rankscad)
    time_Rankscad <- c(time_Rankscad, est_time+obj_Rankreg$time_Rankscad)
    Lam_Rankscad <- c(Lam_Rankscad, obj_Rankreg$best_eta)
  }
  return(list(beta_list_lasso=beta_list_lasso,time_lasso=time_lasso, Lam_lasso=Lam_lasso,
              beta_list_sqrtlasso=beta_list_sqrtlasso,time_sqrtlasso=time_sqrtlasso, Lam_sqrtlasso=Lam_sqrtlasso,
              beta_list_scad=beta_list_scad,time_scad=time_scad, Lam_scad=Lam_scad,
              beta_list_Ranklasso=beta_list_Ranklasso,time_Ranklasso=time_Ranklasso, Lam_Ranklasso=Lam_Ranklasso,
              beta_list_Rankscad=beta_list_Rankscad,time_Rankscad=time_Rankscad, Lam_Rankscad=Lam_Rankscad))
}


Estimation.evaluate <- function(sim.res,true.beta,rho,type){
  p <- length(true.beta)-1
  Sigma <- if(type=="AR") rho^abs(outer(1:p,1:p,"-")) else (1-rho)*diag(p)+rho*matrix(1,p,p)
  Sigma <- cbind(0,Sigma)
  Sigma <- rbind(c(1,rep(0,p)),Sigma)
  
  L1error_lasso <- apply(sim.res$beta_list_lasso,1,function(x)Lp.norm(x-true.beta,p=1))
  L2error_lasso <- apply(sim.res$beta_list_lasso,1,function(x)Lp.norm(x-true.beta))
  FN_lasso <- apply(sim.res$beta_list_lasso,1,function(x) sum(x[true.beta!=0]==0))
  FP_lasso <- apply(sim.res$beta_list_lasso,1,function(x) sum((x[true.beta==0]!=0)[-1]))
  prederr_lasso <- apply(sim.res$beta_list_lasso,1,function(x) (x-true.beta)%*%Sigma%*%(x-true.beta))
  
  L1error_sqrtlasso <- apply(sim.res$beta_list_sqrtlasso,1,function(x)Lp.norm(x-true.beta,p=1))
  L2error_sqrtlasso <- apply(sim.res$beta_list_sqrtlasso,1,function(x)Lp.norm(x-true.beta))
  FN_sqrtlasso <- apply(sim.res$beta_list_sqrtlasso,1,function(x) sum(x[true.beta!=0]==0))
  FP_sqrtlasso <- apply(sim.res$beta_list_sqrtlasso,1,function(x) sum((x[true.beta==0]!=0)[-1]))
  prederr_sqrtlasso <- apply(sim.res$beta_list_sqrtlasso,1,function(x) (x-true.beta)%*%Sigma%*%(x-true.beta))

  L1error_scad <- apply(sim.res$beta_list_scad,1,function(x)Lp.norm(x-true.beta,p=1))
  L2error_scad <- apply(sim.res$beta_list_scad,1,function(x)Lp.norm(x-true.beta))
  FN_scad <- apply(sim.res$beta_list_scad,1,function(x) sum(x[true.beta!=0]==0))
  FP_scad <- apply(sim.res$beta_list_scad,1,function(x) sum((x[true.beta==0]!=0)[-1]))
  prederr_scad <- apply(sim.res$beta_list_scad,1,function(x) (x-true.beta)%*%Sigma%*%(x-true.beta))
  
  L1error_Ranklasso <- apply(sim.res$beta_list_Ranklasso,1,function(x)Lp.norm(x-true.beta,p=1))
  L2error_Ranklasso <- apply(sim.res$beta_list_Ranklasso,1,function(x)Lp.norm(x-true.beta))
  FN_Ranklasso <- apply(sim.res$beta_list_Ranklasso,1,function(x) sum(x[true.beta!=0]==0))
  FP_Ranklasso <- apply(sim.res$beta_list_Ranklasso,1,function(x) sum((x[true.beta==0]!=0)[-1]))
  prederr_Ranklasso <- apply(sim.res$beta_list_Ranklasso,1,function(x) (x-true.beta)%*%Sigma%*%(x-true.beta))
  
  L1error_Rankscad <- apply(sim.res$beta_list_Rankscad,1,function(x)Lp.norm(x-true.beta,p=1))
  L2error_Rankscad <- apply(sim.res$beta_list_Rankscad,1,function(x)Lp.norm(x-true.beta))
  FN_Rankscad <- apply(sim.res$beta_list_Rankscad,1,function(x) sum(x[true.beta!=0]==0))
  FP_Rankscad <- apply(sim.res$beta_list_Rankscad,1,function(x) sum((x[true.beta==0]!=0)[-1]))
  prederr_Rankscad <- apply(sim.res$beta_list_Rankscad,1,function(x) (x-true.beta)%*%Sigma%*%(x-true.beta))
  
  
  pure_mean<-function(x){
    #return(round(mean(x),2))
    inter_dist=IQR(x)
    lower=quantile(x,0.25)-1.5*inter_dist
    upper=quantile(x,0.75)+1.5*inter_dist
    return(round(mean(x[x<=upper & x>=lower]),2))
  }
  pure_sd<-function(x){
    #return(round(sd(x),2))
    inter_dist=IQR(x)
    lower=quantile(x,0.25)-1.5*inter_dist
    upper=quantile(x,0.75)+1.5*inter_dist
    return(round(sd(x[x<=upper & x>=lower])/sqrt(length(x)),2))
  }
  
  res.latex <- function(res){
    res.lat <- res[1] 
    for(i in 1:5){
      char_new <- paste0(res[2*i]," (",res[2*i+1], ")")
      res.lat <- c(res.lat,char_new)
    }
    res.lat <- c(res.lat,res[12])
    names(res.lat) <- c("lambda", "l1 error", "l2 error", 
                        "pred error", "FP", "FN", "time")
    return(res.lat)
  }
  
  Lasso <- c(lambda = pure_mean(sim.res$Lam_lasso),
             ave_L1err=pure_mean(L1error_lasso),
             se_L1err=pure_sd(L1error_lasso),
             ave_L2err=pure_mean(L2error_lasso),
             se_L2err=pure_sd(L2error_lasso),
             ave_prederr=pure_mean(prederr_lasso),
             se_prederr=pure_sd(prederr_lasso),
             ave_FP=pure_mean(FP_lasso),
             se_FP=pure_sd(FP_lasso),
             ave_FN=pure_mean(FN_lasso),
             se_FN=pure_sd(FN_lasso),
             time = pure_mean(sim.res$time_lasso))
  
  sqrtLasso <- c(lambda = pure_mean(sim.res$Lam_sqrtlasso),
                 ave_L1err=pure_mean(L1error_sqrtlasso),
                 se_L1err=pure_sd(L1error_sqrtlasso),
                 ave_L2err=pure_mean(L2error_sqrtlasso),
                 se_L2err=pure_sd(L2error_sqrtlasso),
                 ave_prederr=pure_mean(prederr_sqrtlasso),
                 se_prederr=pure_sd(prederr_sqrtlasso),
                 ave_FP=pure_mean(FP_sqrtlasso),
                 se_FP=pure_sd(FP_sqrtlasso),
                 ave_FN=pure_mean(FN_sqrtlasso),
                 se_FN=pure_sd(FN_sqrtlasso),
                 time = pure_mean(sim.res$time_sqrtlasso))
  
  SCAD <- c(lambda = pure_mean(sim.res$Lam_scad),
                 ave_L1err=pure_mean(L1error_scad),
                 se_L1err=pure_sd(L1error_scad),
                 ave_L2err=pure_mean(L2error_scad),
                 se_L2err=pure_sd(L2error_scad),
                 ave_prederr=pure_mean(prederr_scad),
                 se_prederr=pure_sd(prederr_scad),
                 ave_FP=pure_mean(FP_scad),
                 se_FP=pure_sd(FP_scad),
                 ave_FN=pure_mean(FN_scad),
                 se_FN=pure_sd(FN_scad),
                 time = pure_mean(sim.res$time_scad))
  
  RankLasso <- c(lambda = pure_mean(sim.res$Lam_Ranklasso),
                 ave_L1err=pure_mean(L1error_Ranklasso),
                 se_L1err=pure_sd(L1error_Ranklasso),
                 ave_L2err=pure_mean(L2error_Ranklasso),
                 se_L2err=pure_sd(L2error_Ranklasso),
                 ave_prederr=pure_mean(prederr_Ranklasso),
                 se_prederr=pure_sd(prederr_Ranklasso),
                 ave_FP=pure_mean(FP_Ranklasso),
                 se_FP=pure_sd(FP_Ranklasso),
                 ave_FN=pure_mean(FN_Ranklasso),
                 se_FN=pure_sd(FN_Ranklasso),
                 time = pure_mean(sim.res$time_Ranklasso))
  
  Rankscad <- c(lambda = pure_mean(sim.res$Lam_Rankscad),
                ave_L1err=pure_mean(L1error_Rankscad),
                se_L1err=pure_sd(L1error_Rankscad),
                ave_L2err=pure_mean(L2error_Rankscad),
                se_L2err=pure_sd(L2error_Rankscad),
                ave_prederr=pure_mean(prederr_Rankscad),
                se_prederr=pure_sd(prederr_Rankscad),
                ave_FP=pure_mean(FP_Rankscad),
                se_FP=pure_sd(FP_Rankscad),
                ave_FN=pure_mean(FN_Rankscad),
                se_FN=pure_sd(FN_Rankscad),
                time = pure_mean(sim.res$time_Rankscad))
  
  Lasso.latex <- res.latex(Lasso) 
  sqrtLasso.latex <- res.latex(sqrtLasso) 
  SCAD.latex <- res.latex(SCAD) 
  RankLasso.latex <- res.latex(RankLasso) 
  Rankscad.latex <- res.latex(Rankscad) 
  
  
  final.table <- rbind(Lasso=Lasso.latex,sqrtLasso=sqrtLasso.latex,SCAD=SCAD.latex,
                       RankLasso=RankLasso.latex,Rankscad=Rankscad.latex)
  print(xtable(final.table, type = "latex"), file = filename, append = T) 
  
  
  return(final.table)
}



#################### Example 1: t_4 errors ####################

k <- 3; n <- 100; p <- 400
times <- 200; mc.cores <- 4
 
if(k==3){
  true.beta <- c(rep(sqrt(3),k),rep(0,p-k))
} else if(k==7){
  true.beta <- c(2,1.5,1.25,1,0.75,0.5,0.25,rep(0,p-k))  
} else if(k==15){
  true.beta <- c(2,2,2,1.5,1.5,1.25,1.25,1,1,0.75,0.75,0.5,0.5,0.25,0.25,rep(0,p-k))
} else{
  true.beta <- c(2,2,2,2,1.75,1.75,1.75,1.5,1.5,1.5,1.25,1.25,1.25,1,1,1,0.75,0.75,0.75,0.5,0.5,0.5,0.25,0.25,0.25,rep(0,p-k))
}


set.seed(1234)
name <- paste("Sim_res_k",k,sep="")
filename <- paste(name,".txt",sep="")
write(name, file = filename)
Error_type = "t"

data = data.generate(n,p,true.beta,rho=0.5,type="CS")
write(paste("\nError Type: ",Error_type,sep=""), file = filename, append = T)
sim.res.par <- mclapply(rep(as.integer(times/mc.cores),mc.cores), Simulation, x=data$X,true.y=data$true.y,
                        true.beta=true.beta,err_type=Error_type,n=n,p=p,mc.cores=mc.cores)
write("sim.res.par Complete!", file = filename, append = T)

sim.res.par1 <- do.call(rbind,sim.res.par)
varname <- colnames(sim.res.par1)
for(j in 1:length(varname)){
  assign(varname[j], do.call(rbind,sim.res.par1[,j]))
}
sim.res <- list(beta_list_lasso=beta_list_lasso,time_lasso=time_lasso, Lam_lasso=Lam_lasso,
                beta_list_sqrtlasso=beta_list_sqrtlasso,time_sqrtlasso=time_sqrtlasso, Lam_sqrtlasso=Lam_sqrtlasso,
                beta_list_scad=beta_list_scad,time_scad=time_scad, Lam_scad=Lam_scad,
                beta_list_Ranklasso=beta_list_Ranklasso,time_Ranklasso=time_Ranklasso, Lam_Ranklasso=Lam_Ranklasso,
                beta_list_Rankscad=beta_list_Rankscad,time_Rankscad=time_Rankscad, Lam_Rankscad=Lam_Rankscad) 
est.eva <- Estimation.evaluate(sim.res,c(0,true.beta),rho=0.5,type="CS") 

read_csv(filename)
