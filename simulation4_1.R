rm(list = ls())
library(glmnet)
library(tilting)
library("corpcor")
library(Rfast)
library("latex2exp")
library(sail)
library(caret)

# Simulation corresponding to figures 4,5

N = 2000 # Number of independent simulations


beta_naive_error = c()
beta_spectral_error = c()

# vector of simulated \hat{\beta}_1 using traditional debiasing:
beta_1_naive_DD_vec = c()
# vector of simulated \hat{\beta}_{s+1} using traditional debiasing:
beta_2_naive_DD_vec = c()

# vector of simulated \hat{\beta}_1^{DD} using doubly debiasing:
beta_1_DD_vec = c()
# vector of simulated \hat{\beta}_{s+1}^{DD} using doubly debiasing:
beta_2_DD_vec = c()

# vector of simulated std error of \hat{\beta}_1 using traditional debiasing:
est_sd_naive_1_vec = c()
# vector of simulated std error of \hat{\beta}_{s+1} using traditional debiasing:
est_sd_naive_2_vec = c()

# vector of simulated std error of \hat{\beta}_1^{DD} using doubly debiasing:
est_sd_1_vec = c()
# vector of simulated std error of \hat{\beta}_{j+1}^{DD} using doubly debiasing:
est_sd_2_vec = c()

# vector of simulated std error of \hat{\beta}_1^{DD} using doubly debiasing 
# and alternative estimator:
est_sd_1_fold_vec = c()
# vector of simulated std error of \hat{\beta}_{s+1}^{DD} using doubly debiasing 
# and alternative estimator:
est_sd_2_fold_vec = c()
# vector of simulated \hat{\sigma}_{\nu}^{fold} according to estimator in Guo (2020):
sigma_e_est = c()
# vector of simulated \hat{\sigma}_{\nu}^{fold} according to proposed method:
sigma_e_est_fold = c()

set.seed(2024)
n = 200
p = 1000
q = 30 
s = 5 
beta = matrix(c(rep(1,s), rep(0, p-s)), nrow=p)
sigma_E = 1
delta = rnorm(q, mean=0, sd=1)

Gamma = matrix(rnorm(q*p, mean=0, sd=1), nrow=q,ncol=p)
sigma_nu = 1
nu =  rnorm(n, mean=0, sd=sigma_nu) 

require(doMC)
registerDoMC(cores = 4)

# How much we enlarge the penalty term chosen by CV:
lambda_inflate = 1.2


for (r in 1:N) {
  print(r)
  H = matrix(rnorm(n*q, mean=0, sd=1),nrow=n, ncol=q)
  E = matrix(rnorm(p*n,mean=0, sd = sigma_E), nrow=n)
  X = H%*% Gamma + E 
  Y = X%*% beta + H %*% delta + nu 
  
  # Naive lasso regression:
  cvfit = cv.glmnet(y=Y,x=X, penalty.factor = col.norm(X),parallel = T,intercept=F)
  beta.hat.naive = as.vector(coef(cvfit,s= lambda_inflate*cvfit$lambda.min))[-1]

  # Trim transform:
  svd_res = fast.svd(X)
  d = svd_res$d
  tau = median(d)
  U = svd_res$u
  d_tilde = pmin(d,tau)
  
  F_matrix = U[,]
  for (k in 1:dim(U)[1]) {
    F_matrix[,k] = d_tilde[k]/d[k] * F_matrix[,k]
  }
  F_matrix = Tcrossprod(F_matrix,U)
  
  X_tilde = F_matrix %*% X
  Y_tilde = F_matrix %*% Y
  
  
  # Spectral decofounded lasso
  cvfit = cv.glmnet(y=Y_tilde,x=X_tilde,penalty.factor = col.norm(X_tilde), parallel = T,intercept=F,keep=T)
  beta.hat = as.vector(coef(cvfit,s= lambda_inflate*cvfit$lambda.min))[-1]
  
  beta_spectral_error= c(beta_spectral_error, sum(abs(beta.hat - beta)) )
  beta_naive_error = c(beta_naive_error,sum(abs(beta.hat.naive - beta))) 
  
  # Estimating variance of \nu
  trQ2 =sum(diag(F_matrix%*%F_matrix))
  sigma_e_hat =1/sqrt(trQ2)  *norm(Y_tilde -X_tilde%*%beta.hat,type="2")
  
  # Alterative estimate of varriance of \nu
  sigma_e_folds = norm(Y_tilde -cvfit$fit.preval[,which(cvfit$lambda== cvfit$lambda.min)],type="2") / sqrt(trQ2)
  
  # saving for book keeping
  sigma_e_est = c(sigma_e_est, sigma_e_hat)
  sigma_e_est_fold = c(sigma_e_est_fold, sigma_e_folds)
  
  # Naive debiasing
  
  # Inference on beta_1
  naiveinternallasso1CV = cv.glmnet(y=X[,1],x=X[,-1],penalty.factor = col.norm(X[,-1]),parallel = T,intercept=F)
  naiveinternallasso1_coeff = as.vector(coef(naiveinternallasso1CV,s= naiveinternallasso1CV$lambda.min))
  naive_residual1 = (X[,1]) - X[,-1]%*%naiveinternallasso1_coeff[-1]
  PjZjTPjX = crossprod(naive_residual1,X[,1])
  beta_1_naive_DD = crossprod(naive_residual1,Y - X[,-1]%*%beta.hat.naive[-1]) /(PjZjTPjX)
  beta_1_naive_DD_vec = c(beta_1_naive_DD_vec, beta_1_naive_DD)
  sigma_naive_hat_1 = sqrt(1/n) *norm(Y-X[,-1]%*%beta.hat.naive[-1],type="2")
  est_sd_naive_1 = sigma_naive_hat_1 *norm(naive_residual1, type="2") / (PjZjTPjX)
  est_sd_naive_1
  est_sd_naive_1_vec = c(est_sd_naive_1_vec, est_sd_naive_1)
  
  
  # Inference on beta_1 with doubly debiasing
  svd_1 = fast.svd(X[,-1])
  d1 = svd_1$d
  tau1 = median(d1)
  U1 = svd_1$u
  d_tilde_1 = pmin(d1,tau1)
  
  F1 = U1[,]
  for (i in 1:dim(U1)[1]) {
    F1[,i] = d_tilde_1[i]/d1[i] * F1[,i]
  }
  F1=  Tcrossprod(F1,U1)
  
  X_tilde_1 = F1 %*% X[,-1]
  Y_tilde_1 = F1 %*% X[,1]
  
  internallasso1CV = cv.glmnet(y=Y_tilde_1,x=X_tilde_1,penalty.factor = col.norm(X_tilde_1),parallel = T,intercept=F)
  internallasso1_coeff = as.vector(coef(internallasso1CV,s= internallasso1CV$lambda.min))
  residual1 = (X[,1]) - X[,-1]%*%internallasso1_coeff[-1]
  PjZ = F1%*%residual1
  beta_1_DD = crossprod(PjZ,F1%*%(Y- X[,-1]%*%beta.hat[-1])) /crossprod(PjZ,F1%*%X[,1])
  beta_1_DD_vec = c(beta_1_DD_vec, beta_1_DD)
  
  P2 = F1%*%F1
  P4 = P2%*%P2
  est_sd_1 = sigma_e_hat *sqrt(crossprod(residual1, P4%*%residual1)) / abs(crossprod(residual1,P2%*%X[,1]))
  est_sd_1_vec = c(est_sd_1_vec, est_sd_1)
  
  # Book keeepig
  est_sd_1_fold = sigma_e_folds /sigma_e_hat * est_sd_1
  est_sd_1_fold_vec = c(est_sd_1_fold_vec, est_sd_1_fold)
  
  
}
# computing var of observed error term (\sigma and not \sigma_{\nu})
Sigma = Crossprod(Gamma,Gamma)
for (k in 1:p) {
  Sigma[k,k] = Sigma[k,k]+sigma_E
}
SigmaInv= spdinv(Sigma)
sigma2 = sigma_nu^2 + t(delta)%*%(diag(q)-Gamma %*%SigmaInv%*%t(Gamma))%*%delta


# Saving plots:

# Histogram of sigma_hat
dev.new(width=10,height=5,noRStudioGD = TRUE)
par(mfrow=c(1,2),mar=c(4,4,4,4))
hist(sigma_e_est,main=TeX("Histogram of $\\hat{\\sigma}$"), probability = TRUE,ylab="Density",xlab=TeX("$x$"),xlim=c(min(c(sigma_e_est,sigma_e_est_fold))-0.1,max(c(sigma_e_est,sigma_e_est_fold))))
abline(v = sigma_nu)
abline(v = sqrt(sigma2),lty=2)
hist(sigma_e_est_fold, main=TeX("Histogram of $\\hat{\\sigma}_{fold}$"),probability = TRUE,ylab="Density",xlab=TeX("$x$"),xlim=c(min(c(sigma_e_est,sigma_e_est_fold))-0.1,max(c(sigma_e_est,sigma_e_est_fold))))
abline(v = sigma_nu)
abline(v = sqrt(sigma2),lty=2)
quartz.save(file="/Users/peraugust/Documents/cambridge/essay/simulations/simulation4_1.jpeg",type="jpeg",dpi=300)


# Histogram of t statistics
dev.new(width=10,height=10,noRStudioGD = TRUE)
par(mfrow=c(2,1),mar=c(4,4,4,4))

max_x = max(abs(c((beta_1_DD_vec-beta[1])/est_sd_1_vec,(beta_1_DD_vec-beta[1])/est_sd_1_fold_vec)))
#hist((beta_1_DD_vec), breaks=25, probability =T)
hist((beta_1_DD_vec-beta[1])/est_sd_1_vec, breaks=100, probability = T, ylim=c(0,0.45),xlim=c(-max_x,max_x),
     main = TeX("Histogram of $(\\hat{\\beta}_1- \\beta_1)/ \\hat{V}$"),ylab="Density",xlab="x")
dens1 = density((beta_1_DD_vec-beta[1])/est_sd_1_vec)
lines(dens1$x, dens1$y,col=2,lwd=2)
xs = seq(-5,5, length.out=1000)
lines(xs, dnorm(xs),lwd=2)
legend(x="topright", legend=c(TeX("N$(0,1)$ density"), "KDE"),lty=c(1,1),col=c(1,2),lwd=c(2,2))

hist((beta_1_DD_vec-beta[1])/est_sd_1_fold_vec, breaks=40, probability = T,ylim=c(0,0.45),xlim=c(-max_x,max_x),
     main = TeX("Histogram of $(\\hat{\\beta}_1- \\beta_1)/ \\hat{V}_{fold}$"),ylab="Density",xlab="x")
dens2 = density((beta_1_DD_vec-beta[1])/est_sd_1_fold_vec)
lines(dens2$x, dens2$y,col=2,lwd=2)
xs = seq(-5,5, length.out=1000)
lines(xs, dnorm(xs),lwd=2)
legend(x="topright", legend=c(TeX("N$(0,1)$ density"), "KDE"),lty=c(1,1),col=c(1,2),lwd=c(2,2))

quartz.save(file="/Users/peraugust/Documents/cambridge/essay/simulations/simulation4_2.jpeg",type="jpeg",dpi=300)


#liz=list(beta,beta_1_DD_vec,beta_1_naive_DD_vec,est_sd_1_vec,est_sd_1_fold_vec,sigma_e_est,sigma_e_est_fold)
#save(liz,file="/Users/peraugust/Documents/cambridge/essay/simulations/simulation4_1_data.dat")


# Other diagostics not included in essay

mean(beta_spectral_error)
mean(beta_naive_error)
hist((beta_1_DD_vec), breaks=25, probability =T)
hist((beta_1_DD_vec-beta[1])/est_sd_1_vec, breaks=100, probability = T)
xs = seq(-5,5, length.out=1000)
lines(xs, dnorm(xs))
hist((beta_1_DD_vec-beta[1])/est_sd_1_fold_vec, breaks=100, probability = T)
xs = seq(-5,5, length.out=1000)
lines(xs, dnorm(xs))

hist((beta_1_DD_vec-beta[1])/est_sd_1_vec, breaks=100, probability = T)
xs = seq(-5,5, length.out=1000)
lines(xs, dnorm(xs))

#Comparison of estimates of sigma_{\nu}
hist(sigma_e_est_fold)
hist(sigma_e_est)


# Checking coverage of Doubly Debiased confidence interval (for \beta_1)
CI1_left = beta_1_DD_vec -qnorm(0.975)*est_sd_1_vec
CI1_right = beta_1_DD_vec +qnorm(0.975)*est_sd_1_vec
coverage = 0
for (i in 1:length(CI1_left)) {
  if(beta[1] <=CI1_right[i]){
    if(beta[1] >=CI1_left[i]){
      coverage = coverage+1
    }
  }
}

coverage = coverage / length(CI1_left)


# Checking coverage of Doubly Debiased confidence interval (for \beta_1) with 
# alternative estimator of \sigma_{\nu}
CI2_left = beta_1_DD_vec -qnorm(0.975)*est_sd_1_fold_vec
CI2_right = beta_1_DD_vec +qnorm(0.975)*est_sd_1_fold_vec
coverage_fold= 0
for (i in 1:length(CI2_left)) {
  if(beta[1] <=CI2_right[i]){
    if(beta[1] >=CI2_left[i]){
      coverage_fold = coverage_fold+1
    }
  }
}

coverage_fold = coverage_fold / length(CI2_left)
coverage_fold


