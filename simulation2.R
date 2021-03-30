library(glmnet)
library("corpcor")
library(Rfast)
library("latex2exp")

### Simulation corresponding to Figure 4


# Case 1: q small
set.seed(2023)

# Model parameters
p = 600
q = 6
s = 5 
n = 300
sigma_E = 2
beta = matrix(c(rep(1,s), rep(0, p-s)), nrow=p) 
sigma_nu = 1

# Simulation parameters
cv_scale_values = exp(-seq(-5,1, by=0.1))
N = 2000

# Book keeping 
lasso.error.naive = matrix(data=NA, ncol = length(cv_scale_values), nrow = N)
lasso.error = matrix(data=NA, ncol = length(cv_scale_values), nrow = N)
lambda.min.chosen.naive = c()
lambda.min.chosen = c()


for (j in 1:N) {
  if(j%%10 == 0){
    print("j=")
    print(j)
  }

  # Generating delta and gamma
  delta = rnorm(q, mean=0, sd=1)# 
  Gamma = matrix(rnorm(q*p, mean=0, sd=1), nrow=q,ncol=p)
  
  # Generating data
  nu =  rnorm(n, mean=0, sd=sigma_nu) 
  H = matrix(rnorm(n*q, mean=0, sd=1),nrow=n, ncol=q)
  E = matrix(rnorm(p*n,mean=0, sd = sigma_E), nrow=n,ncol=p)
  X = H%*% Gamma + E 
  Y = X%*% beta + H %*% delta + nu 
  
  # Standard/Naive lasso
  cv.fit.naive = cv.glmnet(y=Y,x=X, intercept = FALSE, standardize=FALSE)
  lambdas = cv_scale_values*(cv.fit.naive$lambda.min)
  cv.fit.naive.multiple = glmnet(y=Y,x=X, intercept = FALSE, standardize=FALSE,lambda = lambdas)
  beta.hat.naive = coef(cv.fit.naive.multiple)[-1,]
  lasso.error.naive[j,] = apply(beta.hat.naive, 2,function(x) sum(abs(x - beta)))
  lambda.min.chosen.naive = c(lambda.min.chosen.naive, cv.fit.naive$lambda.min)
  
  # Trimmed Lasso as in Cevid (2018):
  
  # TRIM transform: 
  svd_X = fast.svd(X)
  d = svd_X$d
  tau = median(d)
  U = svd_X$u
  d_tilde = pmin(d,tau)
  
  F_matrix = U[,]
  for (k in 1:dim(U)[1]) {
    F_matrix[,k] = d_tilde[k]/d[k] *F_matrix[,k]
  }
  F_matrix = Tcrossprod(F_matrix,U)

  X_tilde = mat.mult(F_matrix,X)
  Y_tilde = F_matrix %*% Y
  
  
  cv.fit = cv.glmnet(y=Y_tilde,x=X_tilde,intercept = FALSE, standardize=FALSE)
  lambdas = cv_scale_values*(cv.fit$lambda.min)
  cv.fit.multiple = glmnet(y=Y_tilde,x=X_tilde, intercept = FALSE, standardize=FALSE,lambda = lambdas)
  beta.hat = coef(cv.fit.multiple)[-1,]
  lasso.error[j,] = apply(beta.hat, 2,function(x) sum(abs(x - beta)))
  lambda.min.chosen = c(lambda.min.chosen, cv.fit$lambda.min)
}


# Computing the L1 norm of the bias b
set.seed(2023)
b.norm = c()
for (i in 1:N) {
  delta = rnorm(q, mean=0, sd=1)# len q
  Gamma = matrix(rnorm(q*p, mean=0, sd=1), nrow=q,ncol=p)
  Sigma = Crossprod(Gamma,Gamma)
  for (k in 1:p) {
    Sigma[k,k] = Sigma[k,k]+1
  }
  SigmaInv= spdinv(Sigma)
  b = Tcrossprod(SigmaInv, Gamma) %*%delta
  b.norm = c(b.norm, norm(b, type="1"))
}
b.norm.mean =mean(b.norm)






# Case 2: q large
set.seed(2023)

# Model prameters
q = 100

# Book keeping
lasso.error.naive2 = matrix(data=NA, ncol = length(cv_scale_values), nrow = N)
lasso.error2 = matrix(data=NA, ncol = length(cv_scale_values), nrow = N)
lambda.min.chosen.naive2 = c()
lambda.min.chosen2 = c()


for (j in 1:N) {
  if(j%%10 == 0){
    print("j=")
    print(j)
  }
  # Generating delta and gamma
  delta = rnorm(q, mean=0, sd=1)
  Gamma = matrix(rnorm(q*p, mean=0, sd=1), nrow=q,ncol=p)
  
  # Generating data
  nu =  rnorm(n, mean=0, sd=sigma_nu) 
  H = matrix(rnorm(n*q, mean=0, sd=1),nrow=n, ncol=q)
  E = matrix(rnorm(p*n,mean=0, sd = sigma_E), nrow=n,ncol=p)
  X = H%*% Gamma + E 
  Y = X%*% beta + H %*% delta + nu 
  
  # Standard/Naive lasso
  cv.fit.naive = cv.glmnet(y=Y,x=X, intercept = FALSE, standardize=FALSE)
  lambdas = cv_scale_values*(cv.fit.naive$lambda.min)
  cv.fit.naive.multiple = glmnet(y=Y,x=X, intercept = FALSE, standardize=FALSE,lambda = lambdas)
  beta.hat.naive = coef(cv.fit.naive.multiple)[-1,]
  lasso.error.naive2[j,] = apply(beta.hat.naive, 2,function(x) sum(abs(x - beta)))
  lambda.min.chosen.naive2 = c(lambda.min.chosen.naive2, cv.fit.naive$lambda.min)
  
  # Trimmed Lasso as in Cevid (2018)
  
  # TRIM transform: 
  svd_X = fast.svd(X)
  d = svd_X$d
  tau = median(d)
  U = svd_X$u
  d_tilde = pmin(d,tau)
  
  F_matrix = U[,]
  for (k in 1:dim(U)[1]) {
    F_matrix[,k] = d_tilde[k]/d[k] *F_matrix[,k]
  }
  F_matrix = Tcrossprod(F_matrix,U)
  
  X_tilde = mat.mult(F_matrix,X)
  Y_tilde = F_matrix %*% Y
  
  
  cv.fit = cv.glmnet(y=Y_tilde,x=X_tilde,intercept = FALSE, standardize=FALSE)
  lambdas = cv_scale_values*(cv.fit$lambda.min)
  cv.fit.multiple = glmnet(y=Y_tilde,x=X_tilde, intercept = FALSE, standardize=FALSE,lambda = lambdas)
  beta.hat = coef(cv.fit.multiple)[-1,]
  lasso.error2[j,] = apply(beta.hat, 2,function(x) sum(abs(x - beta)))
  lambda.min.chosen2 = c(lambda.min.chosen2, cv.fit$lambda.min)
}


# Computing the L1 norm of the bias b
set.seed(2023)
b.norm = c()
for (i in 1:N) {
  delta = rnorm(q, mean=0, sd=1)# len q
  Gamma = matrix(rnorm(q*p, mean=0, sd=1), nrow=q,ncol=p)
  Sigma = Crossprod(Gamma,Gamma)
  for (k in 1:p) {
    Sigma[k,k] = Sigma[k,k]+sigma_E
  }
  SigmaInv= spdinv(Sigma)
  b = Tcrossprod(SigmaInv, Gamma) %*%delta
  b.norm = c(b.norm, norm(b, type="1"))
}

b.norm.mean2 =mean(b.norm)




# Plotting
directory = "INSERT HERE"
filename = "simulation2.pdf"
setwd(directory)
pdf(file = filename, width = 10, height = 5)
par(mfrow=c(1,2),mar=c(4,4,4,4))

plot(y=colMeans(lasso.error.naive), x=log(cv_scale_values), ylab ="L1 error", main = TeX("$p=600,n=300,q=6,s=5,\\sigma_{\\nu}=1"), type="l", ylim=c(-0.5,6),xlab=TeX("log$(\\phi)$"))
lines(y=colMeans(lasso.error),x= log(cv_scale_values), col=2, type="l")
abline(h = b.norm.mean,col=3,lty=2)
abline(v = 0,col=1,lty=2)
legend(x="bottomright", legend=c("Standard Lasso", "Trimmed Lasso",TeX("$|| b||_1$")),lty=c(1,1,2),col=c(1,2,3))


plot(y=colMeans(lasso.error.naive2), x=log(cv_scale_values), ylab ="L1 error", main = TeX("$p=600,n=300,q=100,s=5,\\sigma_{\\nu}=1"), type="l", ylim=c(-0.5,6),xlab=TeX("log$(\\phi)$"))
lines(y=colMeans(lasso.error2),x= log(cv_scale_values), col=2, type="l")
abline(h = b.norm.mean2,col=3,lty=2)
abline(v = 0,col=1,lty=2)
legend(x="bottomright", legend=c("Standard Lasso", "Trimmed Lasso",TeX("$|| b||_1$")),lty=c(1,1,2),col=c(1,2,3))
dev.off()
