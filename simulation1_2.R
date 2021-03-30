library(glmnet)
library("corpcor")
library(Rfast)
library("latex2exp")
library(grDevices)
library(tilting)

### Simulation corresponding to Figure 3



# Model parameters
p = 600
q = 6 
s = 5 
sigma_E = 2
sigma_nu = 1

beta = matrix(c(rep(1,s), rep(0, p-s)), nrow=p) 

n_values = seq(30,p, by=50)
N =1000 # Number of independent simulations

set.seed(2023)

# Model parameters
q = 6

lasso.error.naive3 = matrix(data=NA, ncol = length(n_values), nrow = N)
lasso.error3 = matrix(data=NA, ncol = length(n_values), nrow = N)

for (i in 1:length(n_values)) {
  print("i=")
  print(i)
  n = n_values[i]
  for (j in 1:N) {
    if(j%%100 == 0){
      print("j=")
      print(j)
    }
    #print(j)
    
    # Generating gamma
    Gamma = matrix(rnorm(q*p, mean=0, sd=6), nrow=q,ncol=p)
    
    # Generating data
    nu =  rnorm(n, mean=0, sd=sigma_nu)
    H = matrix(rnorm(n*q, mean=0, sd=1),nrow=n, ncol=q)
    E = matrix(rnorm(p*n,mean=0,sd = sigma_E), nrow=n,ncol=p)
    X = H%*% Gamma + E 
    Y = X%*% beta +nu 
    
    # Naive lasso
    cv.fit.naive = cv.glmnet(y=Y,x=X, intercept = FALSE, standardize=FALSE,penalty.factor = col.norm(X))
    beta.hat.naive = as.vector(coef(cv.fit.naive,s= cv.fit.naive$lambda.min))[-1]
    lasso.error.naive3[j,i] = sum(abs(beta.hat.naive - beta))
    
    # Spectral deconfounder Lasso as in Cevid (2018)
    
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
    F_matrix = tcrossprod(F_matrix,U)
    
    X_tilde = F_matrix %*% X
    Y_tilde = F_matrix %*% Y
    
    # Spectral Lasso
    cv.fit = cv.glmnet(y=Y_tilde,x=X_tilde,intercept = FALSE, standardize=FALSE,penalty.factor = col.norm(X_tilde))
    beta.hat = as.vector(coef(cv.fit,s= cv.fit$lambda.min))[-1]
    
    lasso.error3[j,i] = sum(abs(beta.hat - beta))
  }
}

# Computing the L1 norm of the bias b
b.norm.mean3 =0


# Plotting
directory = "INSERT HERE"
setwd(directory)
filename = "simulation1_2.pdf"
pdf(file = filename, width = 10*3/4, height = 8*3/4)

plot(y=colMeans(lasso.error.naive3), x=n_values, ylab ="L1 error", main = TeX("$p=600,q=6$. Only factor model. Large confounder variance."), type="l", ylim=c(-0.5,5),xlab="n")
lines(y=colMeans(lasso.error3),x= n_values, col=2, type="l")
abline(h = b.norm.mean3,col=3,lty=2)
legend(x="topright", legend=c("Standard Lasso", "Trimmed Lasso",TeX("$|| b||_1$")),lty=c(1,1,2),col=c(1,2,3))
dev.off()
