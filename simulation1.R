library(glmnet)
library("corpcor")
library(Rfast)
library("latex2exp")

### Simulation corresponding to Figure 2

set.seed(2023)

# Case 1: small q, confounded case

# Model parameters
p = 600
q = 6 
s = 5 
sigma_E = 2
sigma_nu = 1

beta = matrix(c(rep(1,s), rep(0, p-s)), nrow=p) 

n_values = seq(30,p, by=50)
N = 2000 # Number of independent simulations

# Matrix of L1 errors of "naive" lasso regression:
lasso.error.naive1 = matrix(data=NA, ncol = length(n_values), nrow = N) 
# Matrix of L1 errors of the spectral deconfounder lasso:
lasso.error1 = matrix(data=NA, ncol = length(n_values), nrow = N)

for (i in 1:length(n_values)) {
  #print("i=")
  #print(i)
  n = n_values[i]
  for (j in 1:N) {
    if(j%%100 == 0){
      print("j=")
      print(j)
    }

    # Generating delta and gamma
    delta = rnorm(q, mean=0, sd=1)# len q
    Gamma = matrix(rnorm(q*p, mean=0, sd=1), nrow=q,ncol=p)
    
    # Generating data
    nu =  rnorm(n, mean=0, sd=sigma_nu) #error term in (2.1)
    H = matrix(rnorm(n*q, mean=0, sd=1),nrow=n, ncol=q)# iid rows, gaussian, mean 0, cov matrix I
    E = matrix(rnorm(p*n,mean=0, sd = sigma_E), nrow=n,ncol=p)#error term in (2.2)
    X = H%*% Gamma + E #(2.2)
    Y = X%*% beta + H %*% delta + nu #(2.1)
    
    # Naive lasso
    cv.fit.naive = cv.glmnet(y=Y,x=X, intercept = FALSE, standardize=FALSE)
    beta.hat.naive = as.vector(coef(cv.fit.naive,s= cv.fit.naive$lambda.min))[-1]
    lasso.error.naive1[j,i] = sum(abs(beta.hat.naive - beta))
    
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
    F_matrix = Tcrossprod(F_matrix,U)
    
    X_tilde = mat.mult(F_matrix,X)
    Y_tilde = F_matrix %*% Y
    
    # Spectral Lasso
    cv.fit = cv.glmnet(y=Y_tilde,x=X_tilde,intercept = FALSE, standardize=FALSE)
    beta.hat = as.vector(coef(cv.fit,s= cv.fit$lambda.min))[-1]
    
    lasso.error1[j,i] = sum(abs(beta.hat - beta))
  }
}

# Computing the average L1 norm of the bias b term b
set.seed(2023)
b.norm1 = c()
for (i in 1:N) {
  delta = rnorm(q, mean=0, sd=1)# len q
  Gamma = matrix(rnorm(q*p, mean=0, sd=1), nrow=q,ncol=p)
  Sigma = Crossprod(Gamma,Gamma)
  for (k in 1:p) {
    Sigma[k,k] = Sigma[k,k]+1
  }
  SigmaInv= spdinv(Sigma)
  b = Tcrossprod(SigmaInv, Gamma) %*%delta
  b.norm1 = c(b.norm1, norm(b, type="1"))
}

b.norm.mean1 =mean(b.norm1)





# Case 2: Large q, confounded case
set.seed(2023)

# Model parameters
q = 100 

lasso.error.naive2 = matrix(data=NA, ncol = length(n_values), nrow = N)
lasso.error2 = matrix(data=NA, ncol = length(n_values), nrow = N)

for (i in 1:length(n_values)) {
  #print("i=")
  #print(i)
  n = n_values[i]
  for (j in 1:N) {
    if(j%%100 == 0){
      print("j=")
      print(j)
    }
    # Generating delta and gamma
    delta = rnorm(q, mean=0, sd=1)# len q
    Gamma = matrix(rnorm(q*p, mean=0, sd=1), nrow=q,ncol=p)
    
    # Generating data
    nu =  rnorm(n, mean=0, sd=sigma_nu) #error term in (2.1)
    H = matrix(rnorm(n*q, mean=0, sd=1),nrow=n, ncol=q)# iid rows, gaussian, mean 0, cov matrix I
    E = matrix(rnorm(p*n,mean=0, sd = sigma_E), nrow=n,ncol=p)#error term in (2.2)
    X = H%*% Gamma + E #(2.2)
    Y = X%*% beta + H %*% delta + nu #(2.1)
    
    # Naive lasso
    cv.fit.naive = cv.glmnet(y=Y,x=X, intercept = FALSE, standardize=FALSE)
    beta.hat.naive = as.vector(coef(cv.fit.naive,s= cv.fit.naive$lambda.min))[-1]
    lasso.error.naive2[j,i] = sum(abs(beta.hat.naive - beta))
    
    # Spectral deconfoundede Lasso as in Cevid (2018)
    
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
    
    # Spectral Lasso
    cv.fit = cv.glmnet(y=Y_tilde,x=X_tilde,intercept = FALSE, standardize=FALSE)
    beta.hat = as.vector(coef(cv.fit,s= cv.fit$lambda.min))[-1]
    
    lasso.error2[j,i] = sum(abs(beta.hat - beta))
  }
}

# Computing the average L1 norm of the bias b
set.seed(2023)
b.norm2 = c()
for (i in 1:N) {
  delta = rnorm(q, mean=0, sd=1)# len q
  Gamma = matrix(rnorm(q*p, mean=0, sd=1), nrow=q,ncol=p)
  Sigma = Crossprod(Gamma,Gamma)
  for (k in 1:p) {
    Sigma[k,k] = Sigma[k,k]+1
  }
  SigmaInv= spdinv(Sigma)
  b = Tcrossprod(SigmaInv, Gamma) %*%delta
  b.norm2 = c(b.norm2, norm(b, type="1"))
}

b.norm.mean2 =mean(b.norm2)




# Case 3: small q, no confounding, but X is generated from factor model
set.seed(2023)

# Model parameters
q = 6

lasso.error.naive3 = matrix(data=NA, ncol = length(n_values), nrow = N)
lasso.error3 = matrix(data=NA, ncol = length(n_values), nrow = N)

for (i in 1:length(n_values)) {
  #print("i=")
  #print(i)
  n = n_values[i]
  for (j in 1:N) {
    if(j%%100 == 0){
      print("j=")
      print(j)
    }
    #print(j)
    
    # Generating gamma
    Gamma = matrix(rnorm(q*p, mean=0, sd=1), nrow=q,ncol=p)
    
    # Generating data
    nu =  rnorm(n, mean=0, sd=sigma_nu) #error term in (2.1)
    H = matrix(rnorm(n*q, mean=0, sd=1),nrow=n, ncol=q)# iid rows, gaussian, mean 0, cov matrix I
    E = matrix(rnorm(p*n,mean=0,sd = sigma_E), nrow=n,ncol=p)#error term in (2.2)
    X = H%*% Gamma + E #(2.2)
    Y = X%*% beta +nu #(2.1)
    
    # Naive lasso
    cv.fit.naive = cv.glmnet(y=Y,x=X, intercept = FALSE, standardize=FALSE)
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
    F_matrix = Tcrossprod(F_matrix,U)
    
    X_tilde = mat.mult(F_matrix,X)
    Y_tilde = F_matrix %*% Y
    
    # Spectral Lasso
    cv.fit = cv.glmnet(y=Y_tilde,x=X_tilde,intercept = FALSE, standardize=FALSE)
    beta.hat = as.vector(coef(cv.fit,s= cv.fit$lambda.min))[-1]
    
    lasso.error3[j,i] = sum(abs(beta.hat - beta))
  }
}

# Computing the L1 norm of the bias b
b.norm.mean3 =0




# Case 4: no confounding, no factor model
set.seed(2023)


lasso.error.naive4 = matrix(data=NA, ncol = length(n_values), nrow = N)
lasso.error4 = matrix(data=NA, ncol = length(n_values), nrow = N)

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
    # Generating delta and gamma
    #delta = rnorm(q, mean=0, sd=1)# len q
    #Gamma = matrix(rnorm(q*p, mean=0, sd=1), nrow=q,ncol=p)
    
    # Generating data
    nu =  rnorm(n, mean=0, sd=sigma_nu) #error term in (2.1)
    #H = matrix(rnorm(n*q, mean=0, sd=1),nrow=n, ncol=q)# iid rows, gaussian, mean 0, cov matrix I
    E = matrix(rnorm(p*n,mean=0, sd = sigma_E), nrow=n,ncol=p)#error term in (2.2)
    X = E #(2.2)
    Y = X%*% beta +nu #(2.1)
    
    # Naive lasso
    cv.fit.naive = cv.glmnet(y=Y,x=X, intercept = FALSE, standardize=FALSE)
    beta.hat.naive = as.vector(coef(cv.fit.naive,s= cv.fit.naive$lambda.min))[-1]
    lasso.error.naive4[j,i] = sum(abs(beta.hat.naive - beta))
    
    # Spectral Lasso as in Cevid (2018)
    
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
    
    # Spectral Lasso
    cv.fit = cv.glmnet(y=Y_tilde,x=X_tilde,intercept = FALSE, standardize=FALSE)
    beta.hat = as.vector(coef(cv.fit,s= cv.fit$lambda.min))[-1]
    
    lasso.error4[j,i] = sum(abs(beta.hat - beta))
  }
}

# Computing the L1 norm of the bias b


b.norm.mean4 =0
dev.new(width=10,height=8,noRStudioGD = TRUE)
par(mfrow=c(2,2),mar=c(4,4,4,4))

plot(y=colMeans(lasso.error.naive1), x=n_values, ylab ="L1 error", main = TeX("$p=600,q=6,s=5,\\sigma_{\\\\nu\\}=1$. Confounded model"), type="l", ylim=c(-0.5,5),xlab="n")
lines(y=colMeans(lasso.error1),x= n_values, col=2, type="l")
abline(h = b.norm.mean1,col=3,lty=2)
legend(x="topright", legend=c("Naive Lasso", "Spectral Lasso",TeX("$|| b||_1$")),lty=c(1,1,2),col=c(1,2,3))

plot(y=colMeans(lasso.error.naive2), x=n_values, ylab ="L1 error", main = TeX("$p=600,q=100,s=5,\\sigma_{\\\\nu\\}=1$. Confounded model"), type="l", ylim=c(-0.5,15),xlab="n")
lines(y=colMeans(lasso.error2),x= n_values, col=2, type="l")
abline(h = b.norm.mean2,col=3,lty=2)
legend(x="topright", legend=c("Naive Lasso", "Spectral Lasso",TeX("$|| b||_1$")),lty=c(1,1,2),col=c(1,2,3))


plot(y=colMeans(lasso.error.naive3), x=n_values, ylab ="L1 error", main = TeX("$p=600,q=6,s=5,\\sigma_{\\\\nu\\}=1$. Only factor model."), type="l", ylim=c(-0.5,5),xlab="n")
lines(y=colMeans(lasso.error3),x= n_values, col=2, type="l")
abline(h = b.norm.mean3,col=3,lty=2)
legend(x="topright", legend=c("Naive Lasso", "Spectral Lasso",TeX("$|| b||_1$")),lty=c(1,1,2),col=c(1,2,3))


plot(y=colMeans(lasso.error.naive4), x=n_values,ylab ="L1 error", main = TeX("$p=600,q=6,s=5,\\sigma_{\\\\nu\\}=1$, No confounding, no factor model."), type="l", ylim=c(-0.5,5),xlab="n")
lines(y=colMeans(lasso.error4),x= n_values, col=2, type="l")
abline(h = b.norm.mean4,col=3,lty=2)
legend(x="topright",legend=c("Naive Lasso", "Spectral Lasso",TeX("$|| b||_1$")),lty=c(1,1,2),col=c(1,2,3))
#cex=0.8,y.intersp=0.5,x.intersp=2

quartz.save(file="/Users/peraugust/Documents/cambridge/essay/simulations/simulation1.jpeg",type="jpeg",dpi=300)
