rm(list = ls())
library(glmnet)
library(tilting)
library("corpcor")
library(Rfast)
library("latex2exp")
library(sail)
library(caret)


# Simulation corresponding to Figure 7
ps = seq(250,1250,by=100)
N = 2000
set.seed(2024)
n = 200
q =30 
s = 5 
sigma_E = 1
sigma_nu = 1



lambda_inflate = 1.2



sigma_hats_naive = matrix(data=NA, ncol = length(ps),nrow=N)
sigma_hats_dd = matrix(data=NA, ncol = length(ps),nrow=N)
sigma_hats_dd_fold = matrix(data=NA, ncol = length(ps),nrow=N)

v_1_hats_naive = matrix(data=NA, ncol = length(ps),nrow=N)
v_1_hats_dd = matrix(data=NA, ncol = length(ps),nrow=N)
v_1_hats_dd_fold = matrix(data=NA, ncol = length(ps),nrow=N)

beta_1_hats_naive = matrix(data=NA, ncol = length(ps),nrow=N)
beta_1_hats_dd = matrix(data=NA, ncol = length(ps),nrow=N)

v_2_hats_naive = matrix(data=NA, ncol = length(ps),nrow=N)
v_2_hats_dd = matrix(data=NA, ncol = length(ps),nrow=N)
v_2_hats_dd_fold = matrix(data=NA, ncol = length(ps),nrow=N)

beta_2_hats_naive = matrix(data=NA, ncol = length(ps),nrow=N)
beta_2_hats_dd = matrix(data=NA, ncol = length(ps),nrow=N)

coverage_1_naive = c()
coverage_1_spectral = c()
coverage_1_spectral_fold = c()

coverage_2_naive = c()
coverage_2_spectral = c()
coverage_2_spectral_fold = c()

j = s+1

for (l in 1:length(ps)) {
    
  p = ps[l]
  beta = matrix(c(rep(1,s), rep(0, p-s)), nrow=p)
  beta_naive_error = c()
  beta_spectral_error = c()
  
  beta_1_naive_DD_vec = c()
  beta_2_naive_DD_vec = c()
  beta_1_DD_vec = c()
  beta_2_DD_vec = c()
  
  est_sd_naive_1_vec = c()
  est_sd_naive_2_vec = c()
  est_sd_1_vec = c()
  est_sd_2_vec = c()
  
  est_sd_1_boot_vec = c()
  est_sd_2_boot_vec = c()
  
  est_sd_1_fold_vec = c()
  est_sd_2_fold_vec = c()
  
  sigma_e_est = c()
  sigma_e_bootstrap_vec = c()
  sigma_e_est_fold = c()
  sigma_naive_hat_vec = c()
  
  
  
  
  
  
  for (r in 1:N) {
    
    
    delta = rnorm(q, mean=0, sd=1)# len q
    
    Gamma = matrix(rnorm(q*p, mean=0, sd=1), nrow=q,ncol=p)

    nu =  rnorm(n, mean=0, sd=sigma_nu) #error term in (2.1)
    print(r)
    H = matrix(rnorm(n*q, mean=0, sd=1),nrow=n, ncol=q)# iid rows, gaussian, mean 0, cov matrix I
    E = matrix(rnorm(p*n,mean=0, sd = sigma_E), nrow=n)#error term in (2.2)
    X = H%*% Gamma + E #(2.2)
    Y = X%*% beta + H %*% delta + nu #(2.1)
    
    #Standard/naive regression: here with penalty weights
    

    
    lambda = cv.glmnet(y=Y,x=X, penalty.factor = col.norm(X),parallel = F,intercept=F)
    beta.hat.naive = as.vector(coef(lambda,s= lambda_inflate*lambda$lambda.min))[-1]

    
    #doing TRIM: 
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
    
    

    lambda = cv.glmnet(y=Y_tilde,x=X_tilde,penalty.factor = col.norm(X_tilde), parallel = F,intercept=F,keep=T)
    beta.hat = as.vector(coef(lambda,s= lambda_inflate*lambda$lambda.min))[-1]
    
    beta_spectral_error= c(beta_spectral_error, sum(abs(beta.hat - beta)) )# deconfounded
    beta_naive_error = c(beta_naive_error,sum(abs(beta.hat.naive - beta))) #regular lasso
    
    trQ2 =sum(diag(F_matrix%*%F_matrix))
    sigma_e_hat =1/sqrt(trQ2)  *norm(Y_tilde -X_tilde%*%beta.hat,type="2")
    
    sigma_e_folds = norm(Y_tilde -lambda$fit.preval[,which(lambda$lambda== lambda$lambda.min)],type="2") / sqrt(trQ2)
    
    sigma_e_est = c(sigma_e_est, sigma_e_hat)
    sigma_e_est_fold = c(sigma_e_est_fold, sigma_e_folds)
    
    # Debiasing
    
    
    # Standard/naive debiasing first. 
    # We consider a significant variable (1) and a non-significant one (s+1)
    
    # Inference on beta_1
    naiveinternallasso1CV = cv.glmnet(y=X[,1],x=X[,-1],penalty.factor = col.norm(X[,-1]),parallel = F,intercept=F)
    naiveinternallasso1_coeff = as.vector(coef(naiveinternallasso1CV,s= naiveinternallasso1CV$lambda.min))
    naive_residual1 = (X[,1]) - X[,-1]%*%naiveinternallasso1_coeff[-1]
    PjZjTPjX = crossprod(naive_residual1,X[,1])
    beta_1_naive_DD = crossprod(naive_residual1,Y - X[,-1]%*%beta.hat.naive[-1]) /(PjZjTPjX)
    beta_1_naive_DD_vec = c(beta_1_naive_DD_vec, beta_1_naive_DD)
    sigma_naive_hat_1 = sqrt(1/n) *norm(Y-X%*%beta.hat.naive,type="2")
    sigma_naive_hat_vec = c(sigma_naive_hat_vec,sigma_naive_hat_1)
    est_sd_naive_1 = sigma_naive_hat_1 *norm(naive_residual1, type="2") / (PjZjTPjX)
    est_sd_naive_1
    est_sd_naive_1_vec = c(est_sd_naive_1_vec, est_sd_naive_1)
    
    #Inference on beta_(s+1)
    naiveinternallasso2CV = cv.glmnet(y=X[,j],x=X[,-j],penalty.factor = col.norm(X[,-j]),parallel = F,intercept=F)
    naiveinternallasso2_coeff = as.vector(coef(naiveinternallasso2CV,s= naiveinternallasso2CV$lambda.min))
    naive_residual2 = (X[,j]) - X[,-j]%*%naiveinternallasso2_coeff[-1]
    PjZjTPjX = crossprod(naive_residual2,X[,j])
    beta_2_naive_DD = crossprod(naive_residual2,Y- X[,-j]%*%beta.hat.naive[-j]) /(PjZjTPjX)
    beta_2_naive_DD_vec = c(beta_2_naive_DD_vec, beta_2_naive_DD)
    sigma_naive_hat_2 = sqrt(1/n) *norm(Y-X%*%beta.hat.naive,type="2")
    est_sd_naive_2 = sigma_naive_hat_2 *norm(naive_residual2, type="2") / (PjZjTPjX)
    est_sd_naive_2
    est_sd_naive_2_vec = c(est_sd_naive_2_vec, est_sd_naive_2)
    
    
    
    
    # Now we compute Doubly Debiased estimator 
    
    # Inference on beta_1 with DD
    svd_1 = fast.svd(X[,-1])
    d1 = svd_1$d
    tau1 = median(d1)
    U1 = svd_1$u
    d_tilde_1 = pmin(d1,tau1)
    
    F1 = U1[,]
    for (i in 1:dim(U1)[1]) {
      F1[,i] = d_tilde_1[i]/d1[i] * F1[,i]
    }
    F1=  tcrossprod(F1,U1)
    
    X_tilde_1 = F1 %*% X[,-1]
    Y_tilde_1 = F1 %*% X[,1]
    
    internallasso1CV = cv.glmnet(y=Y_tilde_1,x=X_tilde_1,penalty.factor = col.norm(X_tilde_1),parallel = F,intercept=F)
    internallasso1_coeff = as.vector(coef(internallasso1CV,s= internallasso1CV$lambda.min))
    residual1 = (X[,1]) - X[,-1]%*%internallasso1_coeff[-1]
    PjZ = F1%*%residual1
    beta_1_DD = crossprod(PjZ,F1%*%(Y- X[,-1]%*%beta.hat[-1])) /crossprod(PjZ,F1%*%X[,1])
    beta_1_DD
    beta_1_DD_vec = c(beta_1_DD_vec, beta_1_DD)

    P2 = F1%*%F1
    P4 = P2%*%P2
    est_sd_1 = sigma_e_hat *sqrt(crossprod(residual1, P4%*%residual1)) / abs(crossprod(residual1,P2%*%X[,1]))
    est_sd_1_vec = c(est_sd_1_vec, est_sd_1)
    est_sd_1_fold = sigma_e_folds /sigma_e_hat * est_sd_1
    est_sd_1_fold_vec = c(est_sd_1_fold_vec, est_sd_1_fold)
    
    
    #Inference on beta_(s+1) with DD
    #j=s+1
    svd_2 = fast.svd(X[,-j])
    d2 = svd_2$d
    tau2 = median(d2)
    U2 = svd_2$u
    d_tilde_2 = pmin(d2,tau2)
    
    F2 = U2[,]
    for (i in 1:dim(U2)[1]) {
      U2[,i] = d_tilde_2[i]/d2[i] * U2[,i]
    }
    F2= F2 %*% t(U2)
    
    X_tilde_2 = F2 %*% X[,-j]
    Y_tilde_2 = F2 %*% X[,j]
    
    internallasso2CV = cv.glmnet(y=Y_tilde_2,x=X_tilde_2,penalty.factor = col.norm(X_tilde_2),parallel = F,intercept=F)
    internallasso2_coeff = as.vector(coef(internallasso2CV,s= internallasso2CV$lambda.min))
    residual2 = (X[,j]) - X[,-j]%*%internallasso2_coeff[-1]
    PjZ = F2%*%residual2
    beta_2_DD = crossprod(PjZ,F2%*%(Y- X[,-j]%*%beta.hat[-j])) /crossprod(PjZ,F2%*%X[,j])
    beta_2_DD
    beta_2_DD_vec = c(beta_2_DD_vec, beta_2_DD)

    P2 = F2%*%F2
    P4 = P2%*%P2
    est_sd_2 = sigma_e_hat *sqrt(crossprod(residual2, P4%*%residual2)) / abs(crossprod(residual2,P2%*%X[,j]))
    est_sd_2
    est_sd_2_vec = c(est_sd_2_vec, est_sd_2)

    est_sd_2_fold = sigma_e_folds /sigma_e_hat * est_sd_2
    est_sd_2_fold_vec = c(est_sd_2_fold_vec, est_sd_2_fold)
    
    
  }
  sigma_hats_naive[,l] = sigma_naive_hat_vec
  sigma_hats_dd[,l] = sigma_e_est
  sigma_hats_dd_fold[,l] = sigma_e_est_fold
  
  v_1_hats_naive[,l] = est_sd_naive_1_vec
  v_1_hats_dd[,l] = est_sd_1_vec
  v_1_hats_dd_fold[,l] = est_sd_1_fold_vec
  
  beta_1_hats_naive = beta_1_naive_DD_vec
  beta_1_hats_dd = beta_1_DD_vec
  
  v_2_hats_naive[,l] = est_sd_naive_2_vec
  v_2_hats_dd[,l] = est_sd_2_vec
  v_2_hats_dd_fold[,l] = est_sd_2_fold_vec
  
  beta_2_hats_naive[,l] = beta_2_naive_DD_vec
  beta_2_hats_dd[,l] = beta_2_DD_vec
  
  
  coverage_1_naive = c(coverage_1_naive,
                       mean((beta[1] <=beta_1_naive_DD_vec + qnorm(0.975)*(est_sd_naive_1_vec))&(beta[1] >= beta_1_naive_DD_vec  - qnorm(0.975)*(est_sd_naive_1_vec)) ))
  coverage_1_spectral = c(coverage_1_spectral, 
                      mean((beta[1] <=beta_1_DD_vec + qnorm(0.975)*est_sd_1_vec)&(beta[1] >=beta_1_DD_vec - qnorm(0.975)*est_sd_1_vec) ))
  coverage_1_spectral_fold = c(coverage_1_spectral_fold, 
                      mean((beta[1] <=beta_1_DD_vec + qnorm(0.975)*est_sd_1_fold_vec)&(beta[1] >=beta_1_DD_vec - qnorm(0.975)*est_sd_1_fold_vec) )        )
  
  coverage_2_naive = c(coverage_2_naive,
                       mean((beta[j] <=beta_2_naive_DD_vec + qnorm(0.975)*(est_sd_naive_2_vec))&(beta[j] >= beta_2_naive_DD_vec  - qnorm(0.975)*(est_sd_naive_2_vec)) ))
  coverage_2_spectral = c(coverage_2_spectral, 
                          mean((beta[j] <=beta_2_DD_vec + qnorm(0.975)*est_sd_2_vec)&(beta[j] >=beta_2_DD_vec - qnorm(0.975)*est_sd_2_vec) ))
  coverage_2_spectral_fold = c(coverage_2_spectral_fold, 
                              mean((beta[j] <=beta_2_DD_vec + qnorm(0.975)*est_sd_2_fold_vec)&(beta[j] >=beta_2_DD_vec - qnorm(0.975)*est_sd_2_fold_vec) )        )
  
}





directory = "INSERT HERE"
filename = "simulation4_2.pdf"
setwd(directory)
pdf(file = filename, width = 10, height = 5)
par(mfrow=c(1,2),mar=c(4,4,4,4))

plot(y=coverage_1_naive, x=ps, ylab ="Coverage frequency", main = TeX("Coverage of $\\beta_1$. $\ n=200,q=30,\\sigma_{E}=1$"), type="l", ylim=c(0.1,1.0),xlab=TeX("$p$"))
lines(y=coverage_1_spectral,x= ps, col=2, type="l")
lines(y=coverage_1_spectral_fold,x= ps, col=3, type="l")

#Plotting grid:
for (i in 1:10) {
  abline(h=i/10,col = "lightgray", lty = "dotted",
         lwd = par("lwd"))
}
abline(h = 0.95,col=1,lty=2)
legend(x="bottomleft", legend=c("Debiased Lasso", "DD Lasso","DD Lasso Fold"),lty=c(1,1,1),col=c(1,2,3))


plot(y=coverage_2_naive, x=ps, ylab ="Coverage frequency", main = TeX("Coverage of $\\beta_{s+1}$. $\ n=200,q=30,\\sigma_{E}=1$"), type="l", ylim=c(0.1,1.0),xlab=TeX("$p$"))
lines(y=coverage_2_spectral,x= ps, col=2, type="l")
lines(y=coverage_2_spectral_fold,x= ps, col=3, type="l")

#Plotting grid:
for (i in 1:10) {
  abline(h=i/10,col = "lightgray", lty = "dotted",
         lwd = par("lwd"))
}
abline(h = 0.95,col=1,lty=2)
legend(x="bottomleft", legend=c("Debiased Lasso", "DD Lasso","DD Lasso Fold"),lty=c(1,1,1),col=c(1,2,3))
dev.off()


# Average coverages:
mean(coverage_1_naive) #0.3321818
mean(coverage_1_spectral) #0.8844545
mean(coverage_1_spectral_fold) #0.9301364

mean(coverage_2_naive) #0.9241364
mean(coverage_2_spectral) #0.9341364
mean(coverage_2_spectral_fold) #0.9558182


