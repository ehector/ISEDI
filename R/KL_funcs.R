LDW_meta <- function(X_1, X_2, y_1, y_2, lambda_seq, family){
  n_1 <- length(y_1)
  n_2 <- length(y_2)
  n_2_half <- sqrt(n_2)
  p <- ncol(X_2)
  z_q <- qnorm(0.025, lower.tail=FALSE)
  pars <- paste0("beta",1:p)
  L <- length(lambda_seq)
  
  if(family=="gaussian"){
    Xy2 <- t(X_2)%*%y_2
    G1 <- t(X_1)%*%X_1/n_1
    G2 <- t(X_2)%*%X_2/n_2
    
    hat_beta_1 <- as.vector(solve(G1)%*%t(X_1)%*%y_1)/n_1
    hat_beta_2 <- as.vector(solve(G2)%*%Xy2)/n_2
    Gbeta_1 <- G1%*%hat_beta_1
    
    sigma2_hat <- sum((y_2 - X_2%*%hat_beta_2)^2)/(n_2-p)
    sd <- z_q * sqrt(sigma2_hat)
    
    delta_p <- hat_beta_2 - hat_beta_1
    A <- delta_p%*%t(delta_p) - sigma2_hat*solve(G2)/n_2
    A_eig <- eigen(A)
    if(sum(A_eig$values>0)>=1){
      A_eig$values[A_eig$values < 0] <- 0
      delta2 <- A_eig$vectors %*% diag(A_eig$values) %*% t(A_eig$vectors) 
    } else{
      delta2 <- delta_p%*%t(delta_p)
    }

    lambda_opt <- optim(par=0, fn=gaussian_MSE, G1=G1, G2=G2, delta2=delta2, sigma_2=sigma2_hat, n_2=n_2,
                        method="Brent", lower=0, upper=1e10)$par
    
    estimates <- do.call(rbind, lapply(1:L, function(l){
      lambda <- lambda_seq[l]
      estimate <- solve(G2 + lambda*G1) %*% (
        Xy2/n_2 + lambda*Gbeta_1
      )
      S <- solve(G2 + lambda*G1)
      width <- sd * sqrt(diag(S %*% G2 %*% S)) /n_2_half
      lower <- estimate - width
      upper <- estimate + width
      return(data.frame(parameter=pars, estimate=estimate, lower_CI=lower, upper_CI=upper, lambda=lambda))
    }))
    
    estimate <- solve(G2 + lambda_opt*G1) %*% (
      Xy2/n_2 + lambda_opt*Gbeta_1
    )
    S <- solve(G2 + lambda_opt*G1)
    width <- sd * sqrt(diag(S %*% G2 %*% S)) /n_2_half
    lower <- estimate - width
    upper <- estimate + width
    
    estimates_opt <- data.frame(parameter=pars, estimate=estimate, lower_CI=lower, upper_CI=upper, lambda=lambda_opt)
  }
  if(family=="binomial"){
    glm_1 <- glm(y_1 ~ 0 + X_1, family=family)
    glm_2 <- glm(y_2 ~ 0 + X_2, family=family)
    
    hat_beta_1 <- c(coef(glm_1))
    hat_beta_2 <- c(coef(glm_2))
    
    hat_eta_1 <- X_1%*%hat_beta_1
    hat_mu_1 <- exp(hat_eta_1)/(1+exp(hat_eta_1))
    
    Z_2 <- colSums(X_2[y_2==1,])
    
    A1 <- diag(x=drop(exp(X_1%*%hat_beta_1)/(1+exp(X_1%*%hat_beta_1))^2))
    A12 <- diag(x=drop(exp(X_1%*%hat_beta_2)/(1+exp(X_1%*%hat_beta_2))^2))
    A2 <- diag(x=drop(exp(X_2%*%hat_beta_2)/(1+exp(X_2%*%hat_beta_2))^2))
    Delta <- drop(exp(X_1%*%hat_beta_2)/(1+exp(X_1%*%hat_beta_2)) - exp(X_1%*%hat_beta_1)/(1+exp(X_1%*%hat_beta_1)))
    
    v1 <- t(X_1)%*%A1%*%X_1/n_1
    v2 <- t(X_2)%*%A2%*%X_2/n_2
    v2_inv <- solve(v2)
    vbeta_1 <- v1%*%hat_beta_1
    
    XXDelta2 <- t(X_1)%*%Delta%*%t(Delta)%*%X_1
    
    lambda_opt <- optim(par=0, fn=binomial_MSE, v1=v1, v2=v2, XXDelta2=XXDelta2, n_1=n_1, n_2=n_2,
                        method="Brent", lower=0, upper=1e10)$par
    
    estimates <- do.call(rbind, lapply(1:L, function(l){
      lambda <- lambda_seq[l]
      estimate <- 
        optim(par=hat_beta_2, fn=binom_min_func, gr=binom_min_func_deriv,
              hat_beta_1=hat_beta_1, X_1=X_1, X_2=X_2, Z_2=Z_2, y_1=y_1, y_2=y_2, 
              hat_eta_1=hat_eta_1, hat_mu_1=hat_mu_1, n_1=n_1, n_2=n_2, lambda=lambda, method="BFGS")$par
      S <- solve(v2 + lambda*v1)
      J_inv <- S %*% v2 %*% S
      width <- z_q * sqrt(diag(J_inv))/n_2_half
      lower <- estimate - width
      upper <- estimate + width
      return(data.frame(parameter=pars, estimate=estimate, lower_CI=lower, upper_CI=upper, lambda=lambda))
    }))
    
    estimate <- 
      optim(par=hat_beta_2, fn=binom_min_func, gr=binom_min_func_deriv,
            hat_beta_1=hat_beta_1, X_1=X_1, X_2=X_2, Z_2=Z_2, y_1=y_1, y_2=y_2, 
            hat_eta_1=hat_eta_1, hat_mu_1=hat_mu_1, n_1=n_1, n_2=n_2, lambda=lambda_opt, method="BFGS")$par
    S <- solve(v2 + lambda_opt*v1)
    J_inv <- S %*% v2 %*% S
    width <- z_q * sqrt(diag(J_inv))/n_2_half
    lower <- estimate - width
    upper <- estimate + width
    
    estimates_opt <- data.frame(parameter=pars, estimate=estimate, lower_CI=lower, upper_CI=upper, lambda=lambda_opt)
  }
  return(list(MLEs=cbind(hat_beta_1=hat_beta_1, hat_beta_2=hat_beta_2), estimates=estimates, lambda=lambda_seq, 
              estimates_opt=estimates_opt, lambda_opt=lambda_opt))
}

gaussian_MSE <- function(par, G1, G2, delta2, sigma_2, n_2){
  return(
    (sigma_2/n_2)*sum(diag(solve(G2+par*G1)%*%solve(G2+par*G1)%*%G2)) + 
    par^2*sum(diag(G1%*%solve(G2+par*G1)%*%solve(G2+par*G1)%*%G1%*%delta2))
  )
}

binomial_MSE <- function(par, v1, v2, XXDelta2, n_1, n_2){
  S <- v2 + par*v1
  S2_inv <- solve(S%*%S)
  return(drop(sum(diag(S2_inv%*%v2)) + (n_2/n_1^2)*par^2*sum(diag(S2_inv%*%XXDelta2)))/n_2)
}
