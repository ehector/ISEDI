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
    
    R <- solve(G1)%*%G2
    sigma2_hat <- sum((y_2 - X_2%*%hat_beta_2)^2)/n_2
    sd <- z_q * sqrt(sigma2_hat)
    
    G2_eig <- eigen(G2)
    G2_eigenval <- sort(G2_eig$values, decreasing=FALSE)
    K2 <- diag(G2_eig$values)
    P2 <- t(G2_eig$vectors)
    G2_sq <- t(P2) %*% sqrt(K2) %*% P2
    M_inv_eig <- eigen(G2_sq %*% solve(G1) %*% G2_sq)
    K12 <- M_inv_eig$values
    delta <- hat_beta_2 - hat_beta_1
    
    lambda_opt <- (sigma2_hat/n_2)*min(K12 / G2_eigenval)/max(delta^2)
    
    estimates <- do.call(rbind, lapply(1:L, function(l){
      lambda <- lambda_seq[l]
      estimate <- solve(G2 + lambda*G1) %*% (
        Xy2/n_2 + lambda*Gbeta_1
      )
      S <- solve(G2 + lambda*G1)
      IS <- solve(diag(p) - lambda*S%*%G1)
      D <- IS %*% S %*%G2 %*% S %*% t(IS)
      D_half <- sqrt(diag(D))
      center <- IS%*% (estimate - lambda*S%*%Gbeta_1)
      width <- sd * D_half/n_2_half
      lower <- center - width
      upper <- center + width
      return(data.frame(parameter=pars, estimate=estimate, lower_CI=lower, upper_CI=upper, lambda=lambda))
    }))
    
    estimate <- solve(G2 + lambda_opt*G1) %*% (
      Xy2/n_2 + lambda_opt*Gbeta_1
    )
    S <- solve(G2 + lambda_opt*G1)
    IS <- solve(diag(p) - lambda_opt*S%*%G1)
    D <- IS %*% S %*%G2 %*% S %*% t(IS)
    D_half <- sqrt(diag(D))
    center <- IS%*% (estimate - lambda_opt*S%*%Gbeta_1)
    width <- sd * D_half/n_2_half
    lower <- center - width
    upper <- center + width
    
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
    A2 <- diag(x=drop(exp(X_2%*%hat_beta_2)/(1+exp(X_2%*%hat_beta_2))^2))
    Delta <- drop(exp(X_1%*%hat_beta_2)/(1+exp(X_1%*%hat_beta_2)) - exp(X_1%*%hat_beta_1)/(1+exp(X_1%*%hat_beta_1)))
    
    v1 <- t(X_1)%*%A1%*%X_1/n_1
    v2 <- t(X_2)%*%A2%*%X_2/n_2
    vbeta_1 <- v1%*%hat_beta_1
    v1_inv <- solve(v1)
    v2_inv <- solve(v2)
    # R <- solve(v1)%*%v2
    Delta_tilde <- v1_inv %*% t(X_1) %*%Delta %*%t(Delta)%*% X_1 %*% v1_inv /n_1^2
    
    v2_eig <- eigen(v2)
    v2_eigenval <- sort(v2_eig$values, decreasing=FALSE)
    K2 <- diag(v2_eig$values)
    P2 <- t(v2_eig$vectors)
    v2_sq <- t(P2) %*% sqrt(K2) %*% P2
    M_inv_eig <- eigen(v2_sq %*% solve(v1) %*% v2_sq)
    K12 <- M_inv_eig$values
    
    lambda_opt <- (1/n_2)*min(K12 / v2_eigenval)/max(diag(Delta_tilde))
    
    estimates <- do.call(rbind, lapply(1:L, function(l){
      lambda <- lambda_seq[l]
      estimate <- 
        optim(par=hat_beta_2, fn=binom_min_func, gr=binom_min_func_deriv,
              hat_beta_1=hat_beta_1, X_1=X_1, X_2=X_2, Z_2=Z_2, y_1=y_1, y_2=y_2, 
              hat_eta_1=hat_eta_1, hat_mu_1=hat_mu_1, n_1=n_1, n_2=n_2, lambda=lambda, method="BFGS")$par
      S <- v2 + lambda*v1
      S_inv <- solve(S)
      J_inv <- S_inv%*%v2%*%S_inv
      IS <- solve(diag(p) - lambda*S_inv%*%v1)
      D <- IS %*% J_inv %*% t(IS)
      D_half <- sqrt(diag(D))
      center <- IS%*% (estimate - lambda*S_inv%*%vbeta_1)
      width <- z_q * D_half/n_2_half
      lower <- center - width
      upper <- center + width
      return(data.frame(parameter=pars, estimate=estimate, lower_CI=lower, upper_CI=upper, lambda=lambda))
    }))
    
    estimate <- 
      optim(par=hat_beta_2, fn=binom_min_func, gr=binom_min_func_deriv,
            hat_beta_1=hat_beta_1, X_1=X_1, X_2=X_2, Z_2=Z_2, y_1=y_1, y_2=y_2, 
            hat_eta_1=hat_eta_1, hat_mu_1=hat_mu_1, n_1=n_1, n_2=n_2, lambda=lambda_opt, method="BFGS")$par
    S <- v2 + lambda_opt*v1
    S_inv <- solve(S)
    J_inv <- solve(S%*%v2_inv%*%S)
    IS <- solve(diag(p) - lambda_opt*S_inv%*%v1)
    D <- IS %*% J_inv %*% t(IS)
    D_half <- sqrt(diag(D))
    center <- IS%*% (estimate - lambda_opt*S_inv%*%vbeta_1)
    width <- z_q * D_half/n_2_half
    lower <- center - width
    upper <- center + width
    
    estimates_opt <- data.frame(parameter=pars, estimate=estimate, lower_CI=lower, upper_CI=upper, lambda=lambda_opt)
  }
  return(list(MLEs=cbind(hat_beta_1=hat_beta_1, hat_beta_2=hat_beta_2), estimates=estimates, lambda=lambda_seq, 
              estimates_opt=estimates_opt, lambda_opt=lambda_opt))
}
