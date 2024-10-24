fit_glmm = function(y, X, model_init, libsize, n_iter = 1000, tol = 1e-4){
  
  n = length(y)
  offset = log(libsize)
  
  # initialize 
  eps = rep(0,n)
  beta_init = beta = model_init$coeff
  eta = X %*% beta  + eps + offset
  mu = exp(eta)
  w = (1/mu)*(y-mu) + eta - offset
  tau_init = tau = 0
  
  # iteration 
  for(iter in 1:n_iter) {
    eta = X %*% beta  + eps + offset
    eta = as.vector(eta)
    mu = exp(eta)
    mu[is.infinite(mu)] = median(mu)
    w = (1/mu)*(y-mu) + eta - offset
    vw = 1/mu + tau
    w[abs(w)>1000] = median(w)
    vw[vw>1000] = median(vw)
    
    if(any(vw==Inf)|any(abs(w)==Inf)){
      break
      
    }else{
      XVivX_iv = solve(t(X/vw)%*%X)
      beta = XVivX_iv %*% t(X/vw)%*% w
      res = (w - X %*% beta)/vw
      eps = tau*res    
      
      XViv2X = t(X/vw^2)%*%X
      XViv3X = t(X/vw^3)%*%X
      XViv2X_XVivX = XViv2X%*%XVivX_iv
      XVivres = t(X/vw)%*%res
      trP = sum(1/vw) - sum(diag(XViv2X_XVivX))
      trPP = sum(1/vw^2) - 2*sum(XViv3X*t(XVivX_iv)) + sum(XViv2X_XVivX*t(XViv2X_XVivX))
      yPPy = sum(res^2)
      yPPPy = sum(res^2/vw) - t(XVivres)%*%XVivX_iv%*%XVivres
      G = (yPPy - trP)/2
      H = trPP/2 - yPPPy
      
      GHiv = G/H
      if(abs(GHiv) > max(1, 0.5*tau) ){GHiv = sign(GHiv)*0.01} 
      while (tau - GHiv < 0) {
        GHiv = GHiv*0.5
      }
      tau = as.vector(tau - GHiv)
      
      if(iter > 5 & all(abs(beta-beta_init) < tol) & abs(tau-tau_init) < tol){break}
      beta_init = beta
      tau_init = tau
    }
  }
  
  if(any(vw==Inf)|any(abs(w)==Inf)){
    tau = 0
    beta = model_init$coeff
    converge = model_init$converged
    iter =  model_init$iter
    
    eta = X %*% beta + rep(0,nrow(y)) + offset
    eta = as.vector(eta)
    mu = exp(eta)
    w = (1/mu)*(y-mu) + eta - offset
    vw = 1/mu
    
  }else{
    converge = (iter < n_iter)
  }
  
  beta = t(beta)
  par = cbind(beta, tau, converge, iter)
  
  return(list("par"=par, "w"=w, "vw"=vw, "mu"=mu))
}