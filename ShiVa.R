library(glmnet)
library(ape)
library(MASS)
library(psych)

generate_design_matrix = function(tree, type = "simpX", alpha = 0) 
{
    stopifnot(is.ultrametric(tree))
    stopifnot(sum(1:length(tree$tip.label) %in% tree$edge[, 1]) == 
        0)
    nTips = length(tree$tip.label)
    rNode = nTips + 1
    nEdges = Nedge(tree)
    g = graph.edgelist(tree$edge, directed = TRUE)
    X = matrix(0, nTips, nEdges)
    root2tip = get.shortest.paths(g, rNode, to = 1:nTips, mode = "out", 
        output = "epath")$epath
    stopifnot(all(lapply(root2tip, length) > 0))
    Tval = sum(tree$edge.length[root2tip[[1]]])
    if (type == "orgX") {
        for (i in 1:nTips) {
            lvec = c(0, tree$edge.length[root2tip[[i]]])
            timeVec = Tval - cumsum(lvec)
            timeVec = timeVec[1:length(timeVec) - 1]
            X[i, root2tip[[i]]] = 1 - exp(-alpha * timeVec)
        }
    }
    else if (type == "simpX") {
        for (i in 1:nTips) {
            X[i, root2tip[[i]]] = 1
        }
    }
    else stop("Undefined design matrix type")
    return(X)
}

updatevcv = function(ape.tre, new.rates) {
  vv=vcv(updatetree(ape.tre, new.rates))
  return(vv)
}

updatetree = function(ape.tre, new.rates){
  ape.tre$edge.length=ape.tre$edge.length*new.rates
  return(ape.tre)
}

OU.vcv = function(phy,theta) {
  #theta is STRENGTH OF SELECTION
  times = updatevcv(phy,1)
  V = matrix(nrow=nrow(times),ncol=ncol(times))
  for (i in 1:nrow(V)) {
    for (j in 1:ncol(V)) {
      V[i,j] = 1/(2*theta)*exp(-2*theta*(times[i,i]-times[i,j]))*(1-exp(-2*theta*times[i,j]))
    }
  }
  return(V)
}

plot_shifts = function(tree, s.c,s.c2=NULL,show.tip.label=FALSE,main = NULL,xlab = NULL,no.margin=FALSE){
  #par(mfrow=c(1,1))
  plot(tree,show.tip.label=show.tip.label,main=main,no.margin=no.margin)
  title(xlab=xlab)
  edgelabels(s.c,s.c,cex=1.5)
  if(length(s.c2)>0)edgelabels(s.c2,s.c2,bg='pink',cex=0.8)
}

soft_thresholding = function(z, lambda){
  if(z > lambda){
    new_value = z-lambda
  }else if(z < -lambda){
    new_value = z+lambda
  }else{
    new_value = 0
  }
  return(new_value)
}


 get_mean_var_shifts = function(Y, tree, alpha, lambda1, lambda2, max.steps=1000, t = 0.01, penalty = 'L1', thres = 0.01,sigma2= NULL,measurement_error=FALSE){
   X = generate_design_matrix(tree,type='simpX')
   internal_list = (1:ncol(X))[colSums(X)>1 & colSums(X)!=(nrow(X)-1)]

   if(alpha == 0){
    V = vcv(tree)
   }else{
    V = OU.vcv(tree,alpha)
  }
   n = nrow(X)
   p = ncol(X)
   if(is.null(sigma2)){
     sigma2_0 = 1
   }else{
     sigma2_0 = sigma2
   }
   sigma2_error = if(measurement_error) 1 else 0
   tau_0 = log(sigma2_0)
   tau_error = log(sigma2_error)
   tb = node.depth.edgelength(tree)[tree$edge[,1]]
   q = exp(-2*alpha*max(node.depth.edgelength(tree)))/(2*alpha)*(1-exp(2*alpha*tb))

   b0 = 0
   beta = rep(0,p)
   gamma = rep(0,p)
   r = Y
   #tlist = rep(t,p)
    loss = Inf
#    best_loss = Inf
#    best_sigma2 = sigma2_0
#    best_sigma2_error = sigma2_error
#    best_gamma = gamma
#    best_beta = beta
#    best_b0 = b0
    Sigma = V*sigma2_0 +  (X%*%diag(gamma)%*%t(X))*V+ X%*%diag(gamma*q)%*%t(X)
    InvSigma = solve(Sigma)
    loss_list = c()

   for(s in 1:max.steps){
     
     #print(paste('step:',s,sep=''))
     #print(paste('sigma2:',sigma2_0,sep=''))
     
     last_tau = tau_0
     last_gamma = gamma
     last_beta = beta
     last_b0 = b0
     last_loss = loss
     last_tau_error = tau_error
     
     # update beta
     #for(k in 1:p){
     #  result = update_step_beta(beta[k], X[,k], InvSigma, r, lambda1, t,penalty)
     #  r = r - X[,k]*(result[['beta_k']] - beta[k]) 
     #  beta[k] = result[['beta_k']]
     #}
     svd_result = svd(InvSigma)
     SqrtInvSigma = svd_result$u %*% diag(sqrt(svd_result$d))  %*% t(svd_result$v)
     YY = SqrtInvSigma %*% Y
     XX = SqrtInvSigma %*% X
     beta = as.vector(glmnet(XX,YY,'gaussian',lambda=lambda1,intercept = FALSE)$beta)
     r = Y - X %*% beta
     #b0 = 0
     # update b0
     tmp = t(rep(1,n)) %*% InvSigma
     b0 = as.vector((tmp %*% (r+b0))/(tmp %*% rep(1,n)))
     r = r - (b0-last_b0)
     
     # update gamma
     for(k in internal_list){
       #print(k)
       result = update_step_gamma(gamma[k],X[,k],Sigma,r,lambda2,t,penalty,V,q[k])
       #print(min(diag(Sigma)))
       gamma[k] = result[['gamma_k']]
       Sigma = result[['Sigma']]
       #tlist[k] = result[['t_k']]

     } 

     # update sigma0
     if(is.null(sigma2)){
        InvSigma = solve(Sigma,tol=exp(-100))
        InvSigmar = InvSigma%*%r
        tau_delta =  (-1/2*t(InvSigmar)%*%V%*%InvSigmar + 1/2*tr(V%*%InvSigma))*exp(tau_0)
        #st = abs(2/(t(InvSigmar) %*% V %*% InvSigma %*% V %*% r-tr(V%*%InvSigma%*%V%*%InvSigma)))
        #print(paste("step size:",st,sep=""))
        #while((last_sigma2 - st*sigma2_delta)[1]*V[1,1]+min(gamma)<0){

        #print(paste("value:",min(diag(Sigma))  - st*sigma2_delta[1],sep=""))
        tau_0 = (last_tau - t*tau_delta)[1]
        Sigma = Sigma + (exp(tau_0)-exp(last_tau)) * V
        InvSigma = solve(Sigma,tol=exp(-100))
        InvSigmar = InvSigma%*%r
     }

    # update sigma_error
    if(measurement_error){  
      tau_error_delta = (-1/2*t(InvSigmar)%*%InvSigmar + 1/2*tr(InvSigma))*exp(tau_error)
      tau_error = (last_tau_error - t*tau_error_delta)[1]
      Sigma = Sigma + (exp(tau_error) - exp(last_tau_error)) * diag(1,n)
      InvSigma = solve(Sigma)
    }

      if (any(!is.finite(Sigma))) {
      stop("Sigma has NaN or Inf — likely from unstable gamma or bad matrix ops")
      }
     #print(paste("beta:",paste(beta,collapse=',')))
     #print(paste("gamma:",paste(gamma,collapse=',')))
     #print(paste('b0:',b0))
      max_update = max(abs(c(last_beta-beta,last_gamma-gamma)))
     #print(paste('max_update:',max_update,sep=' '))
      loglik = - 1/2*(n*log(2*pi) + t(r) %*% (InvSigma %*% r) +  determinant(Sigma)$modulu)
      loss = -loglik
      if(lambda1 != Inf) loss = loss + sum(lambda1*abs(beta))
      if(lambda2 != Inf) loss = loss + sum(lambda2*abs(gamma))
      #print(paste('sigma2:',exp(tau_0),sep=''))
      #print(paste('sigma2_error:',exp(tau_error),sep=''))
      #print(paste("det(Sigma):",determinant(Sigma)$modulu,sep=""))
      #print(paste("loss:",loss,sep=""))
#      if(loss < best_loss){
#          best_loss = loss
#          best_sigma2 = sigma2_0
#          best_sigma2_error = sigma2_error
#          best_gamma = gamma
#          best_beta = beta
#          best_b0 = b0
#      }
      loss_list = c(loss_list,loss)
      if(((abs(loglik)!=Inf)&(abs(loss-last_loss)<thres))|loglik==-Inf) {
        print(paste(s,'steps to converge',sep=' '))
        break
      }
    }

    sigma2_0 = exp(tau_0)
    sigma2_error = exp(tau_error)
#    sigma2_0 = best_sigma2
#    sigma2_error = best_sigma2_error
#    gamma = best_gamma 
#    beta = best_beta
#    b0 = best_b0
   
   result = list('sv_mean' = (1:ncol(X))[beta!=0],'sv_var' = (1:ncol(X))[gamma!=0],
  'gamma'=gamma,'beta'=beta,'sigma2'=sigma2_0,'b0'=b0)
  if(measurement_error){
    result[['sigma2_error']] = sigma2_error
  }
  #plot(1:length(loss_list),loss_list,type="l")
  return(result)
}

update_step_gamma <- function(gamma_k, X_k, Sigma, r, lambda2, t, penalty, V, q_k) {
  XXT <- X_k %*% t(X_k)
  M_k <- XXT * V

  # Use Cholesky-based inversion (faster and more stable)
  chol_result <- tryCatch(chol(Sigma), error = function(e) NULL)
  if (is.null(chol_result)) {
    warning("Cholesky failed: skipping gamma update.")
    return(list(gamma_k = gamma_k, Sigma = Sigma))
  }
  InvSigma <- chol2inv(chol_result)

  InvSigmar <- InvSigma %*% r

  # Compute derivative
  delta <- -0.5 * (t(InvSigmar) %*% M_k %*% InvSigmar) -
           0.5 * q_k * (t(X_k) %*% InvSigmar)^2 +
           0.5 * tr(M_k %*% InvSigma) +
           0.5 * q_k * (t(X_k) %*% InvSigma %*% X_k)
  delta <- as.numeric(delta)

  # Skip update if gradient is very small (save time)
  if (abs(delta) < 1e-6 || !is.finite(delta)) {
    return(list(gamma_k = gamma_k, Sigma = Sigma))
  }

  # Initialize update
  z <- gamma_k - t * delta
  if (penalty == "L1") {
    new_value <- soft_thresholding(z, lambda2)
  } else {
    new_value <- z
  }

  # Clamp gamma for stability
  new_value <- max(0, min(new_value, 10))

  # Limit to 5 backtracking steps max
  for (i in 1:5) {
    diag_shift <- (new_value - gamma_k) * (q_k + V[1, 1])
    if (all(diag(Sigma)[X_k != 0] + diag_shift > 0)) {
      break
    }
    t <- t * 0.5
    z <- gamma_k - t * delta
    if (penalty == "L1") {
      new_value <- soft_thresholding(z, lambda2)
    } else {
      new_value <- z
    }
    new_value <- max(0, min(new_value, 10))
  }

  # Final Sigma update
  Sigma <- Sigma + (new_value - gamma_k) * (M_k + q_k * XXT)

  return(list(gamma_k = new_value, Sigma = Sigma))
}

backward_correction = function(tree,Y,alpha,sv_mean,sv_var,criterion='BIC',original_model=NULL,measurement_error = FALSE,
                               max.num.shifts = Inf){

  if(is.null(original_model)){
    OModel = fit_OU_mean_var(tree,Y,alpha,sv_mean,sv_var,measurement_error = measurement_error,
                              max.num.shifts = max.num.shifts)
  }else{
    original_model = OModel
  }
  original_score = OModel[[criterion]][[1]]
  best_score = original_score
  n1 = length(sv_mean)
  n2 = length(sv_var)
  if(n1+n2==0){
    return(OModel)
  }
  after_score_list = rep(0,n1+n2)
  for(i in 1:(n1+n2)){
    if(i <= length(sv_mean)){
      new_Model = fit_OU_mean_var(tree,Y,alpha,setdiff(sv_mean,c(sv_mean[i])),sv_var,measurement_error = measurement_error,
                              max.num.shifts = max.num.shifts)
      after_score_list[i] = new_Model[[criterion]][[1]]
    } else{
      new_Model = fit_OU_mean_var(tree,Y,alpha,sv_mean,setdiff(sv_var,c(sv_var[i-n1])),measurement_error = measurement_error,
                              max.num.shifts = max.num.shifts)
      after_score_list[i] = new_Model[[criterion]][[1]]      
    }
    #print(paste("index:",i,"score:",new_Model[[criterion]][[1]],sep=' '))
    if(after_score_list[i]<best_score){
      OModel = new_Model
      best_score = after_score_list[i]
    }
  }
  index_list = order(after_score_list)[sort(after_score_list)<original_score]
  #print(paste('index_list:',paste(index_list,collapse = ';'),sep=''))
  remove_list = c(index_list[1])
  for(i in index_list[2:length(index_list)]){
    new_Model = fit_OU_mean_var(tree,Y,alpha,
                                setdiff(sv_mean,sv_mean[c(i,remove_list)[c(i,remove_list)<=n1]]),
                                setdiff(sv_var,sv_var[c(i,remove_list)[c(i,remove_list)>n1]-n1]),
                                measurement_error = measurement_error,
                                max.num.shifts = max.num.shifts)
   if(new_Model[[criterion]][[1]]<best_score){
      OModel = new_Model
      best_score = new_Model[[criterion]][[1]]
      remove_list = c(remove_list,i)
    }
  }
  return(OModel)
  }

get_mean_var_shifts_model_selection <- function(Y, tree, alpha, t = 0.01,
                                                lambda1_list = NULL, lambda2_list,
                                                criterion = "BIC", max.steps = 300,
                                                nfolds = 8, top_k = 10,
                                                measurement_error = FALSE,
                                                lambda.type = "lambda.1se",
                                                max.num.shifts = Inf
                                                ) {

  top_k = top_k
  X <- generate_design_matrix(tree, 'simpX')
  Y <- as.vector(Y)
  V <- if (alpha == 0) vcv(tree) else OU.vcv(tree, alpha)

  n <- nrow(X)
  p <- ncol(X)
  tb <- node.depth.edgelength(tree)[tree$edge[,1]]
  q <- exp(-2 * alpha * max(node.depth.edgelength(tree))) / (2 * alpha) * (1 - exp(2 * alpha * tb))

  max.num.shifts = min(max.num.shifts, p)
  score_summary <- data.frame(matrix(nrow = 0, ncol = 8))
  names(score_summary) <- c('lambda1', 'lambda2', 'sv_mean', 'sv_var',
                            'loglik', 'BIC', 'mBIC', 'pBIC')

  Y1 <- Y
  best_score <- Inf
  model_counter <- 1

  for (lambda2 in lambda2_list) {
    cat("====== Model Selection Round", model_counter, "======\n")
    cat("Trying lambda2 =", lambda2, "...\n")

    ret_pre <- get_mean_var_shifts(Y, tree, alpha, Inf, lambda2, max.steps = max.steps)
    Sigma <- V * ret_pre$sigma2 +
             (X %*% diag(ret_pre$gamma) %*% t(X)) * V +
             X %*% diag(ret_pre$gamma * q) %*% t(X)
    InvSigma <- solve(Sigma, tol = exp(-100))
    svd_result <- svd(InvSigma)
    SqrtInvSigma <- svd_result$u %*% diag(sqrt(svd_result$d)) %*% t(svd_result$v)

    YY <- SqrtInvSigma %*% Y
    XX <- SqrtInvSigma %*% X
    lambda1 <- cv.glmnet(XX, YY, lambda = lambda1_list, intercept = FALSE, nfolds = nfolds)[[lambda.type]]
    cat("Selected lambda1 from CV:", lambda1, "\n")

    ret <- get_mean_var_shifts(Y, tree, alpha, lambda1, lambda2,
                               t = t, max.steps = max.steps, measurement_error = measurement_error)

    cat("  sv_mean =", if (length(ret$sv_mean)) paste(ret$sv_mean, collapse = ",") else "none", "\n")
    cat("  sv_var  =", if (length(ret$sv_var)) paste(ret$sv_var, collapse = ",") else "none", "\n")

    OModel <- fit_OU_mean_var(tree, Y1, alpha, ret$sv_mean, ret$sv_var,
                              t = t, measurement_error = measurement_error,
                                                max.num.shifts = max.num.shifts)
    cat("  log-likelihood =", OModel$loglik, "\n\n")

    score_summary[nrow(score_summary) + 1, ] <- c(lambda1, lambda2,
                                                  paste(ret$sv_mean, collapse = ';'),
                                                  paste(ret$sv_var, collapse = ';'),
                                                  OModel$loglik,
                                                  OModel$BIC,
                                                  OModel$mBIC,
                                                  OModel$pBIC)
    if (OModel[[criterion]][[1]] < best_score) {
      best_score <- OModel[[criterion]][[1]]
      best_Model <- OModel
      best_Model$lambda1 <- lambda1
      best_Model$lambda2 <- lambda2
    }

    model_counter <- model_counter + 1
  }

  score_summary <- score_summary[!duplicated(score_summary[, c('sv_mean', 'sv_var')]), ]
  score_summary <- cbind(score_summary, NA, NA, NA, NA, NA, NA)
  names(score_summary)[9:14] <- c('sv_mean_corrected', 'sv_var_corrected',
                                  'loglik_corrected', 'BIC_corrected',
                                  'mBIC_corrected', 'pBIC_corrected')

  cat("\n====== Backward Correction (Top", top_k, ") ======\n")
  for (i in 1:min(top_k, nrow(score_summary))) {
    ind <- order(as.numeric(score_summary[[criterion]]))[i]
    sv_mean <- as.numeric(strsplit(score_summary$sv_mean[ind], ';')[[1]])
    sv_var <- as.numeric(strsplit(score_summary$sv_var[ind], ';')[[1]])

    cat("Correcting model", i, "with sv_mean =", paste(sv_mean, collapse = ","),
        "sv_var =", paste(sv_var, collapse = ","), "...\n")

    OModel <- backward_correction(tree, Y1, alpha, sv_mean, sv_var, measurement_error = measurement_error,
                                                max.num.shifts = max.num.shifts)

    score_summary[ind, 9:14] <- c(paste(which(OModel$beta != 0), collapse = ';'),
                                  paste(which(OModel$gamma != 0), collapse = ';'),
                                  OModel$loglik,
                                  OModel$BIC,
                                  OModel$mBIC,
                                  OModel$pBIC)

    if (OModel[[criterion]][[1]] < best_score) {
      best_score <- OModel[[criterion]][[1]]
      best_Model <- OModel
      best_Model$lambda1 <- score_summary$lambda1[ind]
      best_Model$lambda2 <- score_summary$lambda2[ind]
    }
  }

  cat("====== Selection Finished. Best", criterion, "=", best_score, "======\n")
  return(list('best_model' = best_Model, 'score_summary' = score_summary))
}

fit_OU_mean_var = function(tree, Y, alpha, sv_mean, sv_var,
                            max.steps = 1000, t = 0.01, thres = 0.01,
                            measurement_error = FALSE, max.num.shifts = Inf) {
    X <- generate_design_matrix(tree, "simpX")
    n <- nrow(X)
    p <- ncol(X)
    
    max.num.shifts = min(max.num.shifts, p)
    if (length(sv_mean) > max.num.shifts || length(sv_var) > max.num.shifts) {
      return(list(loglik = -Inf, BIC = Inf, mBIC = Inf, pBIC = Inf))
    }

    if (alpha == 0) {
        V <- vcv(tree)
    } else {
        V <- OU.vcv(tree, alpha)
    }
    
    tb <- node.depth.edgelength(tree)[tree$edge[, 1]]
    q <- exp(-2 * alpha * max(node.depth.edgelength(tree))) / (2 * alpha) * (1 - exp(2 * alpha * tb))
    
    par0 <- c(rep(0, length(sv_mean)), rep(0, length(sv_var)), 0, log(1))
    if (measurement_error) {
        par0 <- c(par0, log(1))
    }
    
    nll_fn <- function(par) {
        m <- length(sv_mean)
        v <- length(sv_var)
        beta <- rep(0, p)
        gamma <- rep(0, p)
        
        if (m > 0) beta[sv_mean] <- par[1:m]
        if (v > 0) gamma[sv_var] <- par[(m+1):(m+v)]
        
        b0 <- par[m + v + 1]
        log_sigma2 <- par[m + v + 2]
        sigma2 <- exp(log_sigma2)
        sigma2_error <- if (measurement_error) exp(par[m + v + 3]) else 0
        
        mu <- X %*% beta + b0
        r <- Y - mu
        
        Sigma <- V * sigma2 +
            (X %*% diag(gamma) %*% t(X)) * V +
            X %*% diag(gamma * q) %*% t(X) +
            sigma2_error * diag(n)
        
        cholSigma <- tryCatch(chol(Sigma), error = function(e) return(NULL))
        if (is.null(cholSigma)) return(1e10)
        
        logdet <- 2 * sum(log(diag(cholSigma)))
        InvSigma_r <- backsolve(cholSigma, forwardsolve(t(cholSigma), r))
        0.5 * (n * log(2 * pi) + logdet + sum(r * InvSigma_r))
    }
    
    opt_result <- optim(par = par0, fn = nll_fn, method = "L-BFGS-B", control = list(maxit = max.steps))
    
    opt_par <- opt_result$par
    m <- length(sv_mean)
    v <- length(sv_var)
    
    beta <- rep(0, p)
    gamma <- rep(0, p)
    if (m > 0) beta[sv_mean] <- opt_par[1:m]
    if (v > 0) gamma[sv_var] <- opt_par[(m+1):(m+v)]
    
    b0 <- opt_par[m + v + 1]
    sigma2 <- exp(opt_par[m + v + 2])
    sigma2_error <- if (measurement_error) exp(opt_par[m + v + 3]) else 0
    
    mu <- X %*% beta + b0
    r <- Y - mu
    Sigma <- V * sigma2 +
        (X %*% diag(gamma) %*% t(X)) * V +
        X %*% diag(gamma * q) %*% t(X) +
        sigma2_error * diag(n)
    InvSigma <- solve(Sigma, tol = exp(-100))
    loglik <- -nll_fn(opt_par)
    
    # Criteria
    k_beta <- sum(beta != 0)
    k_gamma <- sum(gamma != 0)
    BIC <- -2 * loglik + log(n) * (2 * k_beta + 2 * k_gamma + 3)
    
    logn_mean <- if (length(sv_mean) > 1) sum(log(colSums(X[, sv_mean, drop = FALSE])))
    else if (length(sv_mean) == 1) log(sum(X[, sv_mean]))
    else 0
    
    logn_var <- if (length(sv_var) > 1) sum(log(colSums(X[, sv_var, drop = FALSE])))
    else if (length(sv_var) == 1) log(sum(X[, sv_var]))
    else 0
    
    active <- c(sv_mean, sv_var)
    logn0 <- if (length(active) == 0) log(n)
    else if (length(active) == 1) log(n - sum(X[, active] > 0))
    else {
        z <- rowSums(X[, active, drop = FALSE]) > 0
        if (sum(z) < n) log(n - sum(z)) else 0
    }
    
    mBIC <- -2 * loglik + (2 * k_beta + 2 * k_gamma - 1) * log(n) + logn_mean + logn_var + logn0
    
    design_matrix <- cbind(1, X[, sv_mean, drop = FALSE])
    Xt_Si_X <- t(design_matrix) %*% InvSigma %*% design_matrix
    logdet_proj <- as.numeric(determinant(Xt_Si_X, logarithm = TRUE)$modulus)
    
    pBIC <- -2 * loglik + 2 * (k_beta + k_gamma) * log(2 * n - 3) + 2 * log(n) + logdet_proj
    
    result <- list(
        tree = tree,
        shifts_mean = which(beta!=0),
        shifts_var = which(gamma!=0),
        alpha = alpha,
        beta = beta,
        gamma = gamma,
        sigma2 = sigma2,
        b0 = b0,
        loglik = loglik,
        BIC = BIC,
        mBIC = mBIC,
        pBIC = pBIC,
        fitted.values = mu,
        Sigma = Sigma
    )
    
    if (measurement_error) {
        result$sigma2_error <- sigma2_error
    }
    
    return(result)
}


get_test_data = function(tree, Sigma, n_test, alpha,beta){
  X = generate_design_matrix(tree,type='simpX')

  eps = mvrnorm(n = 1, mu = rep(0,nrow(X)), Sigma = Sigma)
  ret = as.data.frame(t(X%*%beta+eps))

  for(i in 2:n_test){
      eps = mvrnorm(n = 1, mu = rep(0,nrow(X)), Sigma = Sigma)
      Y = X %*% beta + eps
      ret[nrow(ret)+1,] = Y
  }
  return(as.matrix(ret))

}

get_prediction_likelihood = function(Y_pred,Sigma,alpha,test_data,sigma2 = 2*alpha){

    #Y_pred = OModel$fitted.values
    #Sigma = OModel$Sigma
    InvSigma = solve(Sigma)
    n = length(Y_pred)

    cmp_likelihood = function(target,prediction,Sigma,InvSigma){
        return(-1/2*n*log(2*pi)-1/2*t(target-prediction)%*%(InvSigma %*% (target-prediction)) - 1/2* determinant(Sigma)$modulu)
            }

    loglik = apply(test_data,1,cmp_likelihood,prediction=Y_pred,Sigma = Sigma,InvSigma = InvSigma)

    return(mean(loglik))

}
