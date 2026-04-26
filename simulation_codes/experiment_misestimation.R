library(ShiVa)         # CRAN release of our method
library(l1ou)
source('../method_helpers/l1ou_fixed.R') # adds estimate_shift_configuration_known_alpha and
                       # patches a few internals so l1ou builds against
                       # current dependency versions (l1ou is no longer
                       # actively maintained on CRAN).
library(PhylogeneticEM)
library(PCMBase)
library(PCMBaseCpp)
library(PCMFit)
library(data.table)
library(MASS)
library(ape)
library(phytools)
library(phylolm)

# ---------------------------------------------------------------------------
# Simulation utilities (not part of the ShiVa package interface)
# ---------------------------------------------------------------------------

updatetree = function(ape.tre, new.rates) {
  ape.tre$edge.length = ape.tre$edge.length * new.rates
  ape.tre
}
updatevcv = function(ape.tre, new.rates) vcv(updatetree(ape.tre, new.rates))

OU.vcv = function(phy, theta) {
  times = updatevcv(phy, 1)
  V = matrix(nrow = nrow(times), ncol = ncol(times))
  for (i in 1:nrow(V)) for (j in 1:ncol(V))
    V[i,j] = 1/(2*theta) * exp(-2*theta*(times[i,i] - times[i,j])) *
             (1 - exp(-2*theta*times[i,j]))
  V
}

get_test_data = function(tree, Sigma, n_test, alpha, beta) {
  X = generate_design_matrix(tree, type = 'simpX')
  ret = matrix(0, nrow = n_test, ncol = nrow(X))
  for (i in 1:n_test) {
    eps = mvrnorm(n = 1, mu = rep(0, nrow(X)), Sigma = Sigma)
    ret[i, ] = X %*% beta + eps
  }
  ret
}

get_prediction_likelihood = function(Y_pred, Sigma, alpha, test_data,
                                      sigma2 = 2 * alpha) {
  InvSigma = solve(Sigma)
  n = length(Y_pred)
  cmp_likelihood = function(target, prediction, Sigma, InvSigma) {
    -1/2 * n * log(2*pi) -
      1/2 * t(target - prediction) %*% (InvSigma %*% (target - prediction)) -
      1/2 * determinant(Sigma)$modulu
  }
  mean(apply(test_data, 1, cmp_likelihood,
             prediction = Y_pred, Sigma = Sigma, InvSigma = InvSigma))
}

# ---------------------------------------------------------------------------
# MISE helpers
# ---------------------------------------------------------------------------

get_effective_sigma2_per_branch = function(tree, sigma2_0, gamma) {
  n_nodes   = max(tree$edge)
  root_node = setdiff(tree$edge[,1], tree$edge[,2])[1]
  cumgamma  = rep(0, n_nodes)
  effective_sigma2 = rep(sigma2_0, nrow(tree$edge))
  queue = root_node
  while (length(queue) > 0) {
    cur = queue[1]; queue = queue[-1]
    child_edges = which(tree$edge[,1] == cur)
    for (ei in child_edges) {
      child = tree$edge[ei, 2]
      cumgamma[child]      = cumgamma[cur] + gamma[ei]
      effective_sigma2[ei] = sigma2_0 + cumgamma[child]
      if (child > length(tree$tip.label)) queue = c(queue, child)
    }
  }
  return(effective_sigma2)
}

compute_MISE = function(tree, sigma2_0_true, gamma_true, sigma2_0_hat, gamma_hat) {
  s2_true = get_effective_sigma2_per_branch(tree, sigma2_0_true, gamma_true)
  s2_hat  = get_effective_sigma2_per_branch(tree, sigma2_0_hat,  gamma_hat)
  sum(tree$edge.length * (s2_true - s2_hat)^2)
}

get_effective_theta_per_branch = function(tree, theta_0, beta) {
  n_nodes   = max(tree$edge)
  root_node = setdiff(tree$edge[,1], tree$edge[,2])[1]
  cumbeta   = rep(0, n_nodes)
  effective_theta = rep(theta_0, nrow(tree$edge))
  queue = root_node
  while (length(queue) > 0) {
    cur = queue[1]; queue = queue[-1]
    child_edges = which(tree$edge[,1] == cur)
    for (ei in child_edges) {
      child = tree$edge[ei, 2]
      cumbeta[child]       = cumbeta[cur] + beta[ei]
      effective_theta[ei]  = theta_0 + cumbeta[child]
      if (child > length(tree$tip.label)) queue = c(queue, child)
    }
  }
  return(effective_theta)
}

compute_MISE_mu = function(tree, theta_0_true, beta_true, theta_0_hat, beta_hat) {
  t_true = get_effective_theta_per_branch(tree, theta_0_true, beta_true)
  t_hat  = get_effective_theta_per_branch(tree, theta_0_hat,  beta_hat)
  sum(tree$edge.length * (t_true - t_hat)^2)
}

get_PCMFit_sigma2_per_branch = function(tree, bestFit) {
  tryCatch({
    model        = bestFit$inferredModel
    n_regimes    = length(model) - 1
    regime_nodes = as.numeric(bestFit$inferredRegimeNodes)

    sigma2_per_regime = sapply(1:n_regimes, function(r) {
      sx = model[[r + 1]]$Sigma_x
      if (!is.null(sx) && length(sx) > 0) sx[1, 1, 1]^2
      else NA_real_
    })

    root_node    = regime_nodes[1]
    change_nodes = if (n_regimes >= 2) regime_nodes[-1] else c()

    n_nodes     = max(tree$edge)
    node_regime = rep(1L, n_nodes)
    for (r in seq_along(change_nodes)) node_regime[change_nodes[r]] = r + 1L

    edge_sigma2     = rep(sigma2_per_regime[1], nrow(tree$edge))
    node_cur_regime = rep(1L, n_nodes)
    queue = root_node

    while (length(queue) > 0) {
      cur = queue[1]; queue = queue[-1]
      cur_regime  = node_cur_regime[cur]
      child_edges = which(tree$edge[,1] == cur)
      for (ei in child_edges) {
        child        = tree$edge[ei, 2]
        child_regime = if (child %in% change_nodes) node_regime[child] else cur_regime
        node_cur_regime[child] = child_regime
        edge_sigma2[ei]        = sigma2_per_regime[child_regime]
        if (child > length(tree$tip.label)) queue = c(queue, child)
      }
    }
    edge_sigma2
  }, error = function(e) {
    cat("Warning: PCMFit sigma2 extraction failed:", conditionMessage(e), "\n")
    rep(NA_real_, nrow(tree$edge))
  })
}

get_PCMFit_theta_per_branch = function(tree, bestFit) {
  tryCatch({
    model        = bestFit$inferredModel
    n_regimes    = length(model) - 1
    regime_nodes = as.numeric(bestFit$inferredRegimeNodes)

    theta_per_regime = sapply(1:n_regimes, function(r) {
      th = model[[r + 1]]$Theta
      if (!is.null(th) && length(th) > 0) th[1, 1]
      else NA_real_
    })

    root_node    = regime_nodes[1]
    change_nodes = if (n_regimes >= 2) regime_nodes[-1] else c()

    n_nodes     = max(tree$edge)
    node_regime = rep(1L, n_nodes)
    for (r in seq_along(change_nodes)) node_regime[change_nodes[r]] = r + 1L

    edge_theta      = rep(theta_per_regime[1], nrow(tree$edge))
    node_cur_regime = rep(1L, n_nodes)
    queue = root_node

    while (length(queue) > 0) {
      cur = queue[1]; queue = queue[-1]
      cur_regime  = node_cur_regime[cur]
      child_edges = which(tree$edge[,1] == cur)
      for (ei in child_edges) {
        child        = tree$edge[ei, 2]
        child_regime = if (child %in% change_nodes) node_regime[child] else cur_regime
        node_cur_regime[child] = child_regime
        edge_theta[ei]         = theta_per_regime[child_regime]
        if (child > length(tree$tip.label)) queue = c(queue, child)
      }
    }
    edge_theta
  }, error = function(e) {
    cat("Warning: PCMFit theta extraction failed:", conditionMessage(e), "\n")
    rep(NA_real_, nrow(tree$edge))
  })
}

# ---------------------------------------------------------------------------
# PCMFit model generator
# ---------------------------------------------------------------------------

generatePCMModelsFunction <- function() {
  set.seed(4, kind = "Mersenne-Twister", normal.kind = "Inversion")
  PCMGenerateModelTypes()
  source('DefineParameterLimits.R', local=FALSE)
}

# ---------------------------------------------------------------------------
# Tree and parameter setup
# ---------------------------------------------------------------------------

data('lizard.tree')
tree1 = lizard.tree
tree1$edge.length = tree1$edge.length / max(node.depth.edgelength(tree1))
X = generate_design_matrix(tree1, type='simpX')

args = commandArgs(T)
k = as.numeric(args[1])
runfile = read.csv('runfile_misestimation.csv')

# True parameters (fixed)
alpha    = 1
sigma2_0 = 2

params           = as.character(runfile[ceiling(k/20), ])
true_shifts_mean = as.numeric(strsplit(params[2], ',')[[1]])
size_mean        = as.numeric(strsplit(params[3], ',')[[1]])
true_shifts_var  = as.numeric(strsplit(params[4], ',')[[1]])
size_var         = as.numeric(strsplit(params[5], ',')[[1]])
alpha_hat        = as.numeric(strsplit(params[6], ',')[[1]])

cat(sprintf("alpha_hat: %.4f | true_shifts_mean: %s | size_mean: %s | true_shifts_var: %s | size_var: %s\n",
    alpha_hat,
    paste(true_shifts_mean, collapse=','), paste(size_mean, collapse=','),
    paste(true_shifts_var,  collapse=','), paste(size_var,  collapse=',')))

# Data generating covariance (with TRUE alpha)
V  = OU.vcv(tree1, alpha)
tb = node.depth.edgelength(tree1)[tree1$edge[,1]]
q  = exp(-2*alpha*max(node.depth.edgelength(tree1)))/(2*alpha)*(1-exp(2*alpha*tb))

gamma = rep(0, ncol(X))
gamma[true_shifts_var] = size_var
Sigma = V*sigma2_0 + X%*%diag(gamma)%*%t(X)*V + X%*%diag(gamma*q)%*%t(X)

beta = rep(0, ncol(X))
beta[true_shifts_mean] = size_mean
test_data = get_test_data(tree1, Sigma, 1000, alpha, beta)

# ---------------------------------------------------------------------------
# Result tables
# ---------------------------------------------------------------------------

r_phyloEM = data.frame(matrix(nrow=0, ncol=ncol(X)))
r_l1ou    = list("pBIC"=data.frame(matrix(nrow=0, ncol=ncol(X))),
                 "BIC" =data.frame(matrix(nrow=0, ncol=ncol(X))))
r_PCMFit  = data.frame(matrix(nrow=0, ncol=ncol(X)))
r_ShiVa   = list("shifts_mean"=data.frame(matrix(nrow=0, ncol=ncol(X))),
                 "shifts_var" =data.frame(matrix(nrow=0, ncol=ncol(X))))

loglik_table = data.frame(matrix(nrow=0, ncol=6))
names(loglik_table) = c('ShiVa','l1ou+pBIC','l1ou+BIC','phyloEM','PCMFit','true')

computing_time_table = data.frame(matrix(nrow=0, ncol=5))
names(computing_time_table) = c('ShiVa','l1ou+pBIC','l1ou+BIC','phyloEM','PCMFit')

mise_table = data.frame(matrix(nrow=0, ncol=6))
names(mise_table) = c('ShiVa','l1ou+pBIC','l1ou+BIC','phyloEM','PCMFit','true')

mise_mu_table = data.frame(matrix(nrow=0, ncol=6))
names(mise_mu_table) = c('ShiVa','l1ou+pBIC','l1ou+BIC','phyloEM','PCMFit','true')

# True sigma2 / theta per branch (fixed)
s2_true_per_branch    = get_effective_sigma2_per_branch(tree1, sigma2_0, gamma)
theta_true_per_branch = get_effective_theta_per_branch(tree1, 0, beta)

# ---------------------------------------------------------------------------
# Simulation loop
# ---------------------------------------------------------------------------

m = 3
i = 0

while(i < m){
  if((i-1<0.1*m) & i>=0.1*m){print('finish: 10%')}
  if((i-1<0.2*m) & i>=0.2*m){print('finish: 20%')}
  if((i-1<0.3*m) & i>=0.3*m){print('finish: 30%')}
  if((i-1<0.4*m) & i>=0.4*m){print('finish: 40%')}
  if((i-1<0.5*m) & i>=0.5*m){print('finish: 50%')}
  if((i-1<0.6*m) & i>=0.6*m){print('finish: 60%')}
  if((i-1<0.7*m) & i>=0.7*m){print('finish: 70%')}
  if((i-1<0.8*m) & i>=0.8*m){print('finish: 80%')}
  if((i-1<0.9*m) & i>=0.9*m){print('finish: 90%')}

  tryCatch({
  eps = mvrnorm(n = 1, mu = rep(0, nrow(X)), Sigma = Sigma)
  Y   = drop(X%*%beta + eps)
  names(Y) = tree1$tip.label

  #### ShiVa — fits with alpha_hat #####
  t       = Sys.time()
  result1 = ShiVa(Y, tree1, alpha = alpha_hat, refit_alpha = FALSE,
                  lambda1_list = NULL, lambda2_list = exp(1:10-6),
                  nfolds = 5, top_k = 3,
                  measurement_error = FALSE, verbose = FALSE)
  ctime_row = c(Sys.time()-t)
  r_ShiVa$shifts_mean[nrow(r_ShiVa$shifts_mean)+1,] = as.numeric(result1$best_model$beta!=0)
  r_ShiVa$shifts_var[ nrow(r_ShiVa$shifts_var) +1,] = as.numeric(result1$best_model$gamma!=0)

  sv_mean = (1:ncol(X))[result1$best_model$beta!=0]
  sv_var  = (1:ncol(X))[result1$best_model$gamma!=0]

  # Refit with TRUE alpha for fair loglik and MISE evaluation
  OModel    = fit_OU_mean_var(tree1, Y, alpha, sv_mean, sv_var)
  loglik_row  = c(get_prediction_likelihood(OModel$fitted.values, OModel$Sigma, alpha, test_data))
  mise_row    = c(compute_MISE(tree1, sigma2_0, gamma, OModel$sigma2, OModel$gamma))
  mise_mu_row = c(compute_MISE_mu(tree1, 0, beta, OModel$b0, OModel$beta))

  #### l1ou — fits with alpha_hat #####
  adj_data  = adjust_data(tree1, Y)
  tree2     = adj_data$tree
  Y2        = adj_data$Y

  l1ou_model = estimate_shift_configuration(tree2, Y2, max.nShifts=20)
  opt        = l1ou_model$l1ou.options
  opt$use.saved.scores = FALSE

  t = Sys.time()
  opt$criterion = 'pBIC'
  l1ou_model1 = estimate_shift_configuration_known_alpha(tree2, Y2, alpha=alpha_hat, opt=opt)
  ctime_row   = c(ctime_row, Sys.time()-t)
  r_l1ou$pBIC[nrow(r_l1ou$pBIC)+1,] = as.numeric(1:ncol(X) %in% l1ou_model1$shift.configuration)
  # Refit with TRUE alpha
  OModel2   = fit_OU_mean_var(tree2, Y2, alpha, l1ou_model1$shift.configuration, c())
  loglik_row  = c(loglik_row, get_prediction_likelihood(OModel2$fitted.values, OModel2$Sigma, alpha, test_data))
  mise_row    = c(mise_row,    compute_MISE(tree1, sigma2_0, gamma, OModel2$sigma2, rep(0, ncol(X))))
  mise_mu_row = c(mise_mu_row, compute_MISE_mu(tree1, 0, beta, OModel2$b0, OModel2$beta))

  t = Sys.time()
  opt$criterion = 'BIC'
  l1ou_model2 = estimate_shift_configuration_known_alpha(tree2, Y2, alpha=alpha_hat, opt=opt)
  ctime_row   = c(ctime_row, Sys.time()-t)
  r_l1ou$BIC[nrow(r_l1ou$BIC)+1,] = as.numeric(1:ncol(X) %in% l1ou_model2$shift.configuration)
  OModel3   = fit_OU_mean_var(tree2, Y2, alpha, l1ou_model2$shift.configuration, c())
  loglik_row  = c(loglik_row, get_prediction_likelihood(OModel3$fitted.values, OModel3$Sigma, alpha, test_data))
  mise_row    = c(mise_row,    compute_MISE(tree1, sigma2_0, gamma, OModel3$sigma2, rep(0, ncol(X))))
  mise_mu_row = c(mise_mu_row, compute_MISE_mu(tree1, 0, beta, OModel3$b0, OModel3$beta))

  #### phyloEM — fits with alpha_hat #####
  t       = Sys.time()
  emmodel = PhyloEM(tree1, as.vector(Y), process='OU', K_max=20, alpha=alpha_hat)
  sv_em   = params_process(emmodel)$shift$edges
  ctime_row = c(ctime_row, Sys.time()-t)
  r_phyloEM[nrow(r_phyloEM)+1,] = as.numeric(1:ncol(X) %in% sv_em)
  OModel4   = fit_OU_mean_var(tree1, Y, alpha, sv_em, c())
  loglik_row  = c(loglik_row, get_prediction_likelihood(OModel4$fitted.values, OModel4$Sigma, alpha, test_data))
  mise_row    = c(mise_row,    compute_MISE(tree1, sigma2_0, gamma, OModel4$sigma2, rep(0, ncol(X))))
  mise_mu_row = c(mise_mu_row, compute_MISE_mu(tree1, 0, beta, OModel4$b0, OModel4$beta))

  #### PCMFit — estimates its own alpha #####
  t = Sys.time()
  fitMGPM_A_F_BC2_RR <- PCMFitMixed(
      X = t(Y), tree = tree1, metaIFun = PCMInfoCpp,
      generatePCMModelsFun = generatePCMModelsFunction,
      maxNumRoundRobins = 2, maxNumPartitionsInRoundRobins = 10, minCladeSizes = 2L,
      tableFitsPrev = NULL,
      doParallel = FALSE)
  bestFit   <- RetrieveBestFitScore(fitMGPM_A_F_BC2_RR)
  ctime_row = c(ctime_row, Sys.time()-t)
  r_PCMFit[nrow(r_PCMFit)+1,] = as.numeric(tree1$edge[,2] %in% as.numeric(bestFit$inferredRegimeNodes))
  loglik_row = c(loglik_row,
                 get_prediction_likelihood(PCMBase::PCMMean(tree1, bestFit$inferredModel),
                                           PCMBase::PCMVar(tree1,  bestFit$inferredModel),
                                           alpha, test_data))
  s2_pcmfit   = get_PCMFit_sigma2_per_branch(tree1, bestFit)
  mise_pcmfit = if (!any(is.na(s2_pcmfit))) {
    sum(tree1$edge.length * (s2_true_per_branch - s2_pcmfit)^2)
  } else NA_real_
  mise_row = c(mise_row, mise_pcmfit)
  theta_pcmfit   = get_PCMFit_theta_per_branch(tree1, bestFit)
  mise_mu_pcmfit = if (!any(is.na(theta_pcmfit))) {
    sum(tree1$edge.length * (theta_true_per_branch - theta_pcmfit)^2)
  } else NA_real_
  mise_mu_row = c(mise_mu_row, mise_mu_pcmfit)

  #### True model #####
  loglik_row  = c(loglik_row, get_prediction_likelihood(X%*%beta, Sigma, alpha, test_data))
  mise_row    = c(mise_row, 0)
  mise_mu_row = c(mise_mu_row, 0)

  loglik_table[nrow(loglik_table)+1,]                = loglik_row
  computing_time_table[nrow(computing_time_table)+1,] = ctime_row
  mise_table[nrow(mise_table)+1,]                    = mise_row
  mise_mu_table[nrow(mise_mu_table)+1,]              = mise_mu_row

  i = i+1

  if(i == 1){ctime = proc.time()[1]*100}

  result = list('ShiVa'        = r_ShiVa,
                'l1ou'         = r_l1ou,
                'phyloEM'      = r_phyloEM,
                'PCMFit'       = r_PCMFit,
                'loglik'       = loglik_table,
                'compute_time' = computing_time_table,
                'MISE'         = mise_table,
                'MISE_mu'      = mise_mu_table)
  saveRDS(result,
          paste('results/experiment4_',
                paste(true_shifts_mean, collapse='-'), '_',
                paste(size_mean,        collapse='-'), '_',
                paste(true_shifts_var,  collapse='-'), '_',
                paste(size_var,         collapse='-'), '_',
                alpha_hat, '_', ctime, '.rds', sep=""))

  }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

result = list('ShiVa'        = r_ShiVa,
              'l1ou'         = r_l1ou,
              'phyloEM'      = r_phyloEM,
              'PCMFit'       = r_PCMFit,
              'loglik'       = loglik_table,
              'compute_time' = computing_time_table,
              'MISE'         = mise_table,
              'MISE_mu'      = mise_mu_table)
saveRDS(result,
        paste('results/experiment4_',
              paste(true_shifts_mean, collapse='-'), '_',
              paste(size_mean,        collapse='-'), '_',
              paste(true_shifts_var,  collapse='-'), '_',
              paste(size_var,         collapse='-'), '_',
              alpha_hat, '_', ctime, '.rds', sep=""))
