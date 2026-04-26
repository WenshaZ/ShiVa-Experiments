library(ShiVa)
library(phylolm)
library(l1ou)
source('../method_helpers/l1ou_fixed.R')
library(PhylogeneticEM)
library(PCMBase)
library(PCMBaseCpp)
library(PCMFit)
library(data.table)
library(stringr)
library(MASS)
library(phytools)

generatePCMModelsFunction <- function() {
  set.seed(4, kind = "Mersenne-Twister", normal.kind = "Inversion")
  PCMGenerateModelTypes()
  source('DefineParameterLimits.R', local=FALSE)
}

######### Define simulation setting ######
file_path = ""
sigma2_0  = 2
m = 2     # Number of simulations per process
n_proc = 50  # Number of independent processes per parameter setting
tree1 = readRDS(paste0(file_path, "/real_data_tree.rds"))   # The tree of real data

args = commandArgs(T)
k = as.numeric(args[1])
model_list = c("ShiVa_model.rds",
                "l1ou_pBIC_model.rds",
                "l1ou_pBIC_alpha_model.rds",
                "l1ou_BIC_model.rds",
                "l1ou_BIC_alpha_model.rds",
                "phyloEM_model.rds",
                "phyloEM_alpha_model.rds",
                "PCMFit_model.rds",
                "PCMFit_alpha_model.rds")


model_file = model_list[ceiling(k / n_proc)]
model_name = str_sub(model_file, 1, -5)
sim_model = readRDS(paste0(file_path, "/", model_file))   # Model to base simulations on
X = generate_design_matrix(tree1, type = 'simpX')
Sigma = sim_model$Sigma
beta  = sim_model$beta
alpha = sim_model$alpha
V = ShiVa:::OU.vcv(tree1, alpha)

test_data = ShiVa:::get_test_data(tree1, Sigma, 1000, alpha, beta)

r_phyloEM = data.frame(matrix(nrow = 0, ncol = ncol(X)))
r_l1ou = list("pBIC" = data.frame(matrix(nrow = 0, ncol = ncol(X))),
              "BIC"  = data.frame(matrix(nrow = 0, ncol = ncol(X))))
r_PCMFit = data.frame(matrix(nrow = 0, ncol = ncol(X)))
r_ShiVa = list("shifts_mean" = data.frame(matrix(nrow = 0, ncol = ncol(X))),
               "shifts_var"  = data.frame(matrix(nrow = 0, ncol = ncol(X))))

loglik_table = data.frame(matrix(nrow = 0, ncol = 6))
names(loglik_table) = c('ShiVa', 'l1ou+pBIC', 'l1ou+BIC', 'phyloEM', 'PCMFit', 'true')

computing_time_table = data.frame(matrix(nrow = 0, ncol = 5))
names(computing_time_table) = c('ShiVa', 'l1ou+pBIC', 'l1ou+BIC', 'phyloEM', 'PCMFit')


i = 0

while (i < m) {
  if ((i - 1 < 0.1*m) & i >= 0.1*m) print('finish: 10%')
  if ((i - 1 < 0.2*m) & i >= 0.2*m) print('finish: 20%')
  if ((i - 1 < 0.3*m) & i >= 0.3*m) print('finish: 30%')
  if ((i - 1 < 0.4*m) & i >= 0.4*m) print('finish: 40%')
  if ((i - 1 < 0.5*m) & i >= 0.5*m) print('finish: 50%')
  if ((i - 1 < 0.6*m) & i >= 0.6*m) print('finish: 60%')
  if ((i - 1 < 0.7*m) & i >= 0.7*m) print('finish: 70%')
  if ((i - 1 < 0.8*m) & i >= 0.8*m) print('finish: 80%')
  if ((i - 1 < 0.9*m) & i >= 0.9*m) print('finish: 90%')

  tryCatch({
    eps = mvrnorm(n = 1, mu = rep(0, nrow(X)), Sigma = Sigma)
    Y = X %*% beta + eps

    adj_data = adjust_data(tree1, Y)
    tree2 = adj_data$tree
    Y2    = adj_data$Y

    #### ShiVa #####
    t0 = Sys.time()
    result1 = ShiVa(Y2, tree2, alpha = NULL, refit_alpha = TRUE,
                    lambda1_list = NULL, lambda2_list = exp(1:10*0.4-6),
                    nfolds = 8, top_k = 6, max.steps = 300,
                    lambda.type = "lambda.min", max.num.shifts = 16,
                    measurement_error = FALSE, verbose = FALSE)
    ctime_row = c(as.numeric(Sys.time() - t0, units = 'secs'))
    r_ShiVa$shifts_mean[nrow(r_ShiVa$shifts_mean) + 1, ] = as.numeric(result1$best_model$beta != 0)
    r_ShiVa$shifts_var[nrow(r_ShiVa$shifts_var) + 1, ]   = as.numeric(result1$best_model$gamma != 0)

    sv_mean = (1:ncol(X))[result1$best_model$beta != 0]
    sv_var  = (1:ncol(X))[result1$best_model$gamma != 0]
    loglik_row = c(ShiVa:::get_prediction_likelihood(result1$best_model$fitted.values,
                                                     result1$best_model$Sigma,
                                                     test_data))

    ############### l1ou #############
    t0 = Sys.time()
    l1ou_model1 = estimate_shift_configuration(tree2, Y2)
    ctime_row = c(ctime_row, as.numeric(Sys.time() - t0, units = 'secs'))
    r_l1ou$pBIC[nrow(r_l1ou$pBIC) + 1, ] = as.numeric(1:ncol(X) %in% l1ou_model1$shift.configuration)
    lmod1 = if (length(l1ou_model1$shift.configuration) > 0)
              phylolm(Y2 ~ X[, l1ou_model1$shift.configuration], phy = tree2, model = 'OUfixedRoot')
            else
              phylolm(Y2 ~ 1, phy = tree2, model = 'OUfixedRoot')
    loglik_row = c(loglik_row, ShiVa:::get_prediction_likelihood(lmod1$fitted.values,
                                                                  V * lmod1$sigma2,
                                                                  test_data))

    t0 = Sys.time()
    l1ou_model2 = estimate_shift_configuration(tree2, Y2, criterion = "BIC")
    ctime_row = c(ctime_row, as.numeric(Sys.time() - t0, units = 'secs'))
    r_l1ou$BIC[nrow(r_l1ou$BIC) + 1, ] = as.numeric(1:ncol(X) %in% l1ou_model2$shift.configuration)
    lmod2 = if (length(l1ou_model2$shift.configuration) > 0)
              phylolm(Y2 ~ X[, l1ou_model2$shift.configuration], phy = tree2, model = 'OUfixedRoot')
            else
              phylolm(Y2 ~ 1, phy = tree2, model = 'OUfixedRoot')
    loglik_row = c(loglik_row, ShiVa:::get_prediction_likelihood(lmod2$fitted.values,
                                                                  V * lmod2$sigma2,
                                                                  test_data))

    ########## phyloEM ###########
    t0 = Sys.time()
    emmodel = PhyloEM(force.ultrametric(tree2), as.vector(Y2), process = 'OU', K_max = 20)
    sv_em = params_process(emmodel)$shift$edges
    ctime_row = c(ctime_row, as.numeric(Sys.time() - t0, units = 'secs'))
    r_phyloEM[nrow(r_phyloEM) + 1, ] = as.numeric(1:ncol(X) %in% sv_em)
    lmod3 = if (length(sv_em) > 0)
              phylolm(Y2 ~ X[, sv_em], phy = tree2, model = 'OUfixedRoot')
            else
              phylolm(Y2 ~ 1, phy = tree2, model = 'OUfixedRoot')
    loglik_row = c(loglik_row, ShiVa:::get_prediction_likelihood(lmod3$fitted.values,
                                                                  V * lmod3$sigma2,
                                                                  test_data))

    ############ PCMFit ###########
    t0 = Sys.time()
    fitMGPM_A_F_BC2_RR <- PCMFitMixed(
        X = t(Y2), tree = tree2, metaIFun = PCMInfoCpp,
        generatePCMModelsFun = generatePCMModelsFunction,
        maxNumRoundRobins = 2, maxNumPartitionsInRoundRobins = 10, minCladeSizes = 2L,
        tableFitsPrev = NULL,
        doParallel = FALSE)
    bestFit <- RetrieveBestFitScore(fitMGPM_A_F_BC2_RR)
    ctime_row = c(ctime_row, as.numeric(Sys.time() - t0, units = 'secs'))
    r_PCMFit[nrow(r_PCMFit) + 1, ] = as.numeric(tree2$edge[, 2] %in% as.numeric(bestFit$inferredRegimeNodes))
    loglik_row = c(loglik_row, ShiVa:::get_prediction_likelihood(PCMBase::PCMMean(tree2, bestFit$inferredModel),
                                                                  PCMBase::PCMVar(tree2, bestFit$inferredModel),
                                                                  test_data))

    loglik_row = c(loglik_row, ShiVa:::get_prediction_likelihood(X %*% beta, Sigma, test_data))

    loglik_table[nrow(loglik_table) + 1, ] = loglik_row
    computing_time_table[nrow(computing_time_table) + 1, ] = ctime_row

    i = i + 1

    if (i == 1) ctime = proc.time()[1] * 100

    result = list('ShiVa' = r_ShiVa, 'l1ou' = r_l1ou, 'phyloEM' = r_phyloEM,
                  'PCMFit' = r_PCMFit, 'loglik' = loglik_table,
                  'compute_time' = computing_time_table)
    saveRDS(result, paste(file_path, '/results/results_', model_name, '_', ctime, '.rds', sep = ""))

  }, error = function(e) cat("ERROR :", conditionMessage(e), "\n"))
}


result = list('ShiVa' = r_ShiVa, 'l1ou' = r_l1ou, 'phyloEM' = r_phyloEM,
              'PCMFit' = r_PCMFit, 'loglik' = loglik_table,
              'compute_time' = computing_time_table)
saveRDS(result, paste(file_path, '/results/results_', model_name, '_', ctime, '.rds', sep = ""))
