source('ShiVa.R')
library(phylolm)
library(l1ou)
source('method_helpers/l1ou_fixed.R')
library(PhylogeneticEM)
library(PCMBase)
library(PCMBaseCpp)
library(PCMFit)
library(data.table)

generatePCMModelsFunction <- function() {
# make results reproducible
set.seed(4, kind = "Mersenne-Twister", normal.kind = "Inversion")

PCMGenerateModelTypes()
# An example DefineParameterLimits.R file can be found in 
# vignettes/DefineParameterLimits.R of the PCMFit package source-code. 
# Note that the path here is relative to the working directory of
# the worker node R-process. 
source('DefineParameterLimits.R', local=FALSE)
}

######### Define simulation setting ######
file_path = ""  
alpha = 1
sigma2_0 = 2
m = 2   # Number of simulations per process
t = 25  # Number of independent processes per parameter setting

data('lizard.tree')
tree1 = lizard.tree
#tree1 = rcoal(5)
tree1$edge.length = tree1$edge.length/max(node.depth.edgelength(tree1))
X = generate_design_matrix(tree1,type='simpX')

args = commandArgs(T)
k = as.numeric(args[1])
runfile = read.csv('runfile.csv')


params = as.character(runfile[ceiling(k/t),])
true_shifts_mean = as.numeric(strsplit(params[2],',')[[1]])
size_mean = as.numeric(strsplit(params[3],',')[[1]])
true_shifts_var = as.numeric(strsplit(params[4],',')[[1]])
size_var = as.numeric(strsplit(params[5],',')[[1]])

V = OU.vcv(tree1,alpha)
tb = node.depth.edgelength(tree1)[tree1$edge[,1]]
q = exp(-2*alpha*max(node.depth.edgelength(tree1)))/(2*alpha)*(1-exp(2*alpha*tb))

gamma = rep(0,ncol(X))
gamma[true_shifts_var]  = size_var
Sigma = V*sigma2_0 +  X%*%diag(gamma)%*%t(X)*V+ X%*%diag(gamma*q)%*%t(X)

beta = rep(0,ncol(X))
beta[true_shifts_mean] = size_mean
test_data =get_test_data(tree1,Sigma,1000,alpha,beta)

r_phyloEM = data.frame(matrix(nrow=0,ncol=ncol(X)))
r_l1ou = list("pBIC"=data.frame(matrix(nrow=0,ncol=ncol(X))),
              "BIC"=data.frame(matrix(nrow=0,ncol=ncol(X))))
r_PCMFit = data.frame(matrix(nrow=0,ncol=ncol(X)))
r_ShiVa = list("shifts_mean"=data.frame(matrix(nrow=0,ncol=ncol(X))),
               "shifts_var"=data.frame(matrix(nrow=0,ncol=ncol(X))))

loglik_table = data.frame(matrix(nrow=0,ncol=6))
names(loglik_table) = c('ShiVa','l1ou+pBIC','l1ou+BIC','phyloEM','PCMFit','true')

computing_time_table = data.frame(matrix(nrow=0,ncol=5))
names(computing_time_table) = c('ShiVa','l1ou+pBIC','l1ou+BIC','phyloEM','PCMFit')


i = 0

while(i<m){


  tryCatch({
  eps = mvrnorm(n = 1, mu = rep(0,nrow(X)), Sigma = Sigma)
  Y = X%*%beta+eps

  adj_data = adjust_data(tree1,Y)
  tree2 = adj_data$tree
  Y2 = adj_data$Y

  #### ShiVa #####
  source('ShiVa_me.R')
  t = Sys.time()
  alpha_hat = phylolm(Y2~1, phy= tree2,model="OUfixedRoot")$optpar
  result1 = get_mean_var_shifts_model_selection(Y,tree1,alpha_hat,NULL,exp(1:10-6),top_k = 3, measurement_error = FALSE)
  ctime_row = c(Sys.time()-t)
  r_ShiVa$shifts_mean[nrow(r_ShiVa$shifts_mean)+1,] = as.numeric(result1$best_model$beta!=0)
  r_ShiVa$shifts_var[nrow(r_ShiVa$shifts_var)+1,] = as.numeric(result1$best_model$gamma!=0)  

  sv_mean = (1:ncol(X))[result1$best_model$beta!=0]
  sv_var = (1:ncol(X))[result1$best_model$gamma!=0]
  loglik_row = c(get_prediction_likelihood(result1$best_model$fitted.values,result1$best_model$Sigma,alpha,test_data))
  
###############l1ou#############

  t = Sys.time()
  l1ou_model1 = estimate_shift_configuration(tree2,Y2,max.nShifts=20,criterion="pBIC")
  
  ctime_row = c(ctime_row,Sys.time()-t)
  r_l1ou$pBIC[nrow(r_l1ou$pBIC)+1,] = as.numeric(1:ncol(X) %in% l1ou_model1$shift.configuration)
  lmod1 = if(length(l1ou_model1$shift.configuration)>0) phylolm(Y2~X[,l1ou_model1$shift.configuration],phy=tree2,model='OUfixedRoot') else phylolm(Y2~1,phy=tree2,model='OUfixedRoot')
  loglik_row = c(loglik_row, get_prediction_likelihood(lmod1$fitted.values,V*lmod1$sigma2,alpha,test_data))
  
  t = Sys.time()
  l1ou_model2 = estimate_shift_configuration(tree2,Y2,max.nShifts=20,criterion="BIC")

  ctime_row = c(ctime_row,Sys.time()-t)
  r_l1ou$BIC[nrow(r_l1ou$BIC)+1,] = as.numeric(1:ncol(X) %in% l1ou_model2$shift.configuration)
  lmod2 = if(length(l1ou_model2$shift.configuration)>0) phylolm(Y2~X[,l1ou_model2$shift.configuration],phy=tree2,model='OUfixedRoot') else phylolm(Y2~1,phy=tree2,model='OUfixedRoot')
  loglik_row = c(loglik_row, get_prediction_likelihood(lmod2$fitted.values,V*lmod2$sigma2,alpha,test_data))
 
 ##########phyloEM###########
  t = Sys.time()
  emmodel = PhyloEM(tree1,as.vector(Y),process='OU',K_max = 20)
  sv_em = params_process(emmodel)$shift$edges
  
  ctime_row = c(ctime_row,Sys.time()-t)
  r_phyloEM[nrow(r_phyloEM)+1,] = as.numeric(1:ncol(X) %in% sv_em)
  lmod3 = if(length(sv_em)>0) phylolm(Y2~X[,sv_em],phy=tree2,model='OUfixedRoot') else phylolm(Y2~1,phy=tree2,model='OUfixedRoot')
  loglik_row = c(loglik_row, get_prediction_likelihood(lmod3$fitted.values,V*lmod3$sigma2,alpha,test_data))

  ############PCMFit###########

  t = Sys.time()
  fitMGPM_A_F_BC2_RR <- PCMFitMixed(
      X = t(Y), tree = tree1, metaIFun = PCMInfoCpp,
      generatePCMModelsFun = generatePCMModelsFunction, 
      maxNumRoundRobins = 2, maxNumPartitionsInRoundRobins = 10,minCladeSizes = 2L,
      tableFitsPrev = NULL,
     # prefixFiles = prefixFiles,
      doParallel = FALSE)
 
  bestFit <- RetrieveBestFitScore(fitMGPM_A_F_BC2_RR)
  ctime_row = c(ctime_row,Sys.time()-t)
  r_PCMFit[nrow(r_PCMFit)+1,] = as.numeric(tree1$edge[,2] %in% as.numeric(bestFit$inferredRegimeNodes))
  loglik_row = c(loglik_row, get_prediction_likelihood(PCMBase::PCMMean(tree1,bestFit$inferredModel),PCMBase::PCMVar(tree1,bestFit$inferredModel),alpha,test_data))

  loglik_row = c(loglik_row, get_prediction_likelihood(X%*%beta,Sigma,alpha,test_data))

  loglik_table[nrow(loglik_table)+1,] = loglik_row
  computing_time_table[nrow(computing_time_table)+1,] = ctime_row

  i = i+1
  
  if(i == 1){ctime = proc.time()[1]*100}

  result = list('ShiVa' = r_ShiVa, 'l1ou' = r_l1ou,'phyloEM' = r_phyloEM, 'PCMFit'=r_PCMFit,'loglik'=loglik_table,'compute_time'=computing_time_table)
  saveRDS(result,paste(file_path,'/results/experiment_',paste(true_shifts_mean,collapse = '-'),'_',paste(size_mean,collapse = '-'),'_',paste(true_shifts_var,collapse = '-'),'_',paste(size_var,collapse = '-'),'_',ctime,'.rds',sep=""))


  }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


result = list('ShiVa' = r_ShiVa, 'l1ou' = r_l1ou,'phyloEM' = r_phyloEM, 'PCMFit'=r_PCMFit, 'loglik'=loglik_table,'compute_time'=computing_time_table)
saveRDS(result,paste(file_path,'/results/experiment_',paste(true_shifts_mean,collapse = '-'),'_',paste(size_mean,collapse = '-'),'_',paste(true_shifts_var,collapse = '-'),'_',paste(size_var,collapse = '-'),'_',ctime,'.rds',sep=""))
