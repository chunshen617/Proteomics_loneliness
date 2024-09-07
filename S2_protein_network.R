# Conducting protein-expression networks
# Chun Shen, 2024

# Load R packages
library(data.table)
library(netboost)
library(impute)

###Load data
#protein data and other covariates
data_prot <- fread('protein_UKB.csv',data.table = F)
#missing value
smiss <- apply(data_prot[,2:2921],1,function(x) length(which(is.na(x)==TRUE))/ncol(data_prot[,2:2921]))
#participants with missing protein data exceeding 50% were excluded
data_prot2 <- data_prot[smiss<0.5,pmiss<0.5]#sample N=46850

#permutation
prot_dat <- as.matrix(t(data_prot2[,2:2921]))
prot_impute <- impute.knn(prot_dat,k = 10, rowmax = 0.5, colmax = 0.5, maxp = 1500, rng.seed=66666)
dat_impute <- as.data.frame(t(prot_impute$data))

#netboost
results <- netboost(datan=dat_impute, stepno=20L,
                    filter_method = 'spearman',progress = 1000L,cores=10L,
                    soft_power=2, min_cluster_size=20L, n_pc=1, robust_PCs = TRUE,
                    scale=TRUE, ME_diss_thres=0.25, qc_plot=TRUE, method = 'spearman',verbose=3)

save(results,file = 'netboostResult.RData')

