# Accuracy evaluation of CINP
# Author: Yang Xiao, PKU, 2025
# @: xiaoyang9604@gmail.com

rm(list=ls())

library(readxl)
library(glmnet)
library(stringr)
library(data.table)
library(caret)
library(dplyr)
library(tidyr)
library(doParallel) 
library(foreach)

# set parameters
mypath <- ' '
weight_idx <- ' '  
site_idx <- " "  
alpha_value <- seq(0, 1,by = 0.1)
rep_time_cv <- 50
rep_time_inner <- 30
n_permutations = 1000
n_cores <- parallel::detectCores() - 2 
registerDoParallel(n_cores)

# load datasets
panss_train <- read_excel(str_c(mypath,'demo_train_delta_panss.xlsx'))
panss_test <- read_excel(str_c(mypath,'demo_test_delta_panss.xlsx'))

#training set
struc_train <- read.csv(str_c(mypath, 'demo_struc_train_df.csv'), header = TRUE)
names<- struc_train[[1]]
struc_train <- struc_train[,-1]
rownames(struc_train) <- names

fiber_train <- read.csv(str_c(mypath, 'demo_fibers_train_df.csv'), header = TRUE)
fiber_train <- fiber_train[,-1]
rownames(fiber_train) <- names

func_train <- read.csv(str_c(mypath, 'demo_function_train_df.csv'), header = TRUE)
func_train <- func_train[,-1]
rownames(func_train) <- names

#testing set
struc_test <- read.csv(str_c(mypath, 'demo_struc_test_df.csv'), header = TRUE)
names<- struc_test[[1]]
struc_test <- struc_test[,-1]
rownames(struc_test) <- names

fiber_test <- read.csv(str_c(mypath, 'demo_fibers_test_df.csv'), header = TRUE)
fiber_test <- fiber_test[,-1]
rownames(fiber_test) <- names

func_test <- read.csv(str_c(mypath, 'demo_function_test_df.csv'),header = TRUE)
func_test <- func_test[,-1]
rownames(func_test) <- names

# load causal factors
causal_factor <- fread(paste0(weight_idx, '_penalty_factors.txt'))
wf <- as.matrix(causal_factor)
wf <- rbind(0, 0, wf) # age and sex occupation

# set datasets
if (site_idx == "train") {
  xdata <- as.matrix(cbind(struc_train, fiber_train, func_train))
  xdata <- cbind(panss_train$sex, panss_train$age, xdata)
  
  ydata <- panss_train$allscores
  
} else if (site_idx == "test") {
  xdata <- as.matrix(cbind(struc_test, fiber_test, func_test))
  xdata <- cbind(panss_test$sex, panss_test$age, xdata)
  
  ydata <- panss_test$allscores

}

# check size match between penalty factor and feature matrix
if (length(wf) != ncol(xdata)) {
  stop("Error: The penalty factor does not match the xdata!")
}

# cv function
run_single_cv <- function(xdata
                         , ydata
                         , wf
                         , alpha_value
                         , rep_time_inner
                         , folds) {
  
  yhat <- numeric(nrow(xdata))
  
  for (k in 1:length(unique(folds))) {
    test_idx <- which(folds == k)
    train_idx <- setdiff(1:nrow(xdata), test_idx)
    
    xtrain <- xdata[train_idx, ]
    ytrain <- ydata[train_idx]
    xtest <- xdata[test_idx, , drop = FALSE]
    
    # Standardize
    xtrain_mean <- colMeans(xtrain)
    xtrain_sd <- apply(xtrain, 2, sd)
    xtrain <- scale(xtrain, center = xtrain_mean, scale = xtrain_sd)
    xtest <- scale(xtest, center = xtrain_mean, scale = xtrain_sd)
    
    # Hyperparameter tuning
    grid_errors <- matrix(NA, nrow = length(alpha_value), ncol = rep_time_inner)
    grid_lambdas <- matrix(NA, nrow = length(alpha_value), ncol = rep_time_inner)
    
    for (i in seq_along(alpha_value)) {
      for (r in 1:rep_time_inner) {
        tryCatch({
          cv_fit <- cv.glmnet(xtrain
                              , ytrain
                              , alpha = alpha_value[i]
                              , family = "gaussian"
                              , type.measure = "mse"
                              , nfolds = 10
                              , penalty.factor = wf
                              , standardize = FALSE)
          grid_errors[i, r] <- min(cv_fit$cvm)
          grid_lambdas[i, r] <- cv_fit$lambda.min
        }, error = function(e) { 
          grid_errors[i, r] <- NA 
          grid_lambdas[i, r] <- NA
          })
      }
    }
    
    best_idx <- which(grid_errors == min(grid_errors), arr.ind = TRUE)
    best_alpha <- alpha_value[best_idx[1,1]]
    best_lambda <- grid_lambdas[best_idx[1,1], best_idx[1,2]]
    
    final_model <- glmnet(xtrain
                          , ytrain
                          , alpha = best_alpha
                          , lambda = best_lambda
                          , family = "gaussian"
                          , penalty.factor = wf
                          , standardize = FALSE)
    
    yhat[test_idx] <- predict(final_model
                              , newx = xtest
                              , s = best_lambda
                              , family = "gaussian"
                              , standardise = FALSE)
  }
  
  return(yhat)
}

# cross-validation setting (10-fold) with 50 replicates
set.seed(123)
folds_list <- replicate(rep_time_cv
                        ,sample(rep(1:10, length.out = nrow(xdata)))
                        ,simplify = FALSE
)

# training model
cat("Training with", rep_time_cv, "CV replicates on", n_cores, "cores...\n")

results_list <- foreach(rep_cv = 1:rep_time_cv
                        , .packages = c("glmnet", "caret", "stringr")
                        , .export = c("run_single_cv"
                                      , "xdata"
                                      , "ydata"
                                      , "wf"
                                      , "alpha_value"
                                      , "rep_time_inner"
                                      , "folds_list")) %dopar% {
  
  yhat <- run_single_cv(xdata
                        , ydata
                        , wf
                        , alpha_value
                        , rep_time_inner
                        , folds_list[[rep_cv]]
                        )
  
  r_val <- cor(yhat, ydata, method = "pearson", use = "complete.obs")
  rmse_val <- RMSE(yhat, ydata)
  
  # save mediate results
  write.table(yhat,
              file.path(mypath
                        , 'prediction'
                        , paste0(site_idx, '_' , weight_idx
                                 , '_yhat_rep', rep_cv, '_all.txt')),
              col.names = FALSE, sep = " ", row.names = FALSE, quote = FALSE)
  
  list(r = r_val, rmse = rmse_val)
}

results_df <- do.call(rbind, lapply(results_list, as.data.frame))
cat("\n================ Summary ================\n")
cat("Site:", site_idx, "| Weight:", weight_idx, "\n")
cat("Mean Pearson R:", mean(unlist(results_df$r)), " (SD:", sd(unlist(results_df$r)), ")\n")
cat("Mean RMSE:", mean(unlist(results_df$rmse)), "\n")

write.csv(results_df,
          file.path(mypath, 'prediction', paste0(site_idx, '_', weight_idx, '_performance_metrics.csv')))

# permutations for estimating statistical significance
set.seed(123)
perm_y_matrix <- replicate(n_permutations, sample(ydata))

cat("\nRunning", n_permutations, "permutations with", rep_time_cv, "CV reps each...\n")

cor_permuted <- foreach(perm = 1:n_permutations
                , .packages = c("glmnet", "caret", "stringr")
                , .combine = 'c'
                , .export = c("run_single_cv"
                              ,"xdata"
                              , "wf"
                              , "alpha_value"
                              , "rep_time_cv"
                              , "rep_time_inner"
                              , "folds_list"
                              , "perm_y_matrix")) %dopar% {
  
  ydata_permuted <- perm_y_matrix[, perm]
  r_reps <- numeric(rep_time_cv)
  
  for (rep_cv in 1:rep_time_cv) {
    
    yhat_perm <- run_single_cv(xdata
                               , ydata_permuted
                               , wf
                               , alpha_value
                               , rep_time_inner
                               , folds_list[[rep_cv]]
                               )
    
    r_reps[rep_cv] <- cor(yhat_perm, ydata_permuted, method = "pearson", use = "complete.obs")
  }
  return(mean(r_reps))                      
}
p_perm <- sum(cor_permuted >= mean(unlist(results_df$r))) / n_permutations

# permutations for penalty factors
set.seed(123)

cat("Running", n_permutations, "permutations for penalty factor with", rep_time_cv, "CV reps each...\n")

cor_permuted <- foreach(perm = 1:n_permutations
                , .packages = c("glmnet", "caret", "stringr")
                , .combine = 'c'
                , .export = c("run_single_cv"
                              ,"xdata"
                              , "wf"
                              , "alpha_value"
                              , "rep_time_cv"
                              , "rep_time_inner"
                              , "folds_list")) %dopar% {
  
  wf_permuted <- as.numeric(wf[-c(1,2), 1])
  wf_permuted <- sample(wf_permuted)
  wf_permuted <- c(0, 0, wf_permuted)
  
  r_reps <- numeric(rep_time_cv)
  
  for (rep_cv in 1:rep_time_cv) {
    
    yhat_perm <- run_single_cv(xdata
                               , ydata
                               , wf_permuted
                               , alpha_value
                               , rep_time_inner
                               , folds_list[[rep_cv]]
                               )
    
    r_reps[rep_cv] <- cor(yhat_perm, ydata, method = "pearson", use = "complete.obs")
  }
  return(mean(r_reps))                      
}

p_causality <- sum(cor_permuted >= mean(unlist(results_df$r))) / n_permutations

stopImplicitCluster()
cat("done!\n")


