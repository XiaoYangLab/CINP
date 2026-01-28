# Generalization evaluation of CINP
# Author: Yang Xiao, PKU, 2025
# @: xiaoyang9604@gmail.com

rm(list=ls())

library(readxl)
library(glmnet)
library(stringr)
library(data.table)
library(caret)
library(doParallel) 
library(foreach)

# set parameters
mypath <- ' '
weight_idx <- ' ' 
test_idx <- " "
alpha_value <- seq(0, 1,by = 0.1)
rep_time <- 30
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

#load causal factors
causal_factor <- fread(paste0(weight_idx, '_penalty_factors.txt'))
wf <- as.matrix(causal_factor)
wf <- rbind(0, 0, wf) # age and sex occupation

# set training and testing sets
if (test_idx == "test") {
  xtrain <- as.matrix(cbind(struc_train, fiber_train, func_train))
  xtrain <- cbind(panss_train$sex, panss_train$age, xtrain)
  
  ytrain <- panss_train$allscores

  xtest <- as.matrix(cbind(struc_test, fiber_test, func_test))
  xtest <- cbind(panss_test$sex, panss_test$age, xtest)
  
  ytest <- panss_test$allscores
  
} else if (test_idx == "train") {
  xtrain <- as.matrix(cbind(struc_test, fiber_test, func_test))
  xtrain <- cbind(panss_test$sex, panss_test$age, xtrain)
  
  ytrain <- panss_test$allscores

  xtest <- as.matrix(cbind(struc_train, fiber_train, func_train))
  xtest <- cbind(panss_train$sex, panss_train$age, xtest)
  
  ytest <- panss_train$allscores
}

# evaluate generalizability using leave-one-site cross-validation  
cv_results <- vector("list", length(alpha_value))
for (i in 1:length(alpha_value)) cv_results[[i]] <- vector("list", rep_time)

min_cv_errors <- matrix(NA
                        , nrow = length(alpha_value)
                        , ncol = rep_time) 

xtrain_mean <- apply(xtrain, 2, mean)
xtrain_sd <- apply(xtrain, 2, sd)
xtrain <- scale(xtrain, center = xtrain_mean, scale = xtrain_sd)
xtest <- scale(xtest, center = xtrain_mean, scale = xtrain_sd)

set.seed(123)
for (rep in 1:rep_time) {
  
  for (i in seq_along(alpha_value)) {
    cv_fit <- cv.glmnet(xtrain
                        , ytrain
                        , alpha = alpha_value[i]
                        , family = "gaussian"
                        , type.measure = "mse"
                        , nfolds = 10
                        , penalty.factor = wf
                        , standardise = FALSE
                        , keep = TRUE
    )
    
    cv_results[[i]][[rep]] <- cv_fit
    min_cv_errors[i, rep] <- min(cv_fit$cvm)
    
  }
}

best_idx <- which(min_cv_errors == min(min_cv_errors), arr.ind = TRUE)

model <- glmnet(xtrain
                , ytrain
                , alpha = alpha_value[best_idx[1,1]]
                , lambda = cv_results[[best_idx[1,1]]][[best_idx[1,2]]]$lambda.min
                , family = "gaussian"
                , penalty.factor = wf
                , standardise = FALSE
                )

yhat <- predict(model
                , newx = xtest
                , s = cv_results[[best_idx[1,1]]][[best_idx[1,2]]]$lambda.min
                , family = "gaussian"
                , standardise = FALSE
                )

pred_cor <- cor.test(yhat, ytest, method = "pearson")
rmse <- RMSE(yhat, ytest)
plot(yhat, ytest)

# save results
write.table(yhat
            , paste0(test_idx, '_', weight_idx,'_yhat_all.txt')
            , col.names = FALSE, sep = " ", row.names = FALSE, quote = FALSE)

write.csv(as.matrix(coef(model))
          , paste0(test_idx, '_', weight_idx, '_coef_all_mr.csv'))

# permutations for estimating statistical significance
set.seed(123)
perm_y_matrix <- replicate(n_permutations, sample(ytrain))

cor_permuted <- foreach(perm = 1:n_permutations
                , .packages = c("glmnet", "caret", "stringr")
                , .combine = 'c'
                , .export = c("xtrain"
                              , "xtest"
                              , "ytest"
                              , "wf"
                              , "alpha_value"
                              , "rep_time"
                              , "perm_y_matrix")) %dopar% {
  
  ytrain_permuted <- perm_y_matrix[, perm]
  cv_results <- vector("list", length(alpha_value))
  for (i in 1:length(alpha_value)) cv_results[[i]] <- vector("list", rep_time)
  
  min_cv_errors <- matrix(NA
                          , nrow = length(alpha_value)
                          , ncol = rep_time) 
  for (rep in 1:rep_time) {
    
    for (i in seq_along(alpha_value)) {
      cv_fit <- cv.glmnet(xtrain
                          , ytrain_permuted
                          , alpha = alpha_value[i]
                          , family = "gaussian"
                          , type.measure = "mse"
                          , nfolds = 10
                          , penalty.factor = wf
                          , standardise = FALSE
                          , keep = TRUE
      )
      
      cv_results[[i]][[rep]] <- cv_fit
      min_cv_errors[i, rep] <- min(cv_fit$cvm)
      
    }
  }
  
  best_idx <- which(min_cv_errors == min(min_cv_errors), arr.ind = TRUE)
  
  model_permuted <- glmnet(xtrain
                           , ytrain_permuted
                           , alpha = alpha_value[best_idx[1,1]]
                           , lambda = cv_results[[best_idx[1,1]]][[best_idx[1,2]]]$lambda.min
                           , family = "gaussian"
                           , penalty.factor = wf
                           , standardise = FALSE
  )
  
  yhat_permuted <- predict(model_permuted
                           , newx = xtest
                           , s = cv_results[[best_idx[1,1]]][[best_idx[1,2]]]$lambda.min
                           , family = "gaussian"
                           , standardise = FALSE
  )

  return(cor(ytest, yhat_permuted, method = "pearson", use = "complete.obs"))
  
}
p_perm <- sum(cor_permuted >= pred_cor$estimate) / n_permutations

# permutations for penalty factors
set.seed(123)

cor_permuted <- foreach(perm = 1:n_permutations
                , .packages = c("glmnet", "caret", "stringr")
                , .combine = 'c'
                , .export = c("xtrain"
                              , "ytrain"
                              , "xtest"
                              , "ytest"
                              , "wf"
                              , "alpha_value"
                              , "rep_time")) %dopar% {
  
  wf_permuted <- as.numeric(wf[-c(1,2), 1])
  wf_permuted <- sample(wf_permuted)
  wf_permuted <- c(0, 0, wf_permuted)
  
  cv_results <- vector("list", length(alpha_value))
  for (i in 1:length(alpha_value)) cv_results[[i]] <- vector("list", rep_time)
  
  min_cv_errors <- matrix(NA
                          , nrow = length(alpha_value)
                          , ncol = rep_time) 
  for (rep in 1:rep_time) {
    
    for (i in seq_along(alpha_value)) {
      cv_fit <- cv.glmnet(xtrain
                          , ytrain
                          , alpha = alpha_value[i]
                          , family = "gaussian"
                          , type.measure = "mse"
                          , nfolds = 10
                          , penalty.factor = wf_permuted
                          , standardise = FALSE
                          , keep = TRUE
      )
      
      cv_results[[i]][[rep]] <- cv_fit
      min_cv_errors[i, rep] <- min(cv_fit$cvm)
      
    }
  }
  
  best_idx <- which(min_cv_errors == min(min_cv_errors), arr.ind = TRUE)
  
  model_permuted <- glmnet(xtrain
                           , ytrain
                           , alpha = alpha_value[best_idx[1,1]]
                           , lambda = cv_results[[best_idx[1,1]]][[best_idx[1,2]]]$lambda.min
                           , family = "gaussian"
                           , penalty.factor = wf_permuted
                           , standardise = FALSE
  )
  
  yhat_permuted <- predict(model_permuted
                           , newx = xtest
                           , s = cv_results[[best_idx[1,1]]][[best_idx[1,2]]]$lambda.min
                           , family = "gaussian"
                           , standardise = FALSE
  )
  
  return(cor(ytest, yhat_permuted, method = "pearson", use = "complete.obs"))
  
}
p_causality <- sum(cor_permuted >= pred_cor$estimate) / n_permutations

stopImplicitCluster()
cat("done!\n")

