# CINP
CINP: a causality-informed neuroimaging prediction

This repository provides the core implementation of the CINP (Causality-Informed Neuroimaging Prediction) framework, together with a trained CINP model and auxiliary materials to facilitate transparent evaluation and replication.

CINP integrates Mendelian randomization (MR)-derived causal priors into neuroimaging-based prediction by introducing feature-specific penalty factors within an elastic net regression framework.

## Workflow
<img width="1480" height="504" alt="image" src="https://github.com/user-attachments/assets/a5b93a4c-767f-4811-b532-80ad50ff78cc" />

## Contents
- Core code
  - Pipeline applied to train a CINP model to neuroimaging data.

- A trained CINP model
  - _CINP_whitematter.rds: a CINP model trained on white matter fiber tract metrics using custom data for predicting antipsychotic treatment response in patients with schizophrenia.
 
- CINP wrapper
   - Functions required to apply the trained CINP model to your neuroimaging data to predict antipsychotic treatment response.

- Demo datasets
  - Demo data and some feature extraction template to facilitate model application on your own datasets.

## Applying the CINP wrapper
The function below applies a trained CINP model to your neuroimaging data to predict treatment response. This wrapper is implemented in **R** and relies on the **glmnet** package:

```r
library(glmnet)

apply_CINP_model <- function(model_pkg, neuroimaging_feas, sex, age) {
  
  xnew <- as.matrix(neuroimaging_feas)
  xnew <- cbind(sex, age, xnew)
  
  xnew <- xnew[, model_pkg$feature_name, drop = FALSE]
  
  xnew <- scale(
    xnew,
    center = model_pkg$x_center,
    scale  = model_pkg$x_scale
  )
  
  yhat <- predict(
    model_pkg$model,
    newx = xnew,
    s    = model_pkg$lambda,
    family = "gaussian",
    standardise = FALSE
  )
  
  return(as.numeric(yhat))
}
```

## Usage
```r
model_pkg <- readRDS("_CINP_whitematter.rds")

y_pred <- apply_CINP_model(
    model_pkg,
    neuroimaging_feas = yourdata$neuroimaging_feas,
    sex = yourdata$sex,
    age = yourdata$age
)
```
## Resource
- A causal prior used to predicting antipsychotic treatment response was provied in our work.

## Contact
Please create an issue on the github repo if you encounter any problems. You can also contact the developers through email: xiaoyang9604@gmail.com.

## Citation




