# CINP
CINP: a causality-informed neuroimaging prediction

This repository provides the core implementation of the CINP (Causality-Informed Neuroimaging Prediction) framework, together with pre-trained models and auxiliary materials to facilitate transparent evaluation and replication.

CINP integrates Mendenlian randomization (MR)-derived causal priors into neuroimaging-based prediction by introducing feature-specific penalty factors within an elastic net regression framework.

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
The function below applies a trained CINP model to your neuroimaging data to predict treatment response. This wrapper is implemented in R and relies on the glmnet package:



