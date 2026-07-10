# IBD-AI-Triage
Code for the paper "Interpretable Artificial Intelligence for Non-Invasive Triage of Inflammatory Bowel Disease: A Cross-Cohort Validation Study“
# Interpretable Artificial Intelligence for Non-Invasive Triage of Inflammatory Bowel Disease: A Cross-Cohort Validation Study

This repository contains the core analytical code and machine learning pipelines for the paper:  
**"Interpretable Artificial Intelligence for Non-Invasive Triage of Inflammatory Bowel Disease: A Cross-Cohort Validation Study"** (Currently under review at *The Lancet Digital Health*).

## 📌 Overview
This repository provides the R scripts used for:
1. **Shapley-AUC Feature Selection**: A novel framework combining Shapley marginal contributions and multi-method consistency (LASSO, Random Forest, XGBoost, SVM) to identify robust predictors.
2. **Reduced Model External Validation**: Strict direct validation using locked coefficients without local refitting, ensuring no data leakage during cross-cohort validation on the HMP2 dataset.

## 🛠 Prerequisites
All statistical analyses and machine learning modeling were performed using R software.
* **R version**: 4.5.2
* **Core Packages**: `xgboost`, `randomForest`, `e1071` (SVM), `glmnet` (LASSO), `shapr`, `pROC`, `dcurves`, `brglm2` (Firth's penalized logistic regression).

## 📂 File Structure
* `01_Shapley_AUC_Selection.R`: Script for calculating marginal AUC gains and cross-method consistency.
* `02_External_Validation_Reduced_Model.R`: Script for applying Z-score standardized internal coefficients to the external cohort.
* `dummy_data.csv`: A synthetic, randomly generated toy dataset to demonstrate how the code runs without compromising patient privacy. 

## 🔒 Data Availability
Due to institutional ethical restrictions and patient privacy protocols, the **de-identified individual participant data (IPD)** from the internal ZJU cohort cannot be uploaded publicly. Researchers seeking access to the real clinical data for independent replication should contact the corresponding author (Yue Cheng; 2026t006@zcmu.edu.cn) to sign a formal data access agreement. 

The external HMP2 validation dataset is publicly available at the IBDMDB project portal (https://ibdmdb.org).
