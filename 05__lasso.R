#!/usr/bin/env Rscript
#
# Wim (Feb-May 2025; w.m.otte@umcutrecht.nl)
#
# AIM: L2/L1 optimalization:  (Lasso), (Ridge), (Elastic)
#
################################################################################

# Load required package
library( 'readr' )
library( 'glmnet' )
library( 'survival' )
library( 'parameters' )
library( 'pec' )
library( 'riskRegression' )
library( 'timeROC' )
library( 'rms' ) # C-stat

################################################################################
# Custom functions
################################################################################



################################################################################
# End custom functions
################################################################################

# set seed
set.seed( 456 )

# output directory
outdir <- 'out.03.optimize'
dir.create( outdir, showWarnings = FALSE )

# load m imputed data.frames
imp <- readRDS( 'out.00.impute/mice_imputations.rds' )

# get first imputed data.frame
m <- 10
data <- mice::complete( imp, m )

all_preds <- c( "semiolmotor", "ageonset", "numaeds",
                "focal", "newgen", "oldgen", "female", "fh",
                "neonatalszs", "fs", "gt9szs", "priordc", "semiolgtc", 
                "benign", "dd", "abnormalexam", "abnormalimaging", "eegabnormal",
                "ageonsetcat", "eegepileptiform", "etiolunknown", "yearsszs", "yearsszfree",
                "dc", "etiolbirthtrauma", "brainsurgery", "semiolmyoclonic", "status", 
                "impairing", "ddd", "semioltc", "semiolabsence", "etioltbi",
                "etiolinfxn", "etiolstructural", "yearsszfreecat", "etiolstroke", "etioltumor",
                "etiolgenetic", "etiolmcd", "etiolstructuralother", "abnormalmri", "abnormalct",
                "etiolencephalomalac", "etiolperiventleuk", "etiolatrophy", "numszstc_cat" )

# TODO lambda optimization does no work properly (i.e., very low number)) and very little Deviance variation
sel <- read.csv( 'out.01.multivariable/significant_predictors.csv' )
sign_preds <- sel$Predictor


######################
# Lasso
######################

# Define survival outcome
y <- Surv( time = data$exit_time, event = data$event )

# Create a formula from the vector of predictor names
formula_str <- paste( "~", paste( sign_preds, collapse = " + " ) )
formula_obj <- as.formula(formula_str)

# Use the formula to create the design matrix
X <- model.matrix(formula_obj, data = data)[, -1]  # Remove intercept column

# Fit Cox model with cross-validation
cvfit_lasso <- cv.glmnet( X, y, family = "cox", alpha = 1 )  # alpha = 1 for Lasso

# Optimal lambda
plot( cvfit_lasso )
best_lambda_lasso <- cvfit_lasso$lambda.min
print( best_lambda_lasso )

# fit final model with optimal lambda
final_model_lasso <- glmnet( X, y, family = "cox", alpha = 1, lambda = best_lambda_lasso )


# Extract nonzero coefficients
round( coef( final_model_lasso ), 3 )





#####################################
# RIDGE
#####################################
# For L2 regularization (Ridge regression), you follow the same procedure as L1, 
# but set alpha = 0 in glmnet(). Ridge penalization shrinks coefficients towards 
# zero but does not set them exactly to zero, so all variables will be included in the model.

set.seed( 210 )  # Ensure reproducibility

# Fit Cox model with cross-validation
cvfit_ridge <- cv.glmnet( X, y, family = "cox", alpha = 0 )  # alpha = 0 for Ridge

# Optimal lambda
plot( cvfit_ridge )
best_lambda_ridge <- cvfit_ridge$lambda.min
print( best_lambda_ridge )

# fit final model with optimal lambda
final_model_ridge <- glmnet( X, y, family = "cox", alpha = 0, lambda = best_lambda_ridge )

# Extract nonzero coefficients
round( coef( final_model_ridge ), 3 )


# ELASTICS
##############################

set.seed( 105 )  # Ensure reproducibility

# Fit Elastic Cox model with cross-validation
cvfit_elastic <- cv.glmnet( X, y, family = "cox", alpha = 0.5 )  # alpha = 0.5 for Elastic

# Optimal lambda
plot( cvfit_elastic )
best_lambda_elastic <- cvfit_elastic$lambda.min
print( best_lambda_elastic )

# fit final model with optimal lambda
final_model_elastic <- glmnet( X, y, family = "cox", alpha = 0.5, lambda = best_lambda_elastic )

# Extract nonzero coefficients
round( coef( final_model_elastic ), 3 )


####### C-stat ###############



# Compute C-index for Lasso
lasso_pred <- predict( final_model_lasso, newx = X, type = "link" )
cindex_lasso <- rcorr.cens( -lasso_pred, y )[ 1 ]

# Compute C-index for Ridge
ridge_pred <- predict( final_model_ridge, newx = X, type = "link" )
cindex_ridge <- rcorr.cens( -ridge_pred, y)[ 1 ]

# Compute C-index for Elastic
elastic_pred <- predict( final_model_elastic, newx = X, type = "link" )
cindex_elastic <- rcorr.cens( -elastic_pred, y)[ 1 ]

# Print results
cat( "C-Index (Lasso):", round( cindex_lasso, 3 ), "\n" )
cat( "C-Index (Ridge):", round( cindex_ridge, 3 ), "\n" )
cat( "C-Index (Elastic):", round( cindex_elastic, 3 ), "\n" )

######################
######### AUC ########
######################

# Get predicted risk scores for Lasso and Ridge
lasso_risk <- predict( final_model_lasso, newx = X, type = "link" )
ridge_risk <- predict( final_model_ridge, newx = X, type = "link" )
elastic_risk <- predict( final_model_elastic, newx = X, type = "link" )


# Define time points for evaluation
time_points <- quantile(data$exit_time, probs = c(0.25, 0.5, 0.75))  # 25%, 50%, 75% percentiles

# Compute AUC for Lasso
auc_lasso <- timeROC(T = data$exit_time,  
                     delta = data$event,  
                     marker = -lasso_risk,  
                     cause = 1,  
                     times = time_points,  
                     iid = TRUE)

# Compute AUC for Ridge
auc_ridge <- timeROC(T = data$exit_time,  
                     delta = data$event,  
                     marker = -ridge_risk,  
                     cause = 1,  
                     times = time_points,  
                     iid = TRUE)

# Compute AUC for Elastic
auc_elastic <- timeROC(T = data$exit_time,  
                     delta = data$event,  
                     marker = -elastic_risk,  
                     cause = 1,  
                     times = time_points,  
                     iid = TRUE)
