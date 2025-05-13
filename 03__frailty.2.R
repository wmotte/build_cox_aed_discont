#!/usr/bin/env Rscript
#
# Wim Otte (May; w.m.otte@umcutrecht.nl)
#
# AIM: Mixed-effect Cox prop hazards model (with correction for center)
#
################################################################################

# Load required package
library('readr')
library('survival')
library('coxme')

################################################################################
# Custom functions
################################################################################

###
# Return best combination of predictors (mixed-effects) model
##
get_best_preds <- function(data)
{
    # Define all potential predictors
    all_preds <- c("semiolmotor", "ageonset", "numaeds",
                   "focal", "newgen", "oldgen", "female", "fh",
                   "neonatalszs", "fs", "gt9szs", "priordc", "semiolgtc", 
                   "benign", "dd", "abnormalexam", "abnormalimaging", "eegabnormal",
                   "ageonsetcat", "eegepileptiform", "etiolunknown", "yearsszs", "yearsszfree",
                   "dc", "etiolbirthtrauma", "brainsurgery", "semiolmyoclonic", "status", 
                   "impairing", "ddd", "semioltc", "semiolabsence", "etioltbi",
                   "etiolinfxn", "etiolstructural", "yearsszfreecat", "etiolstroke", "etioltumor",
                   "etiolgenetic", "etiolmcd", "etiolstructuralother", "abnormalmri", "abnormalct",
                   "etiolencephalomalac", "etiolperiventleuk", "etiolatrophy", "numszstc_cat")
    
    # Create the survival part of the formula
    surv_part <- "Surv(entry_time, exit_time, event)"
    
    # Start with base model (only random effect)
    base_formula_str <- paste(surv_part, "~ (1 | center)")
    base_formula <- as.formula(base_formula_str)
    
    # Initialize base model to calculate initial AIC
    base_model <- tryCatch({
        coxme(base_formula, data = data)
    }, error = function(e) {
        stop("Error fitting base model: ", conditionMessage(e))
    })
    
    # Calculate AIC for base model
    # Note: We use the negative log-likelihood *2 + 2*df as AIC
    base_aic <- -2 * base_model$loglik[2] + 2 * (length(base_model$coefficients) + 1)
    cat("Base model AIC:", base_aic, "\n")
    
    # Initialize variables to track best model
    best_model <- base_model
    best_aic <- base_aic
    successful_preds <- c()
    
    # Try each predictor one by one
    for (pred in all_preds) {
        print(pred)
        # Create formula with current successful predictors plus the new one being tested
        test_preds <- c(successful_preds, pred)
        test_formula_str <- paste(surv_part, "~", paste(test_preds, collapse = " + "), "+ (1 | center)")
        test_formula <- as.formula(test_formula_str)
        
        # Try fitting the model with the new predictor
        model_result <- tryCatch({
            model <- coxme(test_formula, data = data)
            current_aic <- -2 * model$loglik[2] + 2 * (length(model$coefficients) + 1)
            list(success = TRUE, model = model, aic = current_aic)
        }, error = function(e) {
            cat("Error adding predictor", pred, ":", conditionMessage(e), "\n")
            list(success = FALSE, model = NULL, aic = Inf)
        })
        
        # If successful and AIC improves, update the model
        if (model_result$success) {
            cat("Successfully added", pred, "- AIC:", model_result$aic, 
                "(Change: ", model_result$aic - best_aic, ")\n")
            
            if (model_result$aic < best_aic) {
                successful_preds <- test_preds
                best_model <- model_result$model
                best_aic <- model_result$aic
                cat("*** Keeping", pred, "- Better AIC ***\n")
            } else {
                cat("Not keeping", pred, "- AIC did not improve\n")
            }
        }
    }
    
    # Output results
    if (length(successful_preds) > 0) {
        cat("\nFinal model includes the following predictors:\n")
        cat(paste(successful_preds, collapse = ", "), "\n")
        cat("Final model AIC:", best_aic, "\n")
        cat("Improvement over base model:", best_aic - base_aic, "\n\n")
    } else {
        cat("No predictors improved the model. Using only random effect for center.\n")
    }
    
    return(successful_preds)
}

################################################################################
# End custom functions
################################################################################

# set seed
set.seed(456)

# output directory
outdir <- 'out.03.frailty'
dir.create(outdir, showWarnings = FALSE)

# load m imputed data.frames
imp <- readRDS('out.00.impute/mice_imputations.rds')

# Get number of imputations
num_imputations <- imp$m

# Define all potential predictors (same as in the function)
all_preds <- c("semiolmotor", "ageonset", "numaeds",
               "focal", "newgen", "oldgen", "female", "fh",
               "neonatalszs", "fs", "gt9szs", "priordc", "semiolgtc", 
               "benign", "dd", "abnormalexam", "abnormalimaging", "eegabnormal",
               "ageonsetcat", "eegepileptiform", "etiolunknown", "yearsszs", "yearsszfree",
               "dc", "etiolbirthtrauma", "brainsurgery", "semiolmyoclonic", "status", 
               "impairing", "ddd", "semioltc", "semiolabsence", "etioltbi",
               "etiolinfxn", "etiolstructural", "yearsszfreecat", "etiolstroke", "etioltumor",
               "etiolgenetic", "etiolmcd", "etiolstructuralother", "abnormalmri", "abnormalct",
               "etiolencephalomalac", "etiolperiventleuk", "etiolatrophy", "numszstc_cat")

# Create a matrix to store which predictors were selected in each imputation
selected_predictors <- matrix(FALSE, nrow=num_imputations, ncol=length(all_preds))
colnames(selected_predictors) <- all_preds

# Loop through all imputations
for (m in 1:num_imputations) {
    cat("\n\n=========================================\n")
    cat("Processing imputation", m, "of", num_imputations, "\n")
    cat("=========================================\n")
    
    # Get the complete data for this imputation
    data <- mice::complete(imp, m)
    
    # Get best predictors for this imputation
    best_preds <- get_best_preds(data)
    
    # Mark these predictors as selected in our matrix
    if (length(best_preds) > 0) {
        selected_predictors[m, best_preds] <- TRUE
    }
    
    cat("\nSelected predictors for imputation", m, ":\n")
    if (length(best_preds) > 0) {
        cat(paste(best_preds, collapse = ", "), "\n")
    } else {
        cat("None\n")
    }
}

# Calculate the frequency of selection for each predictor
selection_freq <- colSums(selected_predictors) / num_imputations

# Apply majority vote rule (select predictors chosen in > 50% of imputations)
final_predictors <- names(selection_freq[selection_freq > 0.5])

# Output results
cat("\n\n=========================================\n")
cat("FINAL PREDICTORS BASED ON MAJORITY VOTE:\n")
cat("=========================================\n")

if (length(final_predictors) > 0) {
    # Show all predictors with their selection frequencies
    cat("\nSelection frequencies for all predictors:\n")
    pred_freq_sorted <- sort(selection_freq, decreasing = TRUE)
    for (pred in names(pred_freq_sorted)) {
        freq_percent <- pred_freq_sorted[pred] * 100
        selected <- pred %in% final_predictors
        cat(sprintf("%s: %.1f%% %s\n", 
                    pred, freq_percent, 
                    ifelse(selected, "[SELECTED]", "")))
    }
    
    cat("\nFinal model includes the following predictors:\n")
    cat(paste(final_predictors, collapse = ", "), "\n")
} else {
    cat("No predictors were selected by majority vote.\n")
}

# Save the results
results <- list(
    selection_matrix = selected_predictors,
    selection_frequency = selection_freq,
    final_predictors = final_predictors
)
saveRDS(results, file = file.path(outdir, "majority_vote_predictors.rds"))

# Create a CSV with the selection frequencies
pred_freq_df <- data.frame(
    predictor = names(selection_freq),
    frequency = selection_freq,
    selected = names(selection_freq) %in% final_predictors
)

# Print summary message
cat("\nResults saved to", outdir, "\n")
cat("- Detailed results: majority_vote_predictors.rds\n")
cat("- Selection frequencies: predictor_selection_frequencies.csv\n")


# sort and write to disk
rownames( pred_freq_df ) <- NULL

# double sort
df_sorted <- pred_freq_df[order(pred_freq_df$selected, pred_freq_df$frequency, decreasing = c(TRUE, TRUE)), ]

# write to disk
write.csv( df_sorted, file = file.path(outdir, "predictor_selection_frequencies.csv"), row.names = FALSE)



