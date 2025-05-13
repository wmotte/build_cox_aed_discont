#!/usr/bin/env Rscript
#
# Wim Otte (May 2025; w.m.otte@umcutrecht.nl)
#
# Build multivariable prediction model for epilepsy recurrence using
# multicenter data (IPD-meta-analysis).
#
# AIM: select variables for multivariable model, fit model and validate.
#
# It performs a comprehensive leave-one-center-out validation with the following features:
#    
# For each center:
#    
# Train on data from all other centers
# Test on the left-out center
# Handle multiple imputations correctly
#
#
# Performance metrics with 95% CIs:
#   
# C-index (Harrell's concordance)
# Sensitivity and specificity at optimal cutpoints
# Calibration slope and intercept
#
#
# Visualizations:
#
# Forest plots of C-index by center
# AUC plots for a time point with sufficient events
# Calibration slope plots
# 
#
# Summary tables:
#
# Center-specific performance metrics with 95% CIs
# Overall pooled performance with bootstrap 95% CIs
################################################################################

# Load required packages
library('readr')
library('parameters')
library('survival')
library('mice')
library('pROC')
library('rms')
library('timeROC')
library('gt') # tables
library('dplyr')
library('ggplot2')
library('boot')

################################################################################
# Start Custom Functions
################################################################################

###
# Function to perform leave-one-out center validation
##
leave_one_center_out_validation <- function(imp_data, formula_str, centers) 
{
    # Initialize results storage
    results <- list()
    
    center_val <- centers[ 1 ]
    
    # For each center
    for (center_val in centers) {
        cat("\nValidating on center:", center_val, "\n")
        
        # Setup results storage for this center
        results[[as.character(center_val)]] <- list(
            c_index = NULL,
            auc = NULL,
            sensitivity = NULL,
            specificity = NULL,
            calibration_slope = NULL,
            calibration_intercept = NULL
        )
        
        # Get formula object
        formula_obj <- as.formula(formula_str)
        
        # For each imputed dataset
        imp_models <- list()
        imp_preds <- list()
        imp_test_data <- list()
        
        m <- 1
        
        for (m in 1:imp_data$m) {
            # Extract complete data for this imputation
            complete_data <- complete(imp_data, m)
            
            # Split into training (all centers except current) and test (current center only)
            train_data <- complete_data[complete_data$center != center_val, ]
            test_data <- complete_data[complete_data$center == center_val, ]
            
            # Store test data for later
            imp_test_data[[m]] <- test_data
            
            # Fit model on training data
            model <- tryCatch({
                survival::coxph(formula = formula_obj, data = train_data)
            }, error = function(e) {
                message("Error in model fitting for imputation ", m, ", center ", center_val)
                message(e)
                return(NULL)
            })
            
            # Skip if model fitting failed
            if (is.null(model)) next
            
            # Store model
            imp_models[[m]] <- model
            
            # Get linear predictors for test data
            lp <- predict(model, newdata = test_data, type = "lp")
            
            # Store predictions
            imp_preds[[m]] <- data.frame(
                lp = lp,
                risk = exp(lp),
                event = test_data$event,
                entry_time = test_data$entry_time,
                exit_time = test_data$exit_time,
                survival_time = test_data$exit_time - test_data$entry_time
            )
        }
        
        # Skip if all models failed
        if (length(imp_models) == 0) {
            warning("All models failed for center ", center_val)
            next
        }
        
        # Pool results across imputations
        
        # 1. C-index (Harrell's concordance)
        # Initialize vector to store c-indices
        c_indices <- numeric(length(imp_models))
        
        # Loop through each imputation
        for (p in seq_along(imp_models)) {
            # Skip if model is NULL
            if (is.null(imp_models[[p]])) {
                message("Skipping c-index calculation for imputation ", p, " due to NULL model")
                c_indices[p] <- NA
                next
            }
            
            # Initialize c-index
            c_index <- 0
            
            # Get the model and test data for this imputation
            model <- imp_models[[p]]
            
            # Make sure we have test data for this imputation
            if (is.null(imp_test_data[[p]])) {
                message("Missing test data for imputation ", p)
                c_indices[p] <- NA
                next
            }
            
            val_data <- imp_test_data[[p]]
            
            # Try to calculate c-index
            tryCatch({
                # Print debugging info
                message("Calculating c-index for imputation ", p, ", test data rows: ", nrow(val_data))
                
                # Get linear predictor for test data
                lp <- predict(model, newdata = val_data, type = "lp")
                
                # Check if lp calculation worked
                if (length(lp) != nrow(val_data)) {
                    message("LP length mismatch: ", length(lp), " vs ", nrow(val_data))
                    c_indices[p] <- NA
                    next
                }
                
                # Add LP to data for validation
                val_data$lp <- lp
                
                # Fit new model with only LP as predictor
                new_model <- coxph(Surv(time = entry_time, time2 = exit_time, event) ~ lp, 
                                   data = val_data)
                
                # Extract c-index
                c_index <- summary(new_model)$concordance[1]
                message("Imputation ", p, " c-index: ", c_index)
                
                # Store result
                c_indices[p] <- c_index
            }, error = function(e) {
                message("Error calculating c-index for imputation ", p, ": ", conditionMessage(e))
                c_indices[p] <- NA
            })
        }
        
        # Print summary of c-indices for debugging
        message("C-indices for center ", center_val, ":")
        for (m in seq_along(c_indices)) {
            message("Imputation ", m, ": ", ifelse(is.na(c_indices[m]), "NA", round(c_indices[m], 3)))
        }
        
        # Calculate mean and standard error for C-index
        c_index_mean <- mean(c_indices, na.rm = TRUE)
        c_index_se <- sd(c_indices, na.rm = TRUE) / sqrt(sum(!is.na(c_indices)))
        c_index_lower <- c_index_mean - 1.96 * c_index_se
        c_index_upper <- c_index_mean + 1.96 * c_index_se
        
        # Print final results for this center
        message("Center ", center_val, " C-index: ", round(c_index_mean, 3), 
                " (95% CI: ", round(c_index_lower, 3), "-", round(c_index_upper, 3), ")")
        
        
        # Calculate mean and standard error for C-index
        c_index_mean <- mean(c_indices, na.rm = TRUE)
        c_index_se <- sd(c_indices, na.rm = TRUE) / sqrt(sum(!is.na(c_indices)))
        c_index_lower <- c_index_mean - 1.96 * c_index_se
        c_index_upper <- c_index_mean + 1.96 * c_index_se
        
        # Store C-index
        results[[as.character(center_val)]]$c_index <- c(c_index_mean, c_index_lower, c_index_upper)
        
        # 2. AUC at specific time points
        # Need to define time points based on your data
        time_points <- c( 2 ) # years, based on 'exit_time'
        
        # Convert to days if your survival time is in months
        time_points_days <- time_points * 12
   
             
        # Calculate AUC for each time point
        auc_results <- list()
        
        t <- time_points_days[ 1 ]
        
        for (t in time_points_days) {
            auc_values <- numeric(length(imp_preds))
            
            for (m in 1:length(imp_preds)) {
                if (is.null(imp_preds[[m]])) next
                
                tryCatch({
                    # Use timeROC package for time-dependent AUC
                    roc_obj <- timeROC::timeROC(
                        T = imp_preds[[m]]$survival_time,
                        delta = imp_preds[[m]]$event,
                        marker = imp_preds[[m]]$lp,
                        cause = 1,
                        times = t,
                        ROC = TRUE
                    )
                    
                    auc_values[m] <- roc_obj$AUC[2]
                }, error = function(e) {
                    message("Error calculating AUC for imputation ", m, " at time ", t)
                    auc_values[m] <- NA
                })
            }
            
            # Calculate mean and CI for AUC
            auc_mean <- mean(auc_values, na.rm = TRUE)
            auc_se <- sd(auc_values, na.rm = TRUE) / sqrt(sum(!is.na(auc_values)))
            auc_lower <- auc_mean - 1.96 * auc_se
            auc_upper <- auc_mean + 1.96 * auc_se
            
            auc_results[[paste0(time_points[which(time_points_days == t)], "yr")]] <- 
                c(auc_mean, auc_lower, auc_upper)
        }
        
        # Store AUC results
        results[[as.character(center_val)]]$auc <- auc_results
        
        # 3. Sensitivity and Specificity at optimal cutoff for each time point
        sens_spec_results <- list()
        
        for (t_idx in 1:length(time_points)) {
            t <- time_points_days[t_idx]
            t_label <- paste0(time_points[t_idx], "yr")
            
            sens_values <- numeric(length(imp_preds))
            spec_values <- numeric(length(imp_preds))
            
            for (m in 1:length(imp_preds)) {
                if (is.null(imp_preds[[m]])) next
                
                tryCatch({
                    # Calculate status at time t
                    status_t <- rep(0, nrow(imp_preds[[m]]))
                    eligible <- imp_preds[[m]]$survival_time <= t & imp_preds[[m]]$event == 1
                    status_t[eligible] <- 1
                    
                    # Create ROC curve
                    roc_obj <- pROC::roc(status_t, imp_preds[[m]]$lp, quiet = TRUE)
                    
                    # Find optimal cutpoint (Youden's J statistic)
                    opt_cut <- pROC::coords(roc_obj, "best", best.method = "youden")
                    
                    sens_values[m] <- opt_cut$sensitivity
                    spec_values[m] <- opt_cut$specificity
                }, error = function(e) {
                    message("Error calculating sensitivity/specificity for imputation ", m, " at time ", t)
                    sens_values[m] <- NA
                    spec_values[m] <- NA
                })
            }
            
            # Calculate mean and CI for sensitivity
            sens_mean <- mean(sens_values, na.rm = TRUE)
            sens_se <- sd(sens_values, na.rm = TRUE) / sqrt(sum(!is.na(sens_values)))
            sens_lower <- sens_mean - 1.96 * sens_se
            sens_upper <- sens_mean + 1.96 * sens_se
            
            # Calculate mean and CI for specificity
            spec_mean <- mean(spec_values, na.rm = TRUE)
            spec_se <- sd(spec_values, na.rm = TRUE) / sqrt(sum(!is.na(spec_values)))
            spec_lower <- spec_mean - 1.96 * spec_se
            spec_upper <- spec_mean + 1.96 * spec_se
            
            sens_spec_results[[t_label]] <- list(
                sensitivity = c(sens_mean, sens_lower, sens_upper),
                specificity = c(spec_mean, spec_lower, spec_upper)
            )
        }
        
        # Store sensitivity/specificity results
        results[[as.character(center_val)]]$sensitivity <- lapply(sens_spec_results, function(x) x$sensitivity)
        results[[as.character(center_val)]]$specificity <- lapply(sens_spec_results, function(x) x$specificity)
        
        # 4. Calibration (slope and intercept)
        # This is more complex for survival data with multiple imputations
        # We use approach from val.surv in rms package
        
        cal_slopes <- numeric(length(imp_models))
        cal_intercepts <- numeric(length(imp_models))
        
        for (m in 1:length(imp_models)) {
            if (is.null(imp_models[[m]]) || is.null(imp_test_data[[m]])) next
            
            tryCatch({
                # Get predictions
                lp <- predict(imp_models[[m]], newdata = imp_test_data[[m]], type = "lp")
                
                # Fit calibration model
                cal_model <- coxph(
                    Surv(time = entry_time, time2 = exit_time, event) ~ lp, 
                    data = imp_test_data[[m]]
                )
                
                # Extract slope (coefficient of lp)
                cal_slopes[m] <- coef(cal_model)[1]
                
                # For intercept in Cox model, need more complex approach
                # Here we use a simplified approach based on comparing baseline hazards
                baseline_orig <- basehaz(imp_models[[m]])
                baseline_val <- basehaz(cal_model)
                
                # Compare at median time
                median_time <- median(imp_test_data[[m]]$exit_time - imp_test_data[[m]]$entry_time)
                
                # Find closest time points
                orig_idx <- which.min(abs(baseline_orig$time - median_time))
                val_idx <- which.min(abs(baseline_val$time - median_time))
                
                # Calculate log ratio of baseline hazards as intercept
                cal_intercepts[m] <- log(baseline_val$hazard[val_idx] / baseline_orig$hazard[orig_idx])
            }, error = function(e) {
                message("Error calculating calibration for imputation ", m)
                cal_slopes[m] <- NA
                cal_intercepts[m] <- NA
            })
        }
        
        # Calculate mean and CI for calibration slope
        slope_mean <- mean(cal_slopes, na.rm = TRUE)
        slope_se <- sd(cal_slopes, na.rm = TRUE) / sqrt(sum(!is.na(cal_slopes)))
        slope_lower <- slope_mean - 1.96 * slope_se
        slope_upper <- slope_mean + 1.96 * slope_se
        
        # Calculate mean and CI for calibration intercept
        intercept_mean <- mean(cal_intercepts, na.rm = TRUE)
        intercept_se <- sd(cal_intercepts, na.rm = TRUE) / sqrt(sum(!is.na(cal_intercepts)))
        intercept_lower <- intercept_mean - 1.96 * intercept_se
        intercept_upper <- intercept_mean + 1.96 * intercept_se
        
        # Store calibration results
        results[[as.character(center_val)]]$calibration_slope <- c(slope_mean, slope_lower, slope_upper)
        results[[as.character(center_val)]]$calibration_intercept <- c(intercept_mean, intercept_lower, intercept_upper)
    }
    
    return(results)
}

###
# Function to create plots for visualizing validation results
##
plot_validation_results <- function(validation_results) {
    # Extract center names
    centers <- names(validation_results)
    
    # 1. Plot C-index
    c_index_data <- data.frame(
        center = centers,
        c_index = sapply(validation_results, function(x) x$c_index[1]),
        lower = sapply(validation_results, function(x) x$c_index[2]),
        upper = sapply(validation_results, function(x) x$c_index[3])
    )
    
    c_index_plot <- ggplot(c_index_data, aes(x = center, y = c_index)) +
        geom_point(size = 3) +
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
        geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
        labs(title = "C-index by Center", x = "Center", y = "C-index") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # 2. Plot AUC for each time point
    # Example for 2-year AUC
    auc_2yr_data <- data.frame(
        center = centers,
        auc = sapply(validation_results, function(x) x$auc[["2yr"]][1]),
        lower = sapply(validation_results, function(x) x$auc[["2yr"]][2]),
        upper = sapply(validation_results, function(x) x$auc[["2yr"]][3])
    )
    
    auc_2yr_plot <- ggplot(auc_2yr_data, aes(x = center, y = auc)) +
        geom_point(size = 3) +
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
        geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
        labs(title = "2-Year AUC by Center", x = "Center", y = "AUC") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # 3. Plot calibration slope
    cal_slope_data <- data.frame(
        center = centers,
        slope = sapply(validation_results, function(x) x$calibration_slope[1]),
        lower = sapply(validation_results, function(x) x$calibration_slope[2]),
        upper = sapply(validation_results, function(x) x$calibration_slope[3])
    )
    
    cal_slope_plot <- ggplot(cal_slope_data, aes(x = center, y = slope)) +
        geom_point(size = 3) +
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
        labs(title = "Calibration Slope by Center", x = "Center", y = "Slope") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Return the plots
    return(list(
        c_index_plot = c_index_plot,
        auc_2yr_plot = auc_2yr_plot,
        cal_slope_plot = cal_slope_plot
    ))
}

###
# Function to create a summary table of validation results
##
create_summary_table <- function(validation_results) {
    # Extract center names
    centers <- names(validation_results)
    
    # Create data frame for summary
    summary_df <- data.frame(
        Center = centers,
        
        # C-index
        C_Index = sapply(validation_results, function(x) 
            sprintf("%.3f (%.3f-%.3f)", x$c_index[1], x$c_index[2], x$c_index[3])),
        
        # 2-year AUC
        AUC_2yr = sapply(validation_results, function(x) 
            sprintf("%.3f (%.3f-%.3f)", x$auc[["2yr"]][1], x$auc[["2yr"]][2], x$auc[["2yr"]][3])),
        
        
        # 2-year Sensitivity
        Sens_2yr = sapply(validation_results, function(x) 
            sprintf("%.3f (%.3f-%.3f)", x$sensitivity[["2yr"]][1], x$sensitivity[["2yr"]][2], x$sensitivity[["2yr"]][3])),
        
        # 2-year Specificity
        Spec_2yr = sapply(validation_results, function(x) 
            sprintf("%.3f (%.3f-%.3f)", x$specificity[["2yr"]][1], x$specificity[["2yr"]][2], x$specificity[["2yr"]][3])),
        
        # Calibration slope
        Cal_Slope = sapply(validation_results, function(x) 
            sprintf("%.3f (%.3f-%.3f)", x$calibration_slope[1], x$calibration_slope[2], x$calibration_slope[3])),
        
        # Calibration intercept
        Cal_Intercept = sapply(validation_results, function(x) 
            sprintf("%.3f (%.3f-%.3f)", x$calibration_intercept[1], x$calibration_intercept[2], x$calibration_intercept[3]))
    )
    
    return(summary_df)
}

### 
# Backward Elimination using AIC
###
backward_elimination_mi_aic <- function(imp, initial_predictors, trace = FALSE) {
    cat("\nStarting backward elimination using AIC criterion\n")
    
    if(length(initial_predictors) == 0) {
        cat("No predictors provided. Stopping.\n")
        return(NULL)
    }
    
    cat("Initial predictors:", paste(initial_predictors, collapse = ", "), "\n")
    
    # Track iterations
    iteration <- 1
    removed_vars <- c()
    current_predictors <- initial_predictors
    continue <- TRUE
    
    # Store current best AIC
    current_formula_str <- paste0("Surv(time = entry_time, time2 = exit_time, event) ~ ", 
                                  paste(current_predictors, collapse = " + "))
    
    # Fit current full model
    current_fit <- with(imp, {
        formula_obj <- as.formula(current_formula_str)
        survival::coxph(formula = formula_obj)
    })
    
    # Calculate pooled AIC (average across imputations)
    calculate_pooled_aic <- function(mira_obj) {
        # Extract AIC from each imputed model fit
        aics <- sapply(mira_obj$analyses, AIC)
        # Return mean AIC
        mean(aics)
    }
    
    current_aic <- calculate_pooled_aic(current_fit)
    cat("Starting AIC:", current_aic, "\n")
    
    while(continue && length(current_predictors) > 1) {  # Need at least 1 predictor
        cat("\n*** Iteration", iteration, "***\n")
        cat("Current predictors:", paste(current_predictors, collapse = ", "), "\n")
        
        # Track best model from this iteration
        best_aic <- current_aic
        best_formula <- current_formula_str
        best_fit <- current_fit
        predictor_to_remove <- NULL
        
        # Try removing each predictor and check AIC
        for(i in 1:length(current_predictors)) {
            test_predictors <- current_predictors[-i]
            test_formula_str <- paste0("Surv(time = entry_time, time2 = exit_time, event) ~ ", 
                                       paste(test_predictors, collapse = " + "))
            
            if(trace) cat("Testing formula:", test_formula_str, "\n")
            
            # Fit model without this predictor
            test_fit <- tryCatch({
                with(imp, {
                    formula_obj <- as.formula(test_formula_str)
                    survival::coxph(formula = formula_obj)
                })
            }, error = function(e) {
                if(trace) cat("Error fitting model:", conditionMessage(e), "\n")
                return(NULL)
            })
            
            if(!is.null(test_fit)) {
                # Calculate AIC for this reduced model
                test_aic <- calculate_pooled_aic(test_fit)
                
                if(trace) {
                    cat("  Removed predictor:", current_predictors[i], "\n")
                    cat("  AIC:", test_aic, "\n")
                }
                
                # If this is better than our current best, update
                if(test_aic < best_aic) {
                    best_aic <- test_aic
                    best_formula <- test_formula_str
                    best_fit <- test_fit
                    predictor_to_remove <- current_predictors[i]
                }
            }
        }
        
        # Check if we found a better model
        if(!is.null(predictor_to_remove) && best_aic < current_aic) {
            cat("Removing predictor", predictor_to_remove, "improves AIC from", 
                current_aic, "to", best_aic, "\n")
            
            # Update current model
            current_aic <- best_aic
            current_formula_str <- best_formula
            current_fit <- best_fit
            
            # Update predictors
            removed_vars <- c(removed_vars, predictor_to_remove)
            current_predictors <- current_predictors[current_predictors != predictor_to_remove]
        } else {
            cat("No improvement in AIC by removing any predictor. Stopping.\n")
            continue <- FALSE
        }
        
        iteration <- iteration + 1
        
        # Safety check
        if(iteration > 100) {
            cat("Warning: Maximum iterations reached. Stopping.\n")
            break
        }
    }
    
    # Final model
    if(length(current_predictors) > 0) {
        cat("\nFinal model formula after backward elimination (AIC):\n")
        cat(current_formula_str, "\n")
        cat("Final AIC:", current_aic, "\n")
        
        cat("\nFinal model coefficients after backward elimination:\n")
        final_params <- parameters::parameters(current_fit)
        print(final_params)
        
        return(list(
            formula = current_formula_str,
            fit = current_fit,
            parameters = final_params,
            aic = current_aic,
            removed = removed_vars,
            remaining = current_predictors
        ))
    } else {
        cat("No predictors remain in the model after elimination.\n")
        return(NULL)
    }
}

###
# Calculate overall performance by averaging across centers
##
calculate_overall_performance <- function(validation_results)
{
    # Initialize storage
    overall <- list(
        c_index = NULL,
        auc_2yr = NULL,
        
        sensitivity_2yr = NULL,
        specificity_2yr = NULL,
        calibration_slope = NULL,
        calibration_intercept = NULL
    )
    
    # Extract values for each metric
    c_index_values <- sapply(validation_results, function(x) x$c_index[1])
    auc_2yr_values <- sapply(validation_results, function(x) x$auc[["2yr"]][1])
    sens_2yr_values <- sapply(validation_results, function(x) x$sensitivity[["2yr"]][1])
    spec_2yr_values <- sapply(validation_results, function(x) x$specificity[["2yr"]][1])
    cal_slope_values <- sapply(validation_results, function(x) x$calibration_slope[1])
    cal_int_values <- sapply(validation_results, function(x) x$calibration_intercept[1])
    
    # Calculate means and bootstrap 95% CIs
    bootstrap_ci <- function(values, R = 1000) {
        boot_means <- replicate(R, mean(sample(values, replace = TRUE)))
        c(mean(values), quantile(boot_means, 0.025), quantile(boot_means, 0.975))
    }
    
    overall$c_index <- bootstrap_ci(c_index_values)
    overall$auc_2yr <- bootstrap_ci(na.omit( auc_2yr_values ))
    
    overall$sensitivity_2yr <- bootstrap_ci(sens_2yr_values)
    overall$specificity_2yr <- bootstrap_ci(spec_2yr_values)
    overall$calibration_slope <- bootstrap_ci(cal_slope_values)
    overall$calibration_intercept <- bootstrap_ci(cal_int_values)
    
    # Create summary table
    overall_summary <- data.frame(
        Metric = c("C-index", "2-Year AUC", 
                   "2-Year Sensitivity", "2-Year Specificity", 
                   "Calibration Slope", "Calibration Intercept"),
        
        Value = c(
            sprintf("%.3f (%.3f-%.3f)", overall$c_index[1], overall$c_index[2], overall$c_index[3]),
            sprintf("%.3f (%.3f-%.3f)", overall$auc_2yr[1], overall$auc_2yr[2], overall$auc_2yr[3]),
            
            sprintf("%.3f (%.3f-%.3f)", overall$sensitivity_2yr[1], overall$sensitivity_2yr[2], overall$sensitivity_2yr[3]),
            sprintf("%.3f (%.3f-%.3f)", overall$specificity_2yr[1], overall$specificity_2yr[2], overall$specificity_2yr[3]),
            sprintf("%.3f (%.3f-%.3f)", overall$calibration_slope[1], overall$calibration_slope[2], overall$calibration_slope[3]),
            sprintf("%.3f (%.3f-%.3f)", overall$calibration_intercept[1], overall$calibration_intercept[2], overall$calibration_intercept[3])
        )
    )
    
    return(overall_summary)
}

###
# Save parameters to RTF
##
save_parameters_to_rtf <- function( final_models )
{
    # Plot to html: Using gt package (for publication-quality tables)
    parameters::parameters(final_models) %>%
        gt() %>%
        tab_header(
            title = "Fixed Effects Model Parameters"
        ) %>%
        # Format numeric columns with appropriate decimal places
        fmt_number(
            columns = c("Coefficient", "SE"),
            decimals = 3
        ) %>%
        fmt_number(
            columns = c("CI_low", "CI_high"),  # Format confidence interval columns
            decimals = 2
        ) %>%
        fmt_number(
            columns = "Statistic",  # Format test statistic
            decimals = 2
        ) %>%
        fmt_number(
            columns = "df_error",  # Format degrees of freedom
            decimals = 1
        ) %>%
        fmt_number(
            columns = "p",
            decimals = 3
        ) %>%
        # Highlight significant values
        tab_style(
            style = list(
                cell_fill(color = "#E8F4F8"),
                cell_text(weight = "bold")
            ),
            locations = cells_body(
                rows = p < 0.05
            )
        ) %>%
        # Rename columns for clarity
        cols_label(
            Parameter = "Predictor",
            Coefficient = "Coef.",
            SE = "Std. Error",
            CI_low = "CI Lower",
            CI_high = "CI Upper",
            df_error = "degr.freedom"
        ) %>%
        tab_footnote(
            footnote = "Values in bold indicate statistical significance (p < 0.05)",
            locations = cells_column_labels(
                columns = "p"
            )
        ) %>% 
        gtsave( filename = paste0( outdir, "/model_parameters_table.rtf" ) )
    
}

################################################################################
# End custom functions
################################################################################

# Set output directory
outdir <- 'out.02.multivariable'
dir.create(outdir, showWarnings = FALSE)

# read imputed data
imp <- readRDS( 'out.00.impute/mice_imputations.rds' )

set.seed( 321 )

# only run selection of not yet done
reduced_preds_file <- 'out.01.multivariable/remaining_preds_after_AIC.csv'

# check
if( !file.exists( reduced_preds_file ) )
{

    ########################
    # Univariable prediction
    ########################
    
    # List of potential predictors
    preds <- c("semiolmotor", "ageonset", "numaeds", "focal", 
           "newgen", "oldgen", "female", "fh", "neonatalszs", "fs", "gt9szs", 
           "priordc", "semiolgtc", "benign", "dd", "abnormalexam", "abnormalimaging", "eegabnormal", 
           "ageonsetcat", "eegepileptiform", "etiolunknown", "yearsszs", "yearsszfree", "dc", 
           "etiolbirthtrauma", "brainsurgery", "semiolmyoclonic", "numszstc_cat", "status", "impairing", "ddd",
           "semioltc", "semiolabsence", "etioltbi", "etiolinfxn", "etiolstructural", "yearsszfreecat", "etiolstroke", 
           "etioltumor", "etiolgenetic", "etiolmcd", "etiolstructuralother", "abnormalmri", "abnormalct", 
           "etiolencephalomalac", "etiolperiventleuk", "etiolatrophy")

    # Create a list to store pooled results
    significant_preds <- c()
    
    # Loop through each predictor for univariable analysis
    for( pred in preds ) 
    {
        print(paste("Running analysis for predictor:", pred))
        
        # Create the formula string dynamically
        formula_str <- paste0("survival::Surv(time = entry_time, time2 = exit_time, event) ~ ", pred)
        print(paste("Using formula string:", formula_str))
        
        # Fit univariable Cox models within imputed sets
        fit_models <- tryCatch({
            with(imp, {
                formula_obj <- as.formula(formula_str)
                survival::coxph(formula = formula_obj)
            })
        }, error = function(e) {
            print(paste("Error during with() for predictor:", pred))
            message("Error details:")
            message(e)
            return(NULL)
        })
        
        # Check if model fitting was successful before pooling
        if(!is.null(fit_models) && inherits(fit_models, "mira")) {
            # Pool the results from the models fitted to each imputed dataset
            pooled_summary <- tryCatch({
                summary(pool(fit_models))
            }, error = function(e) {
                print(paste("Error during pool() for predictor:", pred))
                message("Pooling error details:")
                message(e)
                return(NULL)
            })
            
            if(!is.null(pooled_summary)) {
                
                # Check for significance
                if(any(pooled_summary$p.value < 0.2)) {
                    significant_preds <- c(significant_preds, pred)
                }
            } else {
                print(paste("Skipping results storage for predictor:", pred, "due to pooling error."))
            }
        } else {
            print(paste("Skipping pooling for predictor:", pred, "due to error during model fitting."))
        }
    }
    
    
    # Print summary of univariable results
    cat("Predictors with p < 0.2:", paste(significant_preds, collapse = ", "), "\n")
    
    # Save the significant predictors for later use
    write.csv(data.frame(Predictor = significant_preds), file = paste0(outdir, '/significant_predictors.csv'), row.names = FALSE)
    




    #############################
    # Multivariable analysis
    #############################
    
    # Build formula for multivariable model
    multi_formula_str <- paste0("Surv(time = entry_time, time2 = exit_time, event) ~ ", paste(significant_preds, collapse = " + "))
    multi_formula_obj <- as.formula(multi_formula_str)
    
    # Fit multivariable Cox models
    multi_fit_models <- tryCatch({
        with(imp, {
            multi_formula_obj <- as.formula(multi_formula_str)
            survival::coxph(formula = multi_formula_obj)
        })
    }, error = function(e) {
        print("Error during multivariable model fitting")
        message("Error details:")
        message(e)
        return(NULL)
    })
    
    # Print parameters of multivariable model
    cat("\nFinal multivariable Cox model with predictors where p < 0.2 in univariable analysis:\n")
    multi_params <- parameters::parameters(multi_fit_models)
    print(multi_params)
    
    # Save the multivariable model parameters
    write.csv(as.data.frame(multi_params), file = paste0(outdir, '/multivariable_model_parameters.csv'), row.names = FALSE)
    
    #############################
    # Multivariable model pruning
    #############################
    
    # remove problematic preds
    sel_preds <- significant_preds[ !significant_preds %in% c( "gt9szs" ) ]
    
    # prune model (backward selection with AIC)
    pruned_model_aic <- backward_elimination_mi_aic(imp, sel_preds, trace = FALSE)
    
    # print remaining predictors
    print( remaining_preds <- pruned_model_aic$remaining )
    
    # Save the final predictors for next script
    write.csv( data.frame( pred = remaining_preds ), file = paste0(outdir, '/remaining_preds_after_AIC.csv'), row.names = FALSE)

} else {
    df <- read.csv( 'out.01.multivariable/remaining_preds_after_AIC.csv' ) 
    remaining_preds <- df$pred
}


#############################
# FINAL Multivariable model
#############################

# Build formula for multivariable model
final_formula_str <- paste0("Surv(time = entry_time, time2 = exit_time, event) ~ ", paste(remaining_preds, collapse = " + "))
final_formula_obj <- as.formula(final_formula_str)

# Fit multivariable Cox models
final_models <- tryCatch({
    with(imp, {
        final_formula_obj <- as.formula(final_formula_str)
        survival::coxph(formula = final_formula_obj)
    })
}, error = function(e) {
    print("Error during final multivariable model fitting")
    message("Error details:")
    message(e)
    return(NULL)
})

# check
parameters::parameters( final_models )

# save to nice table (RTF)
save_parameters_to_rtf( final_models )


#############################################
# Now run the leave-one-center-out validation
#############################################

# Get unique center values
centers <- unique(complete(imp, 1)$center)

# TODO remove (debug)
#imp_data <- imp
#formula_str <- final_formula_str
# TODO

# Run validation
validation_results <- leave_one_center_out_validation(
    imp_data = imp,
    formula_str = final_formula_str,
    centers = centers
)

# Create summary table
summary_table <- create_summary_table(validation_results)

# Print summary table
print(summary_table)

# Create plots
validation_plots <- plot_validation_results(validation_results)

# Display plots
print(validation_plots$c_index_plot)
print(validation_plots$auc_2yr_plot)
print(validation_plots$cal_slope_plot)

# save to disk
ggsave( plot = validation_plots$c_index_plot, file = paste0( outdir, '/c_index.png' ), dpi = 300, height = 6, width = 9, bg = 'white' )
ggsave( plot = validation_plots$auc_2yr_plot, file = paste0( outdir, '/auc_2yr.png' ), dpi = 300, height = 6, width = 9, bg = 'white' )
ggsave( plot = validation_plots$cal_slope_plot, file = paste0( outdir, '/cal_slope.png' ), dpi = 300, height = 6, width = 9, bg = 'white' )


# Calculate overall performance
overall_performance <- calculate_overall_performance(validation_results)
print(overall_performance)

# Save results
write.csv(summary_table, paste0( outdir, "/center_validation_results.csv" ), row.names = FALSE)
write.csv(overall_performance, paste0( outdir, "/overall_performance.csv" ), row.names = FALSE)

