#!/usr/bin/env Rscript
#
# Wim Otte (Apr/May 2025; w.m.otte@umcutrecht.nl)
#
# Build multivariable prediction model for epilepsy recurrence using
# multicenter data (IPD-meta-analysis).
#
# AIM: this script imputes the data.
#
################################################################################

# Load required packages
library('readr')
library('mice')
library('survival')
library('parameters')
library('dplyr') # Added for better data manipulation

################################################################################
# Start Custom Functions
################################################################################

###
# Get binary variables
###
is_binary_numeric <- function(df) {
    binary_cols <- sapply(df, function(col) {
        # Check if column is numeric
        if (!is.numeric(col)) {
            return(FALSE)
        }
        
        # Get unique values, excluding NA
        unique_vals <- unique(col[!is.na(col)])
        
        # Check if all values are 0 or 1
        all(unique_vals %in% c(0, 1)) && length(unique_vals) <= 2
    })
    
    # Return names of binary columns
    names(df)[binary_cols]
}

###
# Convert to factor (0,1)
###
convert_binary_to_factor <- function(data) {
    # Loop through each column in the data frame
    for (col_name in names(data)) {
        # Get unique values in the column (excluding NA)
        unique_vals <- unique(data[[col_name]][!is.na(data[[col_name]])])
        
        # Check if the column only contains 0 and 1
        if (length(unique_vals) <= 2 && all(unique_vals %in% c(0, 1))) {
            # Convert the column to a factor with appropriate labels
            data[[col_name]] <- factor(data[[col_name]], 
                                       levels = c(0, 1), 
                                       labels = c("no", "yes"))
        }
    }
    
    return(data)
}

###
# Remove single values columns
###
remove_single_value_columns <- function(data) {
    # Store column names to remove
    columns_to_remove <- character(0)
    
    # Loop through each column in the data frame
    for (col_name in names(data)) {
        # Get unique values including NA
        # We need to handle NA specially since unique() doesn't consider NAs identical
        if (all(is.na(data[[col_name]]))) {
            # Column contains only NAs
            columns_to_remove <- c(columns_to_remove, col_name)
            cat("Column '", col_name, "' contains only NA values and will be removed.\n", sep = "")
        } else {
            # Get unique non-NA values
            unique_vals <- unique(data[[col_name]][!is.na(data[[col_name]])])
            
            # Check if there's only one unique value
            if (length(unique_vals) == 1) {
                columns_to_remove <- c(columns_to_remove, col_name)
                cat("Column '", col_name, "' contains only the value: ", unique_vals, " and will be removed.\n", sep = "")
            }
        }
    }
    
    # Remove the identified columns
    if (length(columns_to_remove) > 0) {
        data <- data[, !(names(data) %in% columns_to_remove), drop = FALSE]
    }
    
    return(data)
}

###
# Remove imbalanced binary variables
###
remove_imbalanced_binary <- function(data, min_threshold = 10) {
    # Store column names to remove
    columns_to_remove <- character(0)
    
    # Loop through each column in the data frame
    for (col_name in names(data)) {
        # Skip non-binary columns - only process numeric, logical, or factor columns that could be binary
        if (!(is.numeric(data[[col_name]]) || is.logical(data[[col_name]]) || is.factor(data[[col_name]]))) {
            next
        }
        
        # Get frequency table of values
        val_counts <- table(data[[col_name]], useNA = "no")  # Exclude NAs from the check
        
        # Check if it's a binary variable (has exactly 2 unique values)
        if (length(val_counts) == 2) {
            # Check if any value appears less than the threshold
            if (any(val_counts < min_threshold)) {
                columns_to_remove <- c(columns_to_remove, col_name)
                
                # Format output message
                imbalanced_value <- names(val_counts)[which.min(val_counts)]
                count <- min(val_counts)
                
                cat("Column '", col_name, "' is imbalanced: value '", 
                    imbalanced_value, "' appears only ", count, 
                    " times (threshold: ", min_threshold, "). Column will be removed.\n", sep = "")
            }
        }
    }
    
    # Remove the identified columns
    if (length(columns_to_remove) > 0) {
        data <- data[, !(names(data) %in% columns_to_remove), drop = FALSE]
    }
    
    return(data)
}

###
# Get and preprocess data
###
get_data <- function() {
    # read CSV
    data <- read.csv('data/stackedshort_wim__May_01-2025.csv', stringsAsFactors = TRUE)
    
    # very high correlation (85-95%!) with ageonset, so remove (to prevent fitting issues later on)
    data$agestartfu <- NULL
    data$agelastsz <- NULL
   
    # Create a new categorical variable as numsztc has a difficult distribution
    data$numszstc_cat <- cut(data$numszstc, 
                             breaks = c(0, 5, 10, Inf), 
                             labels = c("0-5", "5-10", "10+"),
                             right = FALSE,
                             include.lowest = TRUE)
    
    data$numszstc_cat <- as.factor(data$numszstc_cat)
    data$numszstc <- NULL
    
    # Process EEG abnormal
    data$eegabnormal <- as.character(data$eegabnormal)
    data$eegabnormal[data$eegabnormal == "0. No"] <- "0"
    data$eegabnormal[data$eegabnormal == "1. Yes"] <- "1"
    data$eegabnormal[data$eegabnormal == ""] <- NA
    data$eegabnormal[!data$eegabnormal %in% c("0", "1")] <- NA
    data$eegabnormal <- as.numeric(data$eegabnormal)
    
    # Process EEG epileptiform
    data$eegepileptiform <- as.character(data$eegepileptiform)
    data$eegepileptiform[data$eegepileptiform == "0. No"] <- "0"
    data$eegepileptiform[data$eegepileptiform == "1. Yes"] <- "1"
    data$eegepileptiform[data$eegepileptiform == ""] <- NA
    data$eegepileptiform[!data$eegepileptiform %in% c("0", "1")] <- NA
    data$eegepileptiform <- as.numeric(data$eegepileptiform)
    
    # Process abnormal MRI
    data$abnormalmri <- as.character(data$abnormalmri)
    data$abnormalmri[data$abnormalmri == "0. No"] <- "0"
    data$abnormalmri[data$abnormalmri == "1. Yes"] <- "1"
    data$abnormalmri[data$abnormalmri == ""] <- NA
    data$abnormalmri[!data$abnormalmri %in% c("0", "1")] <- NA
    data$abnormalmri <- as.numeric(data$abnormalmri)
    
    # Process abnormal CT
    data$abnormalct <- as.character(data$abnormalct)
    data$abnormalct[data$abnormalct == "0. No"] <- "0"
    data$abnormalct[data$abnormalct == "1. Yes"] <- "1"
    data$abnormalct[data$abnormalct == ""] <- NA
    data$abnormalct[!data$abnormalct %in% c("0", "1")] <- NA
    data$abnormalct <- as.numeric(data$abnormalct)
    
    # Process age onset category
    data$ageonsetcat <- as.character(data$ageonsetcat)
    data$ageonsetcat[data$ageonsetcat == "0. Child"] <- "0_child"
    data$ageonsetcat[data$ageonsetcat == "1. Adolescent"] <- "1_adolescent"
    data$ageonsetcat[data$ageonsetcat == "2. Adult"] <- "2_adult"
    data$ageonsetcat[data$ageonsetcat == ""] <- NA
    data$ageonsetcat <- as.factor(data$ageonsetcat)
    
    # Remove negative exit_times [n=3]
    data <- data[data$exit_time > data$entry_time, ]
    
    # Fix focal variable that has n=4 '2' values
    data[!is.na(data$focal) & data$focal == 2, 'focal'] <- NA
    
    # Fix abnormalimaging variable that has n=2 '2' values
    data[!is.na(data$abnormalimaging) & data$abnormalimaging == 2, 'abnormalimaging'] <- NA
    
    # Remove single value columns
    data <- remove_single_value_columns(data)
    
    # Convert to binary factor
    data <- convert_binary_to_factor(data)
    
    # Remove imbalanced binary variables
    data <- remove_imbalanced_binary(data)
    
    # Add entry time
    data$entry_time <- 0
    
    # Convert event to 0/1
    data$event <- as.numeric(data$event) - 1
    
    # Rename
    data$center <- data$pub
    data$pub <- NULL
    
    return(data)
}

###
# Check if no NA's remain
###
all_complete <- function(df) {
    # Check all columns at once using complete.cases()
    all(complete.cases(df))
}


################################################################################
# End custom functions
################################################################################

# Set output directory
outdir <- 'out.00.impute'
dir.create(outdir, showWarnings = FALSE)

# Get preprocessed data
data <- get_data()

# Write data to output file
readr::write_tsv( data, file = paste0( outdir, '/data.tsv' ), quote = 'all' )

# Run mice for multiple imputation
m <- 10  # Number of imputations
imp <- mice(data, m = m, seed = 123)

# Save the imputation to disk
saveRDS( imp, file = paste0( outdir, "/mice_imputations.rds" ) )

# Check if NAs are gone in all imputed datasets
for(i in 1:m) {
    print(paste("Checking imputation", i))
    if(!all_complete(mice::complete(imp, i))) {
        stop(paste0("*** ERROR ***: NA's detected after imputation! -> imputation ", i))
    }
}

# Save complete datasets from multiple imputations for later processing
# Fixed to save each imputation to its corresponding file
dfc1 <- mice::complete(imp, 1)
dfc2 <- mice::complete(imp, 2)
dfc3 <- mice::complete(imp, 3)

# Write to disk with correct filenames
readr::write_tsv(dfc1, file = paste0(outdir, '/complete_01.tsv'), quote = 'all')
readr::write_tsv(dfc2, file = paste0(outdir, '/complete_02.tsv'), quote = 'all')
readr::write_tsv(dfc3, file = paste0(outdir, '/complete_03.tsv'), quote = 'all')

# check highly-correlated vars

# get complete set nr. 1 (optional, for checking structure)
data1 <- mice::complete( imp, 1 )

# Get numeric columns
nums <- unlist( lapply( data1, is.numeric ) )
num_data <- data1[, nums]

# Calculate correlation matrix
corr_matrix <- cor(num_data, use = "pairwise.complete.obs")

# Find high correlations
high_corrs <- which(abs(corr_matrix) > 0.75 & abs(corr_matrix) < 1, arr.ind = TRUE)

# conclusion: remove 'agestartfu' and 'agelastsz' (but keep ageonset)
#         var1       var2 correlation
# 1 agestartfu   ageonset   0.7988052
# 2  agelastsz   ageonset   0.8598004
# 3   ageonset agestartfu   0.7988052
# 4  agelastsz agestartfu   0.9547565
# 5   ageonset  agelastsz   0.8598004
# 6 agestartfu  agelastsz   0.9547565
var_pairs <- data.frame(
    var1 = rownames(corr_matrix)[high_corrs[,1]],
    var2 = colnames(corr_matrix)[high_corrs[,2]],
    correlation = corr_matrix[high_corrs]
)

# print high correlations
print( var_pairs )

