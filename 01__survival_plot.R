# Load required libraries
library('survival')
library('survminer')
library('dplyr')
library('ggplot2')
library('RColorBrewer')
library('forestplot')

###############################################################################

# outdir 
outdir <- 'out.01.suvival.plot'
dir.create( outdir, showWarnings = FALSE )

# get (non-imputed) data
data <- readr::read_tsv( 'out.00.impute/data.tsv', show_col_types = FALSE )

# Convert center to a factor (if it's not already)
data$center <- factor(data$center)

# Fit survival model for all centers combined
surv_all <- survfit(Surv(exit_time, event) ~ 1, data = data)
print( surv_all )

# Fit survival model by center
surv_by_center <- survfit(Surv(exit_time, event) ~ center, data = data)

# Create a nice color palette for the 15 centers
my_colors <- colorRampPalette(brewer.pal(12, "Set3"))( length( unique( data$center ) )) # n= 15

 
###

# First, get the data in the right format using survfit and ggsurv
fit <- survfit(Surv(exit_time, event) ~ center, data = data)
surv_data <- surv_summary(fit, data = data)

# Now create the plot with ggplot2 and facet_wrap
p1 <- ggplot(surv_data, aes(x = time, y = surv, color = strata)) +
    geom_vline(xintercept = 12, linetype = "dashed", color = "gray80" ) + # 1 year
    geom_vline(xintercept = 48, linetype = "dashed", color = "gray80" ) + # 4 years
    geom_step() +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = strata), alpha = 0.3, color = NA) +
    facet_wrap(~ strata, ncol = 5 ) +
    labs( x = "Time (months)", y = "Survival Probability") +
    theme_bw() +
    scale_x_continuous( breaks = c( 12, 48, 100, 200, 300 ) ) +
    theme(legend.position = "none")  # Remove legend since facets have labels

# Save plots (uncomment and modify as needed)
ggsave( paste0( outdir, "/survival_curves_per_center.png" ), plot = p1, dpi = 600, width = 9, height = 6, bg = 'white')

# Forest plot for comparing centers
# First, create a Cox model
cox_model <- coxph(Surv(entry_time, event) ~ center, data = data)
summary(cox_model)

# Extract the model information
cox_summary <- summary(cox_model)
hr_data <- as.data.frame(cox_summary$conf.int)
hr_data$p.value <- cox_summary$coefficients[, "Pr(>|z|)"]
hr_data$variable <- rownames(hr_data)

# Prepare the forestplot data
tabletext <- cbind(
    c("Variable", hr_data$variable),
    c("HR (95% CI)", 
      sprintf("%.2f (%.2f-%.2f)", 
              hr_data$`exp(coef)`, 
              hr_data$`lower .95`, 
              hr_data$`upper .95`)),
    c("p-value", 
      sprintf("%.3f", hr_data$p.value))
)




# Create the tabletext with proper alignment
tabletext <- cbind(
    c("Variable", hr_data$variable),
    c("HR (95% CI)", 
      sprintf("%.2f (%.2f-%.2f)", 
              hr_data$`exp(coef)`, 
              hr_data$`lower .95`, 
              hr_data$`upper .95`)),
    c("p-value", 
      sprintf("%.3f", hr_data$p.value))
)

# Create dummy data for the header row
header_row <- data.frame(
    mean = NA,
    lower = NA,
    upper = NA
)

# Combine with the real data
plot_data <- rbind(
    header_row,
    data.frame(
        mean = hr_data$`exp(coef)`,
        lower = hr_data$`lower .95`,
        upper = hr_data$`upper .95`
    )
)

# Set up the PNG device
png(filename = paste0( outdir, "/forest_plot_per_center.png" ), 
    width = 1700,     # Width in pixels
    height = 1400,     # Height in pixels
    res = 300)        # Resolution in ppi (pixels per inch)

# Now create the forestplot
forestplot(tabletext,
           mean = plot_data$mean,
           lower = plot_data$lower,
           upper = plot_data$upper,
           is.summary = c(TRUE, rep(FALSE, nrow(hr_data))),
           zero = 1,
           boxsize = 0.2,
           col = fpColors(box = "royalblue",
                          line = "darkblue",
                          summary = "royalblue"),
           xlab = "Hazard Ratio",
           title = "Forest Plot - Hazard Ratios by Center",
           hrzl_lines = list("2" = gpar(lwd = 1, col = "#444444")))

# Close the device and save the file
dev.off()


# For publication-quality plots, you may want to add tables with Cox proportional hazards results
# Extract hazard ratios and confidence intervals
cox_summary <- summary(cox_model)
hr_table <- data.frame(
    Center = levels(data$center),
    HR = c(1, exp(cox_summary$coefficients[, "coef"])),
    Lower_CI = c(1, exp(cox_summary$coefficients[, "coef"] - 1.96 * cox_summary$coefficients[, "se(coef)"])),
    Upper_CI = c(1, exp(cox_summary$coefficients[, "coef"] + 1.96 * cox_summary$coefficients[, "se(coef)"])),
    P_value = c(NA, cox_summary$coefficients[, "Pr(>|z|)"])
)
print(hr_table)

# Create a table with formatted values
hr_table$HR_CI <- paste0(
    sprintf("%.2f", hr_table$HR), " (",
    sprintf("%.2f", hr_table$Lower_CI), "-",
    sprintf("%.2f", hr_table$Upper_CI), ")"
)
hr_table$P_formatted <- ifelse(
    is.na(hr_table$P_value), "Reference",
    ifelse(hr_table$P_value < 0.001, "<0.001",
           sprintf("%.3f", hr_table$P_value)
    )
)

rownames( hr_table ) <- NULL
outd <- hr_table[, c("Center", "HR_CI", "P_formatted")]
readr::write_tsv( outd, file = paste0( outdir, '/center_hazards.tsv' ), quote = 'all' )

