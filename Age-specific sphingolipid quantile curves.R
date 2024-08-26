#==========================================================================================
# Age-specific sphingolipid quantile curves
# Author: Denis Infanger
#==========================================================================================

#------------------------------------------------------------------------------------------
# Load packages
#------------------------------------------------------------------------------------------

library(gamlss)
library(ggplot2)
library(reshape2)
library(directlabels)
library(readxl)
library(car)
library(broom)
library(dplyr)
library(tidyr)
library(ggdendro)
library(egg)
library(readxl)
library(readr)
library(tibble) 
library(writexl)
library(data.table)

rstudioapi::writeRStudioPreference("data_viewer_max_columns", 1000L)

#------------------------------------------------------------------------------------------
# Define paths
#------------------------------------------------------------------------------------------

# Paths
setwd("") # Main path

data_path <- "./data" # Path for data
graphics_path <- "./output/graphics" # Path for graphics
text_path <- "./output/text" # Path for text-output
diagnostics_path <- "./output/graphics/diagnostics" # Path for graphics

#------------------------------------------------------------------------------------------
# Import data
#------------------------------------------------------------------------------------------

#import_data
dat <- as.data.frame(read_excel(paste0(data_path, "/", "xxx.xlsx")))

# View(dat)
head(dat)
str(dat)

#------------------------------------------------------------------------------------------
# Reshape to long
#------------------------------------------------------------------------------------------
(vars_to_plot <- names(dat)[c(15:43)])

dat_long <- reshape2::melt(
  dat
  , id.vars = c("Age", "ID", "Sex")
  , measure.vars = vars_to_plot
)

#------------------------------------------------------------------------------------------
# 1st PART : MALES ONLY
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
# Males: fit GAMLSS-model
#------------------------------------------------------------------------------------------

# Assuming you have already created the 'dat_long' data frame as you did before

# Step 1: Get unique variable names (dependent variables)
dependent_vars <- unique(dat_long$variable)

# Step 2: Fit the GAMLSS models using the loop
# Create an empty list to store the fitted GAMLSS models
fitted_models_list_m <- list()

# Loop through each dependent variable and fit the GAMLSS model
for (dep_var in dependent_vars) {
  
  cat("\nSphingolipid:", dep_var, "\n")
  
  # Subset data frame
  
  dat_nona_m <- subset(dat_long, variable %in% dep_var & Sex == 0)
  dat_nona_m <- dat_nona_m[complete.cases(dat_nona_m$value, dat_nona_m$Age), ]
  
  k <- log(nrow(dat_nona_m))  # Using the BIC for model selection. Use k = 2 for AIC.
  
  mod_m <- lms(
    y = value,  # Dependent variable for the current iteration
    x = Age,    # Independent variable for the current iteration
    data = dat_nona_m,  # Subset the data for the current dependent variable
    method.pb = "GAIC",
    k = k,
    cent = c(),  # Percentiles to display (does make no difference for the model fit)
    families = .realplus[!(.realplus %in% c("GIG", "IG", "LNO"))],        # Try all distributions for positive variables (for other possibilities, see ?fitDist)
    # families = "BCT",
    trans.x = FALSE,  # Set this to TRUE to transform the x-variable.
    calibration = FALSE
  )
  
  # Add the fitted model to the list
  fitted_models_list_m[[dep_var]] <- mod_m
}

# Now fitted_models_list contains the fitted GAMLSS models for each of the dependent variables.
# You can access them using fitted_models_list[['variable_name']], where 'variable_name' is the name of the dependent variable.

#------------------------------------------------------------------------------------------
# Males: loop through all species, save diagnostics and plots for each separately 
#------------------------------------------------------------------------------------------

#Introduce the correct sphingolipid names
sphingo_names <- c(
  "SM(d18:1/18:1)",
  "Cer(d18:1/16:0)",
  "SM(d18:1/16:0)",
  "Cer(d18:1/24:1)",
  "CerP(d18:1/16:0)",
  "Sph(d18:1)",
  "Cer(d18:0/16:0)",
  "Sph1P(d18:1)",
  "Cer(d18:1/18:0)", 
  "Spa1P(d18:0)",
  "SM(d18:1/12:0)", 
  "Cer(d18:1/10:0)", 
  "Cer(d18:2/24:1)",
  "Cer(d18:1/18:1)", 
  "GlcSph(18:1)",
  "Cer(d18:2/24:0)",
  "Cer(d18:1/20:0)",
  "SM(d18:1/18:0)",
  "Cer(d18:1/24:0)", 
  "GlcCer(d18:1/16:0)", 
  "Cer(d18:0/22:0) or Cer(m18:12/2:0)", 
  "LacCer(d18:1/16:0)",
  "GlcCer(d18:1/24:1)", 
  "GlcCer(d18:1/18:0)",
  "Cer(d18:1/12:0)",
  "Cer(d18:1/22:0)",
  "Cer(d18:1/16:0) / Cer(d18:1/24:0)",
  "Cer(d18:1/18:0) /  Cer(d18:1/24:0)",
  "Cer(d18:1/24:1) / Cer(d18:1/24:0)"
)


for (i in seq_along(fitted_models_list_m)) {
  
  # Save lipid name for file names etc.
  
  species_name <- names(fitted_models_list_m)[i]
  
  # Clean up name
  
  species_name_clean <- gsub("/", "_", species_name) 
  species_name_clean <- gsub("\\(", "_", species_name_clean) 
  species_name_clean <- gsub("\\)", "_", species_name_clean) 
  species_name_clean <- gsub("\\:", "_", species_name_clean) 
  
  cat("\nSphingolipid:", species_name, "\n")
  
  #------------------------------------------------------------------------------------------
  # Males: diagnostics
  #------------------------------------------------------------------------------------------
  
  png(paste0(diagnostics_path, "/",species_name_clean, "diagnostics_males.png"), width = 10, height = 7, units = "in", res = 300)
  wp(fitted_models_list_m[[i]])
  dev.off()
  
  #------------------------------------------------------------------------------------------
  # Males: Save centiles as Excel-table
  #------------------------------------------------------------------------------------------
  
  mod_m <- fitted_models_list_m[[i]]
  
  # Which centiles to save
  centiles_pred <- c(3, 15, 50, 85, 97)
  
  # Predict centiles and save them in a data frame
  predmat_mod_m <- centiles.pred(
    obj = mod_m
    , type = "centiles"
    , xname = "Age"
    , xvalues = seq(min(dat$Age), max(dat$Age), by = 1) # Only for each age in 1-year-steps
    , cent = centiles_pred
    , calibration = TRUE
  )
  
  # Rename the variables in the data frame
  names(predmat_mod_m)[1] <- "Age"
  names(predmat_mod_m)[2:(length(centiles_pred) + 1)] <- paste0(centiles_pred, " P")
  
  # Save Excel-sheet with centiles
  write_xlsx(predmat_mod_m, path = paste0(text_path, "/", species_name_clean,  "_males.xlsx"))
  
  # Clean up
  rm(predmat_mod_m)
  
  #------------------------------------------------------------------------------------------
  # Males: plot centiles
  #------------------------------------------------------------------------------------------
  
  # Which centiles to plot
  centiles_pred <- c(3, 15, 50, 85, 97)
  
  # Predict centiles and save them in a data frame
  predmat_mod_m <- centiles.pred(
    obj = mod_m
    , type = "centiles"
    , xname = "Age"
    , xvalues = seq(min(dat$Age), max(dat$Age), by = 0.1) # For which x-values (age)
    , cent = centiles_pred
    , calibration = TRUE
  )
  
  # Rename the variables in the data frame
  names(predmat_mod_m)[1] <- "Age"
  names(predmat_mod_m)[2:(length(centiles_pred) + 1)] <- paste0(centiles_pred, " P")
  
  # Reshape centile-dataframe from wide to long format for plotting
  
  predmat_mod_long_m <- reshape2::melt(
    predmat_mod_m
    , id.vars = c("Age")
  )
  
  # Create the plot
  
  theme_set(theme_bw())
  p <- ggplot(predmat_mod_long_m, aes(x = Age, y = value, group = variable)) +
    geom_point(data = subset(dat_long, variable %in% species_name & Sex == 0), aes(y = value, x = Age), inherit.aes = FALSE, size = 2, alpha = 0.5) + # Subset data to include "SM(d18:1/16:0)" only
    geom_line(aes(colour = variable), linewidth = 1) +
    xlab("Age (years)") +
    ylab(paste0(sphingo_names[i])) +
    ggtitle(paste0("N = ", nobs(mod_m), ", distribution = ", mod_m$family[2])) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_x_continuous(breaks = seq(0, 100, 5)) +
    scale_color_brewer(palette = "Dark2", name = "Percentile") +
    theme(
      axis.title.y = element_text(size = 17, hjust = 0.5)
      ,  axis.title.x = element_text(size = 17, hjust = 0.5)
      , legend.position = "right"
      , axis.text.x = element_text(colour = "black", size=14)
      , axis.text.y = element_text(colour = "black", size=15)
      , legend.text = element_text(size= 15)
      , plot.title = element_text(face = "bold")
      , legend.title = element_text(size = 15)
    )
  
  # p
  
  # Save graphic
  ggsave(paste0(graphics_path, "/", species_name_clean, "model_percentiles_males.png"), p, width = 20*0.6, height = 14*0.6, units = "in", dpi = 300, type = "cairo-png")
  
  #------------------------------------------------------------------------------------------
  # Males: plot publication-ready graph
  #------------------------------------------------------------------------------------------
  
  # Set the line indicators so that different centiles are plotted with different linetypes
  predmat_mod_long_m$line_ind <- NA
  predmat_mod_long_m$line_ind[predmat_mod_long_m$variable %in% c("50 P")] <- 0 # Plot 50 P with a continuous line
  predmat_mod_long_m$line_ind[predmat_mod_long_m$variable %in% c("15 P", "85 P")] <- 1 # Plot 15 P and 85 P with a dashed line
  predmat_mod_long_m$line_ind[predmat_mod_long_m$variable %in% c("3 P", "97 P")] <- 2 # Plot 3 P and 97 P with a dotted line
  
  predmat_mod_long_m$line_ind <- factor(predmat_mod_long_m$line_ind)
  
  # Rename levels of percentiles
  levels(predmat_mod_long_m$variable) <- gsub(" P", "th", gsub("3 P", "3rd", levels(predmat_mod_long_m$variable)))
  
  my.dl <- list(box.color=NA, "draw.rects")
  
  theme_set(theme_bw())
  p_publication <- ggplot(predmat_mod_long_m, aes(x = Age, y = value, group = variable)) +
    geom_line(aes(linetype = line_ind), linewidth = 1) +
    geom_dl(aes(label = variable), method = list(list("last.points", "calc.boxes", "enlarge.box", my.dl), dl.trans(x = x + 0.2), cex = 1.2)) +
    xlab("Age (years)") +
    ylab(paste0(sphingo_names[i])) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_x_continuous(breaks = seq(0, 90, 5)) + # Format x-axis
    scale_color_brewer(palette = "Dark2", name = "Percentile") +
    theme(
      axis.title.y.left = element_text(size = 17, hjust = 0.5, margin = margin(t = 0, r = 10, b = 0, l = 0))
      , axis.title.y.right = element_text(size = 17, hjust = 0.5, margin = margin(t = 0, r = 0, b = 0, l = 15))
      ,  axis.title.x = element_text(size = 17, hjust = 0.5)
      , legend.position = "none"
      , axis.text.x = element_text(colour = "black", size=14)
      , axis.text.y = element_text(colour = "black", size=15)
      , legend.text = element_text(size= 15)
      , plot.title = element_text(face = "bold")
      , legend.title = element_text(size = 15)
      , panel.grid.minor = element_line(colour = "gray80")
      , panel.grid.major = element_line(colour = "gray80")
    )
  
  # Save graphic
  ggsave(paste0(graphics_path, "/", species_name_clean, "model_percentiles_final_males.png"), p_publication, width = 20*0.6, height = 14*0.6, units = "in", dpi = 300, type = "cairo-png")
  
}

#------------------------------------------------------------------------------------------
# 2nd PART : FEMALES ONLY
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
# Females: fit GAMLSS-model
#------------------------------------------------------------------------------------------

# Assuming you have already created the 'dat_long' data frame as you did before

# Step 1: Get unique variable names (dependent variables)
dependent_vars <- unique(dat_long$variable)

# Step 2: Fit the GAMLSS models using the loop
# Create an empty list to store the fitted GAMLSS models
fitted_models_list_f <- list()

# Loop through each dependent variable and fit the GAMLSS model
for (dep_var in dependent_vars) {
  
  cat("\nSphingolipid:", dep_var, "\n")
  
  # Subset data frame
  
  dat_nona_f <- subset(dat_long, variable %in% dep_var & Sex == 1)
  dat_nona_f <- dat_nona_f[complete.cases(dat_nona_f$value, dat_nona_f$Age), ]
  
  k <- log(nrow(dat_nona_f))  # Using the BIC for model selection. Use k = 2 for AIC.
  
  mod_f <- lms(
    y = value,  # Dependent variable for the current iteration
    x = Age,    # Independent variable for the current iteration
    data = dat_nona_f,  # Subset the data for the current dependent variable
    method.pb = "GAIC",
    k = k,
    cent = c(),  # Percentiles to display (does make no difference for the model fit)
    families = .realplus[!(.realplus %in% c("GIG", "IG", "LNO"))],        # Try all distributions for positive variables (for other possibilities, see ?fitDist)
    # families = "BCT",
    trans.x = FALSE,  # Set this to TRUE to transform the x-variable.
    calibration = FALSE
  )
  
  # Add the fitted model to the list
  fitted_models_list_f[[dep_var]] <- mod_f
}

# Now fitted_models_list contains the fitted GAMLSS models for each of the dependent variables.
# You can access them using fitted_models_list[['variable_name']], where 'variable_name' is the name of the dependent variable.

#------------------------------------------------------------------------------------------
# Females: loop through all species, save diagnostics and plots for each separately 
#------------------------------------------------------------------------------------------

#Introduce the correct sphingolipid names
sphingo_names <- c(
  "SM(d18:1/18:1)",
  "Cer(d18:1/16:0)",
  "SM(d18:1/16:0)",
  "Cer(d18:1/24:1)",
  "CerP(d18:1/16:0)",
  "Sph(d18:1)",
  "Cer(d18:0/16:0)",
  "Sph1P(d18:1)",
  "Cer(d18:1/18:0)", 
  "Spa1P(d18:0)",
  "SM(d18:1/12:0)", 
  "Cer(d18:1/10:0)", 
  "Cer(d18:2/24:1)",
  "Cer(d18:1/18:1)", 
  "GlcSph(18:1)",
  "Cer(d18:2/24:0)",
  "Cer(d18:1/20:0)",
  "SM(d18:1/18:0)",
  "Cer(d18:1/24:0)", 
  "GlcCer(d18:1/16:0)", 
  "Cer(d18:0/22:0) or Cer(m18:12/2:0)", 
  "LacCer(d18:1/16:0)",
  "GlcCer(d18:1/24:1)", 
  "GlcCer(d18:1/18:0)",
  "Cer(d18:1/12:0)",
  "Cer(d18:1/22:0)",
  "Cer(d18:1/16:0) / Cer(d18:1/24:0)",
  "Cer(d18:1/18:0) /  Cer(d18:1/24:0)",
  "Cer(d18:1/24:1) / Cer(d18:1/24:0)"
)


for (i in seq_along(fitted_models_list_f)) {
  
  # Save lipid name for file names etc.
  
  species_name <- names(fitted_models_list_f)[i]
  
  # Clean up name
  
  species_name_clean <- gsub("/", "_", species_name) 
  species_name_clean <- gsub("\\(", "_", species_name_clean) 
  species_name_clean <- gsub("\\)", "_", species_name_clean) 
  species_name_clean <- gsub("\\:", "_", species_name_clean) 
  
  cat("\nSphingolipid:", species_name, "\n")
  
  #------------------------------------------------------------------------------------------
  # Females: diagnostics
  #------------------------------------------------------------------------------------------
  
  png(paste0(diagnostics_path, "/",species_name_clean, "diagnostics_females.png"), width = 10, height = 7, units = "in", res = 300)
  wp(fitted_models_list_f[[i]])
  dev.off()
  
  #------------------------------------------------------------------------------------------
  # Females: Save centiles as Excel-table
  #------------------------------------------------------------------------------------------
  
  mod_f <- fitted_models_list_f[[i]]
  
  # Which centiles to save
  centiles_pred <- c(3, 15, 50, 85, 97)
  
  # Predict centiles and save them in a data frame
  predmat_mod_f <- centiles.pred(
    obj = mod_f
    , type = "centiles"
    , xname = "Age"
    , xvalues = seq(min(dat$Age), max(dat$Age), by = 1) # Only for each age in 1-year-steps
    , cent = centiles_pred
    , calibration = TRUE
  )
  
  # Rename the variables in the data frame
  names(predmat_mod_f)[1] <- "Age"
  names(predmat_mod_f)[2:(length(centiles_pred) + 1)] <- paste0(centiles_pred, " P")
  
  # Save Excel-sheet with centiles
  write_xlsx(predmat_mod_f, path = paste0(text_path, "/", species_name_clean,  "_females.xlsx"))
  
  # Clean up
  rm(predmat_mod_f)
  
  #------------------------------------------------------------------------------------------
  # Females: plot centiles
  #------------------------------------------------------------------------------------------
  
  # Which centiles to plot
  centiles_pred <- c(3, 15, 50, 85, 97)
  
  # Predict centiles and save them in a data frame
  predmat_mod_f <- centiles.pred(
    obj = mod_f
    , type = "centiles"
    , xname = "Age"
    , xvalues = seq(min(dat$Age), max(dat$Age), by = 0.1) # For which x-values (age)
    , cent = centiles_pred
    , calibration = TRUE
  )
  
  # Rename the variables in the data frame
  names(predmat_mod_f)[1] <- "Age"
  names(predmat_mod_f)[2:(length(centiles_pred) + 1)] <- paste0(centiles_pred, " P")
  
  # Reshape centile-dataframe from wide to long format for plotting
  
  predmat_mod_long_f <- reshape2::melt(
    predmat_mod_f
    , id.vars = c("Age")
  )
  
  # Create the plot
  
  theme_set(theme_bw())
  p <- ggplot(predmat_mod_long_f, aes(x = Age, y = value, group = variable)) +
    geom_point(data = subset(dat_long, variable %in% species_name & Sex == 1), aes(y = value, x = Age), inherit.aes = FALSE, size = 2, alpha = 0.5) + # Subset data to include "SM(d18:1/16:0)" only
    geom_line(aes(colour = variable), linewidth = 1) +
    xlab("Age (years)") +
    ylab(paste0(sphingo_names[i])) +
    ggtitle(paste0("N = ", nobs(mod_f), ", distribution = ", mod_f$family[2])) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_x_continuous(breaks = seq(0, 100, 5)) +
    scale_color_brewer(palette = "Dark2", name = "Percentile") +
    theme(
      axis.title.y = element_text(size = 17, hjust = 0.5)
      ,  axis.title.x = element_text(size = 17, hjust = 0.5)
      , legend.position = "right"
      , axis.text.x = element_text(colour = "black", size=14)
      , axis.text.y = element_text(colour = "black", size=15)
      , legend.text = element_text(size= 15)
      , plot.title = element_text(face = "bold")
      , legend.title = element_text(size = 15)
    )
  
  # p
  
  # Save graphic
  ggsave(paste0(graphics_path, "/", species_name_clean, "model_percentiles_females.png"), p, width = 20*0.6, height = 14*0.6, units = "in", dpi = 300, type = "cairo-png")
  
  #------------------------------------------------------------------------------------------
  # Females: plot publication-ready graph
  #------------------------------------------------------------------------------------------
  
  # Set the line indicators so that different centiles are plotted with different linetypes
  predmat_mod_long_f$line_ind <- NA
  predmat_mod_long_f$line_ind[predmat_mod_long_f$variable %in% c("50 P")] <- 0 # Plot 50 P with a continuous line
  predmat_mod_long_f$line_ind[predmat_mod_long_f$variable %in% c("15 P", "85 P")] <- 1 # Plot 15 P and 85 P with a dashed line
  predmat_mod_long_f$line_ind[predmat_mod_long_f$variable %in% c("3 P", "97 P")] <- 2 # Plot 3 P and 97 P with a dotted line
  
  predmat_mod_long_f$line_ind <- factor(predmat_mod_long_f$line_ind)
  
  # Rename levels of percentiles
  levels(predmat_mod_long_f$variable) <- gsub(" P", "th", gsub("3 P", "3rd", levels(predmat_mod_long_f$variable)))
  
  my.dl <- list(box.color=NA, "draw.rects")
  
  theme_set(theme_bw())
  p_publication <- ggplot(predmat_mod_long_f, aes(x = Age, y = value, group = variable)) +
    geom_line(aes(linetype = line_ind), linewidth = 1) +
    geom_dl(aes(label = variable), method = list(list("last.points", "calc.boxes", "enlarge.box", my.dl), dl.trans(x = x + 0.2), cex = 1.2)) +
    xlab("Age (years)") +
    ylab(paste0(sphingo_names[i])) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_x_continuous(breaks = seq(0, 90, 5)) + # Format x-axis
    scale_color_brewer(palette = "Dark2", name = "Percentile") +
    theme(
      axis.title.y.left = element_text(size = 17, hjust = 0.5, margin = margin(t = 0, r = 10, b = 0, l = 0))
      , axis.title.y.right = element_text(size = 17, hjust = 0.5, margin = margin(t = 0, r = 0, b = 0, l = 15))
      ,  axis.title.x = element_text(size = 17, hjust = 0.5)
      , legend.position = "none"
      , axis.text.x = element_text(colour = "black", size=14)
      , axis.text.y = element_text(colour = "black", size=15)
      , legend.text = element_text(size= 15)
      , plot.title = element_text(face = "bold")
      , legend.title = element_text(size = 15)
      , panel.grid.minor = element_line(colour = "gray80")
      , panel.grid.major = element_line(colour = "gray80")
    )
  
  # Save graphic
  ggsave(paste0(graphics_path, "/", species_name_clean, "model_percentiles_final_females.png"), p_publication, width = 20*0.6, height = 14*0.6, units = "in", dpi = 300, type = "cairo-png")
  
}
