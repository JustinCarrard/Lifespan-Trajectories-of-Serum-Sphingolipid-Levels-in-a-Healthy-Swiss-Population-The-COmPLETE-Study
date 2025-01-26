#==========================================================================================
# Age-specific sphingolipid quantile curves
# Author: Denis Infanger & Justin Carrard
#==========================================================================================

#------------------------------------------------------------------------------------------
# Load packages
#------------------------------------------------------------------------------------------

library(readxl)
library(car)
library(broom)
library(dplyr)
library(tidyr)
library(ggdendro)
library(egg)
library(ggplot2)
library(reshape2)
library(readxl)
library(readr)
library(tibble) 
library(splines)
library(writexl)
library(kableExtra)
library(flextable) 
library(knitr)  
library(writexl)
library(data.table)
library(chron)
library(nlme)
library(mgcv)

rstudioapi::writeRStudioPreference("data_viewer_max_columns", 1000L)

#------------------------------------------------------------------------------------------
# Define paths
#------------------------------------------------------------------------------------------

# Paths
setwd("") # Main path

data_path <- "./data" # Path for data
graphics_path <- "./output/graphics" # Path for graphics
text_path <- "./output/text" # Path for text-output

#------------------------------------------------------------------------------------------
# Import data
#------------------------------------------------------------------------------------------

#import_data
dat <- as.data.frame(read_excel(paste0(data_path, "/", "dat.xlsx")))

# View(dat)
head(dat)
str(dat)

#------------------------------------------------------------------------------------------
# Categorise sampling time
#------------------------------------------------------------------------------------------

#remove wrong date from sampling time column
dat$"Sampling time" = gsub("1899-12-31","",dat$"Sampling time")

#Group sampling time into 5 categories
dat$Sampling_time_cat <- cut(
  chron::times(dat$"Sampling time")
  , breaks = chron::times(c(
    "08:00:00"
    , "10:00:00"
    , "12:00:00"
    , "14:00:00"
    , "16:00:00"
    , "23:59:59"
  ))
  , labels = c("8-9.59","10-11.59","12-13.59","14-15.59","16-18.59")
  , include.lowest = TRUE
  , right = TRUE
)

dat[, c("Sampling time", "Sampling_time_cat")]

#------------------------------------------------------------------------------------------
# Calculate Total Daily PA and add a column
#------------------------------------------------------------------------------------------

#Add column containing total Daily PA
dat$Total_PA <- dat$Total_LPA + dat$Total_MPA + dat$Total_VPA

#------------------------------------------------------------------------------------------
# Log2 transform data
#------------------------------------------------------------------------------------------

#log2_transformation_lipid_subclasses_and_clinical_lipids
(variables_to_transform <- names(dat)[15:43]) # Name of variables to transform with log2

for (i in variables_to_transform) {
  dat[, paste0(i, "_log2")] <- log2(dat[, i])
}

names(dat)

head(dat)

#------------------------------------------------------------------------------------------
# Descriptive scatterplots
#------------------------------------------------------------------------------------------

# Create scatterplot matrix to inspect correlation among clinical data
png(paste(graphics_path, "/", "scatterplot_matrix.png", sep = ""), width = 9*2.5, height = 9*2.5, units = "in", res = 300)

car::scatterplotMatrix(~Age + Sex + PBF + VO2peak_mlkgmin + Total_VPA + Statins + SLJ_PP_AVG_2_kN,  data = dat
                       , diagonal = list(method = "boxplot")
                       , smooth = list(method = loessLine, spread = FALSE, col.smooth = "red")
                       , by.groups = TRUE
                       , use = "pairwise.complete.obs"
                       , regLine  = FALSE
                       # , lwd = 2
                       , col = c("#00ABCE", "violet")
                       # , pch = c(15, 1)
                       # , cex = 1.5
)

dev.off()

#------------------------------------------------------------------------------------------
# Data standardisation
#------------------------------------------------------------------------------------------

#z_standardisation_dependent_variables
(variables_to_standardize <- names(dat)[c(5:6, 9:14, 45:74)])

for (i in variables_to_standardize) {
  dat[, paste0(i, "_std")] <- c(scale(dat[, i], center = TRUE, scale = TRUE))
}

names(dat)
head(dat)

#------------------------------------------------------------------------------------------
# Convert categorical data to factors
#------------------------------------------------------------------------------------------

dat$Sex <- factor(dat$Sex, levels = 0:1, labels = c("male", "female"))
dat$Statins <- factor(dat$Statins)
dat$Sampling_time_cat <- factor(dat$Sampling_time_cat)

#------------------------------------------------------------------------------------------
# Plot descriptive graphics
#------------------------------------------------------------------------------------------

# Some descriptive graphics
(vars_to_plot <- names(dat)[c(84:112)])

# Reshape to long
dat_long <- reshape2::melt(
  dat
  , id.vars = c("Sex", "Age", "Statins")
  , measure.vars = vars_to_plot
)

#theme_set(theme_bw())
p <- ggplot(data = dat_long, aes(x = Sex, y = value, group = Sex)) +
  geom_violin(aes(fill = Sex), trim = TRUE) +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 3) +
  # geom_point(size = 4, colour = "#00ABCE") +
  xlab("Sex") +
  ylab("Value") +
  # scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, 2)) +
  # scale_y_continuous(limits = c(0, 90), breaks = seq(0, 90, 20)) +
  scale_fill_manual(breaks = c(0, 1), values = c("#00ACEE", "#F4A749")) +
  facet_wrap(~variable, scales = "free_y") +
  theme(
    axis.title.y=element_text(colour = "black", size = 17, hjust = 0.5, margin=margin(0,12,0,0)),
    axis.title.x=element_text(colour = "black", size = 17),
    # axis.title.y=element_text(size=15,hjust=0.5, vjust=1),
    axis.text.x=element_text(colour = "black", size=15),
    axis.text.y=element_text(colour = "black", size=15),
    # plot.margin=unit(c(2,2,2,2,2),"line"),
    legend.position="right",
    legend.text=element_text(size=12.5),
    # panel.grid.minor = element_blank(),
    # panel.grid.major = element_line(colour=grey(0.8), size=0.5),
    legend.key=element_blank(),
    plot.title = element_text(face = "bold"),
    # legend.title=element_text(size=15),
    # legend.key.width=unit(.01,"npc"),
    # legend.key.height=unit(.025,"npc"),
    # strip.background=element_rect(fill="white")
    strip.text.x=element_text(size=15)
  )

p

ggsave(paste(graphics_path, paste("violin_sex.png", sep = ""), sep = "/"), p, width = 18*1, height = 14*1, units = "in", dpi = 300)

#------------------------------------------------------------------------------------------
# Count the missing values by column wise
#------------------------------------------------------------------------------------------

cat("\nCount of missing values by column wise\n")
sapply(dat, function(x) sum(is.na(x)))

#------------------------------------------------------------------------------------------
# Multiple linear regressions
#------------------------------------------------------------------------------------------

#univariate_multiple_linear_regressions
my_lms <- lapply(dat[,c(84:112)], function(x) lm(x ~ Age_std * Sex + VO2peak_mlkgmin_std + Total_PA_std + Statins + PBF_std + Sampling_time_cat + Fasting_std, data = dat))

#------------------------------------------------------------------------------------------
# Summaries & plot
#------------------------------------------------------------------------------------------

# summaries
lapply(my_lms, summary)

#------------------------------------------------------------------------------------------
# Check residuals
#------------------------------------------------------------------------------------------

# Get resid for all models
list_resid <- lapply(my_lms, resid)

# If you want them in a data.frame instead of list
df_resid <- do.call(cbind.data.frame, list_resid)

#plot residuals
p <- par(mfrow = c(2, 2))
lapply(my_lms, plot)
p

#------------------------------------------------------------------------------------------
# Check multicollinearity
#------------------------------------------------------------------------------------------

vif <- lapply(my_lms, function(model) car::vif(model, type = "predictor"))
vif

#------------------------------------------------------------------------------------------
# Extract variable from my_lms
#------------------------------------------------------------------------------------------

#extract Age results
res_frame_Age <- data.frame(
  term = names(my_lms)
  , estimate = NA
  , std.error = NA
  , statistic = NA
  , p.value = NA
  , response = "NA"
)

for (i in seq_along(names(my_lms))) {
  
  sum_tmp <- summary(my_lms[[i]])
  
  res_frame_Age$estimate[which(res_frame_Age$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Estimate"]["Age_std"]
  res_frame_Age$std.error[which(res_frame_Age$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Std. Error"]["Age_std"]
  res_frame_Age$statistic[which(res_frame_Age$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "t value"]["Age_std"]
  res_frame_Age$p.value[which(res_frame_Age$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Pr(>|t|)"]["Age_std"]
  res_frame_Age$response[which(res_frame_Age$term %in% names(my_lms)[i])] <- "Age (years)"
  
  rm(sum_tmp)
  
}

#extract Sex results
res_frame_Sex <- data.frame(
  term = names(my_lms)
  , estimate = NA
  , std.error = NA
  , statistic = NA
  , p.value = NA
  , response = "NA"
)

for (i in seq_along(names(my_lms))) {
  
  sum_tmp <- summary(my_lms[[i]])
  
  res_frame_Sex$estimate[which(res_frame_Sex$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Estimate"]["Sexfemale"]
  res_frame_Sex$std.error[which(res_frame_Sex$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Std. Error"]["Sexfemale"]
  res_frame_Sex$statistic[which(res_frame_Sex$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "t value"]["Sexfemale"]
  res_frame_Sex$p.value[which(res_frame_Sex$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Pr(>|t|)"]["Sexfemale"]
  res_frame_Sex$response[which(res_frame_Sex$term %in% names(my_lms)[i])] <- "Female"
  
  rm(sum_tmp)
  
}

#extract Age:Sex results
res_frame_Age_Sex <- data.frame(
  term = names(my_lms)
  , estimate = NA
  , std.error = NA
  , statistic = NA
  , p.value = NA
  , response = "NA"
)

for (i in seq_along(names(my_lms))) {
  
  sum_tmp <- summary(my_lms[[i]])
  
  res_frame_Age_Sex$estimate[which(res_frame_Age_Sex$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Estimate"]["Age_std:Sexfemale"]
  res_frame_Age_Sex$std.error[which(res_frame_Age_Sex$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Std. Error"]["Age_std:Sexfemale"]
  res_frame_Age_Sex$statistic[which(res_frame_Age_Sex$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "t value"]["Age_std:Sexfemale"]
  res_frame_Age_Sex$p.value[which(res_frame_Age_Sex$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Pr(>|t|)"]["Age_std:Sexfemale"]
  res_frame_Age_Sex$response[which(res_frame_Age_Sex$term %in% names(my_lms)[i])] <- "Interaction age:female"
  
  rm(sum_tmp)
  
}

#extract VO2peak results
res_frame_VO2peak <- data.frame(
  term = names(my_lms)
  , estimate = NA
  , std.error = NA
  , statistic = NA
  , p.value = NA
  , response = "NA"
)

for (i in seq_along(names(my_lms))) {
  
  sum_tmp <- summary(my_lms[[i]])
  
  res_frame_VO2peak$estimate[which(res_frame_VO2peak$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Estimate"]["VO2peak_mlkgmin_std"]
  res_frame_VO2peak$std.error[which(res_frame_VO2peak$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Std. Error"]["VO2peak_mlkgmin_std"]
  res_frame_VO2peak$statistic[which(res_frame_VO2peak$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "t value"]["VO2peak_mlkgmin_std"]
  res_frame_VO2peak$p.value[which(res_frame_VO2peak$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Pr(>|t|)"]["VO2peak_mlkgmin_std"]
  res_frame_VO2peak$response[which(res_frame_VO2peak$term %in% names(my_lms)[i])] <- "VO2peak (ml/min/kg)"
  
  rm(sum_tmp)
  
}

#extract Statins results
res_frame_Statins <- data.frame(
  term = names(my_lms)
  , estimate = NA
  , std.error = NA
  , statistic = NA
  , p.value = NA
  , response = "NA"
)

for (i in seq_along(names(my_lms))) {
  
  sum_tmp <- summary(my_lms[[i]])
  
  res_frame_Statins$estimate[which(res_frame_Statins$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Estimate"]["Statins1"]
  res_frame_Statins$std.error[which(res_frame_Statins$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Std. Error"]["Statins1"]
  res_frame_Statins$statistic[which(res_frame_Statins$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "t value"]["Statins1"]
  res_frame_Statins$p.value[which(res_frame_Statins$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Pr(>|t|)"]["Statins1"]
  res_frame_Statins$response[which(res_frame_Statins$term %in% names(my_lms)[i])] <- "Statin intake"
  
  rm(sum_tmp)
  
}

#extract PBF results
res_frame_PBF <- data.frame(
  term = names(my_lms)
  , estimate = NA
  , std.error = NA
  , statistic = NA
  , p.value = NA
  , response = "NA"
)

for (i in seq_along(names(my_lms))) {
  
  sum_tmp <- summary(my_lms[[i]])
  
  res_frame_PBF$estimate[which(res_frame_PBF$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Estimate"]["PBF_std"]
  res_frame_PBF$std.error[which(res_frame_PBF$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Std. Error"]["PBF_std"]
  res_frame_PBF$statistic[which(res_frame_PBF$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "t value"]["PBF_std"]
  res_frame_PBF$p.value[which(res_frame_PBF$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Pr(>|t|)"]["PBF_std"]
  res_frame_PBF$response[which(res_frame_PBF$term %in% names(my_lms)[i])] <- "Body fat (%)"
  
  rm(sum_tmp)
  
}

#extract Total PA results
res_frame_Total_PA <- data.frame(
  term = names(my_lms)
  , estimate = NA
  , std.error = NA
  , statistic = NA
  , p.value = NA
  , response = "NA"
)

for (i in seq_along(names(my_lms))) {
  
  sum_tmp <- summary(my_lms[[i]])
  
  res_frame_Total_PA$estimate[which(res_frame_Total_PA$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Estimate"]["Total_PA_std"]
  res_frame_Total_PA$std.error[which(res_frame_Total_PA$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Std. Error"]["Total_PA_std"]
  res_frame_Total_PA$statistic[which(res_frame_Total_PA$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "t value"]["Total_PA_std"]
  res_frame_Total_PA$p.value[which(res_frame_Total_PA$term %in% names(my_lms)[i])] <- sum_tmp$coefficients[, "Pr(>|t|)"]["Total_PA_std"]
  res_frame_Total_PA$response[which(res_frame_Total_PA$term %in% names(my_lms)[i])] <- "Total daily physical activity"
  
  rm(sum_tmp)
  
}

#------------------------------------------------------------------------------------------
# Extract and add confidence intervals to your result frames
#------------------------------------------------------------------------------------------

for (i in seq_along(names(my_lms))) {
  
  sum_tmp <- summary(my_lms[[i]])
  ci_tmp <- confint(my_lms[[i]], level = 0.95)  # Calculate confidence intervals

  # Extracting CI for Age
  lower_ci_age <- ci_tmp["Age_std", 1]  # Extract lower CI limit for Age_std
  upper_ci_age <- ci_tmp["Age_std", 2]  # Extract upper CI limit for Age_std
  res_frame_Age$lower_ci[which(res_frame_Age$term %in% names(my_lms)[i])] <- lower_ci_age
  res_frame_Age$upper_ci[which(res_frame_Age$term %in% names(my_lms)[i])] <- upper_ci_age

  # Extracting CI for Sex
  lower_ci_sex <- ci_tmp["Sex1", 1]
  upper_ci_sex <- ci_tmp["Sex1", 2]
  res_frame_Sex$lower_ci[which(res_frame_Sex$term %in% names(my_lms)[i])] <- lower_ci_sex
  res_frame_Sex$upper_ci[which(res_frame_Sex$term %in% names(my_lms)[i])] <- upper_ci_sex

  # Extracting CI for the interaction Age and Sex
  lower_ci_age_sex <- ci_tmp["Age_std:Sex1", 1]  
  upper_ci_age_sex <- ci_tmp["Age_std:Sex1", 2]  
  res_frame_Age_Sex$lower_ci[which(res_frame_Age_Sex$term %in% names(my_lms)[i])] <- lower_ci_age_sex
  res_frame_Age_Sex$upper_ci[which(res_frame_Age_Sex$term %in% names(my_lms)[i])] <- upper_ci_age_sex
  
  # Extracting CI for VO2peak
  lower_ci_VO2peak <- ci_tmp["VO2peak_mlkgmin_std", 1]  
  upper_ci_VO2peak <- ci_tmp["VO2peak_mlkgmin_std", 2]  
  res_frame_VO2peak$lower_ci[which(res_frame_VO2peak$term %in% names(my_lms)[i])] <- lower_ci_VO2peak
  res_frame_VO2peak$upper_ci[which(res_frame_VO2peak$term %in% names(my_lms)[i])] <- upper_ci_VO2peak
  
  # Extracting CI for Statins
  lower_ci_statins <- ci_tmp["Statins1", 1]  
  upper_ci_statins <- ci_tmp["Statins1", 2]  
  res_frame_Statins$lower_ci[which(res_frame_Statins$term %in% names(my_lms)[i])] <- lower_ci_statins
  res_frame_Statins$upper_ci[which(res_frame_Statins$term %in% names(my_lms)[i])] <- upper_ci_statins
  
  # Extracting CI for PBF
  lower_ci_pbf <- ci_tmp["PBF_std", 1]  
  upper_ci_pbf <- ci_tmp["PBF_std", 2]  
  res_frame_PBF$lower_ci[which(res_frame_PBF$term %in% names(my_lms)[i])] <- lower_ci_pbf
  res_frame_PBF$upper_ci[which(res_frame_PBF$term %in% names(my_lms)[i])] <- upper_ci_pbf
  
  # Extracting CI for Total PA
  lower_ci_pa <- ci_tmp["Total_PA_std", 1]  
  upper_ci_pa <- ci_tmp["Total_PA_std", 2]  
  res_frame_Total_PA$lower_ci[which(res_frame_Total_PA$term %in% names(my_lms)[i])] <- lower_ci_pa
  res_frame_Total_PA$upper_ci[which(res_frame_Total_PA$term %in% names(my_lms)[i])] <- upper_ci_pa
  
  # Clear temporary variables
  rm(sum_tmp, ci_tmp, lower_ci_age, upper_ci_age, lower_ci_sex, upper_ci_sex, lower_ci_age_sex, upper_ci_age_sex,
    lower_ci_VO2peak, upper_ci_VO2peak, lower_ci_statins, upper_ci_statins, lower_ci_pbf, upper_ci_pbf, lower_ci_pa, upper_ci_pa)
}

#------------------------------------------------------------------------------------------
# Combine extracted data into a data frame
#------------------------------------------------------------------------------------------

#Combine_data_frames
Overall <- rbind(
  res_frame_Age
  , res_frame_Sex
  , res_frame_Age_Sex
  , res_frame_VO2peak
  , res_frame_Statins
  , res_frame_Total_PA
  , res_frame_PBF
  )

#------------------------------------------------------------------------------------------
# Adjust p-values for multiple testing
#------------------------------------------------------------------------------------------

#Adjust p-values
Overall <- Overall[order(Overall$p.value), ]
Overall$BH <- p.adjust(Overall$p.value, method = "BH")

#Categorise BH p-values
Overall$BH_cat <- NA

Overall$BH_cat[Overall$BH > 0.05] <- "> 0.05"
Overall$BH_cat[Overall$BH <= 0.05 & Overall$BH > 0.01] <- "≤ 0.05"
Overall$BH_cat[Overall$BH <= 0.01 & Overall$BH > 0.001] <- "≤ 0.01"
Overall$BH_cat[Overall$BH <= 0.001 & Overall$BH > 0.0001] <- "≤ 0.001"
Overall$BH_cat[Overall$BH <= 0.0001] <- "≤ 0.0001"

Overall$BH_cat <- factor(Overall$BH_cat, levels = c("> 0.05", "≤ 0.05", "≤ 0.01", "≤ 0.001", "≤ 0.0001"))

#Order by decreasing estimates
Overall <- Overall[order(-Overall$estimate), ]

#remove suffixes from lipid subclasses' name
Overall$term=gsub("_log2_std","",Overall$term)

Overall_final <- Overall

#rename columns
names(Overall_final)[1] <- "Dependent variable"
names(Overall_final)[2] <- "β coefficient"
names(Overall_final)[3] <- "standard error"
names(Overall_final)[5] <- "p-value"
names(Overall_final)[6] <- "Independent variables"
names(Overall_final)[7] <- "95% CI lower bound for β coefficient"
names(Overall_final)[8] <- "95% CI higher bound for β coefficient"
names(Overall_final)[9] <- "BH p-value"
names(Overall_final)[10] <- "Categorical BH p-value"

#reorder columns
Overall_final <- Overall_final[, c(1, 6, 2, 7, 8, 3, 4, 5, 9, 10)]
head(Overall_final)

# write_xlsx(Overall_final, "results.xlsx")

#------------------------------------------------------------------------------------------
# Start rain plot
#------------------------------------------------------------------------------------------

#Import data
plot_data <- Overall

# Define Theme and Palette

## Palette

palette <-
  # Blue
  c("#053061",
             "#313695",
             "#4575b4",
             "#74add1",
             "#abd9e9",
             "#e0f3f8",
             "#fee090",
             "#fdae61",
             "#f46d43",
             "#d73027",
             "#a50026",
             '#67001f')
             # Red

# Calculate symmetric limits based on most extreme value
max_abs_estimate <- max(abs(plot_data$estimate))

max_lim <- max_abs_estimate
min_lim = -1 * max_lim

## theme

thm <-
  # Good starting theme + set text size
  theme_light(base_size = 7) +
  theme(
    # Remove axis ticks and titles
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    # Remove gridlines and boxes
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    legend.key = element_blank(),
    
    # White backgrounds
    panel.background = element_rect(fill = 'white'),
    plot.background = element_rect(fill = 'white'),
    legend.background = element_rect(fill = 'white'),
    
    # Angle text
    axis.text.x.top = element_text(angle = 45, hjust = 0)
  )

#------------------------------------------------------------------------------------------
# Bare-bones rainplot
#------------------------------------------------------------------------------------------

plot_data$response <- factor(plot_data$response, levels = c("Age (years)", "Female", "Interaction age:female", "Total daily physical activity", "Body fat (%)", "Statin intake", "VO2peak (ml/min/kg)"))

rainplot <-
  ggplot(plot_data) +
  geom_point(aes(x = response, y = term, colour = estimate, size = BH_cat))

print(rainplot)

#------------------------------------------------------------------------------------------
# Basic Rainplot
#------------------------------------------------------------------------------------------

rainplot <-
  ggplot(plot_data) +
  geom_point(aes(x = response, y = term, colour = estimate, size = BH_cat))  +
  scale_x_discrete(position = 'top') +
  scale_size_manual(name = expression("BH p-value"), values = c(2, 4, 6, 8, 10), drop = FALSE) +
  scale_color_gradientn(
    'Effect Size\n(β coefficient)',
    colors = palette,
    limits = c(min_lim, max_lim),
    breaks = c(min_lim, min_lim / 2, 0 , max_lim / 2, max_lim),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  guides(colour = guide_colourbar(order = 1),
         size = guide_legend(order = 2)) +
  thm

print(rainplot)

#------------------------------------------------------------------------------------------
# Ordering by Cluster
#------------------------------------------------------------------------------------------

# Convert to matrix and reshape for clustering.
cluster_data <-
  plot_data %>%
  select(response, term, estimate) %>%
  spread(response, estimate)

rnms <-
  cluster_data$term

cluster_data <-
  cluster_data %>%
  select(-term) %>%
  as.matrix()

rownames(cluster_data) <- rnms

# cluster dependent variable terms
clust <- hclust(dist(cluster_data), method = 'ward.D2')

# `clust$order` orders `term` into clusters
term_order <-
  clust$labels[clust$order]

# Convert term to a factor, ordered by `term_order`
plot_data_clo <-
  plot_data %>%
  mutate(term = factor(term, levels = term_order))

rainplot <-
  # Use cluter ordered data
  ggplot(plot_data_clo) +
  geom_point(aes(x = response, y = term, colour = estimate, size = BH_cat)) +
  scale_x_discrete(position = 'top') +
  scale_size_manual(name = expression("BH p-value"), values = c(2, 4, 6, 8, 10), drop = FALSE) +
  scale_color_gradientn(
    'Effect Size\n(β coefficient)',
    colors = palette,
    limits = c(min_lim, max_lim),
    breaks = c(min_lim, min_lim / 2, 0 , max_lim / 2, max_lim),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  guides(colour = guide_colourbar(order = 1),
         size = guide_legend(order = 2)) +
  thm

print(rainplot)

#------------------------------------------------------------------------------------------
# Adding dendrograms
#------------------------------------------------------------------------------------------

dendro_dat <- segment(dendro_data(clust))


dendro <-
  # Empty ggplot with same y-scale as rainplot
  ggplot() +
  geom_blank(aes(y = term), data = plot_data) +
  theme_dendro() +
  # 'expand' controls whitespace around the dendrogram. The non-zero argument
  # may need to be increasesed if the line thickness of the dendrogram is
  # increased to make sure the entire dendrogram is plotted
  scale_x_discrete(position = 'top', expand = c(0, 0.03, 0, 0)) +
  # Draw dendrogram
  geom_segment(aes(x = -y, y = x, xend = -yend, yend = xend),
               colour = 'black',
               data = dendro_dat)


p <- ggarrange(dendro, rainplot, ncol = 2, widths = c(1, 5))

p

ggsave(paste(graphics_path, paste("rain_plot.png", sep = ""), sep = "/"), p, width = 6*1, height = 10*1, units = "in", dpi = 300)
