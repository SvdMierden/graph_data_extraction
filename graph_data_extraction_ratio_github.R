################################################################################
## R script for comparing data extracted from graphs with the orginal data
## The script needs the 6 excel files with the extracted data (graph_extraction_X.xlsx)
## as well as the excel sheet with the orginal data (Orignal Data.xlsx)
################################################################################


rm(list = ls())

library(tidyverse)
library(readxl)
library(naniar)
library(epiR)
library(metafor)
library(viridis)
library(ggthemes)


#set theme
theme_set(theme_minimal() +
          theme(axis.text = element_text(size = 12),
                axis.title = element_text(size = 12)))


##
# Import the data
##

# Set working directory to where the filled in extraction sheets
# setwd()

# select and import data extraction files
excel_files <- list.files(pattern = ".xlsx")
files <- lapply(excel_files, read_xlsx, col_types = c(
    rep("guess", 4),
    rep("skip", 2),
    rep("guess", 2),
    rep("skip", 2),
    "guess",
    "skip",
    "guess",
    rep("skip", 2),
    "guess",
    "skip",
    "guess",
    rep("skip", 4)
    )
  )
data <- bind_rows(files, .id = "extractor")


# fill the missing data
# fill missing data for $Figure for the Gurfein studies (special cases)
data <- data %>%
    mutate(Group = case_when(
      Figure == "S 1A" ~ "Cntl",
      Figure == "S 1B" ~ "Calm",
      TRUE ~ Group)
    )

# fill explicitly for $Author, $title, and $Figure
data <- data %>%
  fill(Author, Title, Figure)

# Create unique ID per graph data point, this can be combined with $extractor for unqiue ID per data point, and arrange data
data <- data %>%
  unite(col = data_point, Figure, Group, sep = "_", remove = TRUE) %>%
  group_by(extractor) %>%
  arrange(data_point, .by_group = TRUE) %>%
  select(-Title, study = Author) %>%
  rename(ex_conc_val = "Concentration Value",
         ex_err_val = "Error Value")


# Import original values
# point the pathway to the excel file with the orginal data
data_orig <- read_xlsx("Original Data.xlsx") %>%
  select(Author:Error) %>%
  rename(orig_conc_val = Concentration,
         orig_err_val = Error)

# fill explicitly for $Author, $title, and $Figure
data_orig <- data_orig %>%
  fill(Author, Title, Figure, Unit)

# Create unique ID per graph data point and arrange data
data_orig <- data_orig %>%
  unite(col = data_point, Figure, Group, sep = "_", remove = TRUE) %>%
  arrange(data_point) %>%
  select(-Title, study = Author)



##
# calculate CCCs
##

# Create data_wide which can be used for CCC
data_wide <- data %>%
  pivot_wider(id_cols = c(data_point, study), names_from = extractor, values_from = c(ex_conc_val, ex_err_val), names_sep = "_") %>%
  bind_cols(data_orig[,c(-1)]) %>%
  select(-data_point...15) %>%
  arrange(study) %>%
  rename(data_point = data_point...1)


# Add numerical Study identifiers to data_wide for forest plots
unique_id_wide <- vector()
for (i in 1:length(unique(data_wide$study))) {
  temp <- rep(i, times = length(which(data_wide$study == unique(data_wide$study)[i])))
  unique_id_wide <- c(unique_id_wide, temp)
}

data_wide <- cbind(data_wide, unique_id_wide)

  
# calculate concordance correlation coefficients for extractors compared to original values, for both concentrations and errors
# epi.ccc for some reason cannot use tibbles, so first we need to convert data_wide back to a data.frame
data_wide <- as.data.frame(data_wide)

# Calculate CCC for concentration
CCC_conc_li <- list()
for(i in 1:6){
  CCC_conc_li[[i]] <- epi.ccc(data_wide[, 2 + i], data_wide$orig_conc_val)$rho.c
}
ccc_conc <- bind_rows(CCC_conc_li, .id = "extractor")

# Calculate ccc for error 
CCC_err_li <- list()
for(i in 1:6){
  CCC_err_li[[i]] <- epi.ccc(data_wide[, grep("ex_err_val_1", colnames(data_wide)) -1 + i], data_wide$orig_err_val)$rho.c
}
ccc_err <- bind_rows(CCC_err_li, .id = "extractor")


#ccc plots
#create indiviudal plots
#ccc plot for concentration
(ccc_conc_plot <- ggplot(ccc_conc, aes(x = extractor, y = est)) +
    geom_hline(yintercept = mean(ccc_conc$est), colour = "blue") +
    geom_point(size = 3) +
    geom_linerange(aes(ymin = lower, ymax = upper), size = 1.3) +
    ylab("ccc outcome values") +
    xlab(NULL) +
    theme(axis.text.x = element_blank()) +
    scale_y_continuous(breaks = seq(from = 0.995, to = 1, length.out = 6),
                       minor_breaks = seq(from = 0.995, to = 1, length.out = 6))
)

ggsave("ccc_conc_plot.tiff",
       dpi = 600)


# ccc plot for error
(ccc_err_plot <- ggplot(ccc_err, aes(x = extractor, y = est)) +
    geom_hline(yintercept = mean(ccc_err$est), colour = "blue") +
    geom_point(size = 3) +
    geom_linerange(aes(ymin = lower, ymax = upper), size = 1.3) +
    ylab("ccc error") +
    scale_y_continuous(breaks = seq(from = 0.850, to = 1, length.out = 6),
                       minor_breaks =  seq(from = 0.850, to = 1, length.out = 6))
)

ggsave("ccc_err_plot.tiff",
       dpi = 600)


#combine conc plot & error plot into one plot
(ccc_both_plot <- ggarrange(ccc_conc_plot, ccc_err_plot, ncol = 1, nrow = 2, common.legend = TRUE, legend = "bottom")
)
ggsave("ccc_comb_plot.tiff",
       dpi = 600)


##
# Calculate overall ccc (occc) between extractors only
##
# occc for concentration
occc_conc <- data_wide %>%
  select(ex_conc_val_1:ex_conc_val_6) %>%
  epi.occc(pairs = TRUE)

# occc for errors
occc_err <- data_wide %>%
  select(ex_err_val_1:ex_err_val_6) %>%
  epi.occc(pairs = TRUE)



##
# bland altmann plot
##
## Standardise mean: Divide by SEM
extractor_conc_stand <- data_wide[, grep("ex_conc|orig_conc", colnames(data_wide))]/data_wide[, grep("ex_err|orig_err", colnames(data_wide))]

## Standardise SEM: method 1 (see, link)
min_extractor <- apply(data_wide[, grep("ex_err|orig_err", colnames(data_wide))], 2, min)
max_extractor <- apply(data_wide[, grep("ex_err|orig_err", colnames(data_wide))], 2, max)

extractor_SEM_stand <- data_wide[, grep("ex_err|orig_err", colnames(data_wide))]

# Add 0.001 to the standardised error because the standardising introduces 0 by defintion, and this gives error later on when calculating ratios between errors.
for(i in 1:length(data_wide[, 1])){
  extractor_SEM_stand[i, ] <- (data_wide[i, grep("ex_err|orig_err", colnames(data_wide))] - min_extractor)/(max_extractor - min_extractor)  + 0.001
}


# Create long formats for plotting
data_ba_conc_stand <- extractor_conc_stand %>%
  cbind(data_wide$study) %>%
  rename(study = "data_wide$study") %>%
  pivot_longer(-c(orig_conc_val, study), names_to = "extractor", names_prefix = "ex_conc_val_",values_to = "stand_conc") %>%
  mutate(mean_conc = (orig_conc_val + stand_conc) / 2,
         log_ratio_conc = log(orig_conc_val / stand_conc),
         ratio_conc = orig_conc_val / stand_conc)

data_ba_sem_stand <- extractor_SEM_stand %>%
  pivot_longer(-orig_err_val, names_to = "extractor", names_prefix = "ex_err_val_",values_to = "stand_err") %>%
  mutate(mean_err = (orig_err_val + stand_err) / 2,
         log_ratio_err = log(orig_err_val / stand_err),
         ratio_err = orig_err_val /stand_err
         )


# Create Auxilarry data frame for hlines in plot
data_ba <- bind_cols(data_ba_conc_stand, data_ba_sem_stand) %>%
  select(-extractor...9) %>%
  rename(extractor = extractor...3)

# Create numeric identiiers for the study variable.
unique_id <- vector()
for (i in 1:length(unique(data_ba$study))) {
  temp <- rep(i, times = length(which(data_ba$study == unique(data_ba$study)[i])))
  unique_id <- c(unique_id, temp)
}

data_ba <- cbind(data_ba, as.character(unique_id))

# Create axuilary data table for the BA graphs
data_ba_aux <- data_ba %>%
  summarise(mean_diff_conc = mean(ratio_conc), 
            ci_up_con = exp(mean(log_ratio_conc) + (1.96 * sd(log_ratio_conc))),
            ci_do_con = exp(mean(log_ratio_conc) - (1.96 * sd(log_ratio_conc))),
            mean_diff_err = mean(ratio_err),
            ci_up_err = exp(mean(log_ratio_err) + (1.96 * sd(log_ratio_err))),
            ci_do_err = exp(mean(log_ratio_err) - (1.96 * sd(log_ratio_err))))

# Remove the largest outlier from data_ba_sem_stand to see the influence on the bias
data_ba_sem_outl <- data_ba_sem_stand %>%
  filter(ratio_err < 15) %>%
  group_by(extractor) %>%
  summarise(mean(ratio_err))

# Plot
# ba plot for concentration
# Change y = ratio_conc for y = diff_conc to plot difference on y-axis
(ba_conc <- ggplot(data_ba, aes(x = mean_conc, y = ratio_conc, colour = as.character(unique_id))) +
  geom_point(size = 3, alpha = 0.65) +
  scale_color_viridis(discrete = TRUE, 
                      end = 0.90, 
                      option = "A", 
                      limits = 1:22,
                      name = "Study") +
  geom_hline(yintercept = 1, 
             linetype = "dashed") +
  geom_hline(data = data_ba_aux, aes(yintercept = mean_diff_conc),
             colour = "blue",
             size = 0.7) +
  geom_hline(data = data_ba_aux, aes(yintercept = ci_up_con),
             colour = "red",
             size = 0.5) +
  geom_hline(data = data_ba_aux, aes(yintercept = ci_do_con),
             colour = "red",
             size = 0.5) +
  # scale_y_log10() +
  facet_wrap(vars(extractor)) +
  labs(x = "Means of orignal & extracted standardised outcome values",
       y = "Ratio of original to extracted \nstandardised outcome values") +
  theme(panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0.9, "lines")
    )
)

ggsave("ba_conc.tiff",
       dpi = 600)

# ba plot for error
# change y = ratio_err for y = diff_err to plot difference in stead of ration on y-axis.
# nb, hlines cannot be plotted when data_ba_aux is created using the error ratio in stead of error difference.

(ba_err <- ggplot(data_ba, aes(x = mean_err, y = ratio_err, colour = as.character(unique_id))) +
    geom_point(size = 3, alpha = 0.65) +
    scale_color_viridis(discrete = TRUE, 
                        end = 0.90, 
                        option = "A", 
                        limits = 1:22,
                        name = "Study") +
    geom_hline(yintercept = 1, 
               linetype = "dashed") +
    geom_hline(data = data_ba_aux, aes(yintercept = mean_diff_err),
               colour = "blue",
               size = 0.7) +
    geom_hline(data = data_ba_aux, aes(yintercept = ci_up_err),
               colour = "red",
               size = 0.5) +
    geom_hline(data = data_ba_aux, aes(yintercept = ci_do_err),
               colour = "red",
               size = 0.5) +
    facet_wrap(vars(extractor)) +
    scale_y_log10() +
    labs(x = "Means of orignal & extracted standardised error",
         y = "Ratio of original to extracted \nstandardised error") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.spacing.x = unit(0.9, "lines")
          )
)

ggsave("ba_err.tiff",
       dpi = 600)

    
# Combining the two ba plots into one figure.
(save_plots <- ggarrange(ba_conc, ba_err, 
                         ncol = 1, 
                         nrow = 2, 
                         common.legend = TRUE, 
                         legend = "right"
                         )
  )

ggsave("ba_both.tiff",
       dpi = 600)




### START CODE FOR RANDOM-EFFECT META-ANALYSIS
## Apply RE-MA with REML estimator of between-trial variance
corr <- seq(0.0, 0.99, 0.05)
trial.MD <- trial.var <- est.RE <- std.err.RE <- est.FE <- std.err.FE <- AIC <- list()

for(i in 1:length(corr)){
  
  trial.var[[i]] <- est.RE[[i]]  <- std.err.RE[[i]] <- est.FE[[i]]  <- std.err.FE[[i]] <- AIC[[i]] <- rep(list(), 6)
  
  for(j in 1:6){
    # Calculate within trial mean difference (between original and extrcted data) for each extractor
    trial.MD[[j]] <- data_wide[, 15] - data_wide[, j + 2]
    
    # Calculate within trial variance for each extractor assuming positive correlation (mimic split-mouth study)
    trial.var[[i]][[j]] <- data_wide[, 17]^2 + data_wide[, 8 + j]^2 - 2*corr[i]*data_wide[, 17]*data_wide[, 8 + j]
    
    # Obtain the pooled MD per extractor - RANDOM-EFFECTS MA
    est.RE[[i]][[j]] <- rma(yi = trial.MD[[j]], vi = trial.var[[i]][[j]], method = "REML")$beta[, 1]
    
    # Obtain the pooled MD per extractor - FIXED-EFFECT MA
    est.FE[[i]][[j]] <- rma(yi = trial.MD[[j]], vi = trial.var[[i]][[j]], method = "FE")$beta[, 1]
    
    # Obtain the standard error of the pooled MD per extractor (RE-MA)
    std.err.RE[[i]][[j]] <- rma(yi = trial.MD[[j]], vi = trial.var[[i]][[j]], method = "REML")$se
    
    # Obtain the standard error of the pooled MD per extractor (FE-MA)
    std.err.FE[[i]][[j]] <- rma(yi = trial.MD[[j]], vi = trial.var[[i]][[j]], method = "FE")$se
    
    AIC[[i]][[j]] <- fitstats(rma(yi = trial.MD[[j]], vi = trial.var[[i]][[j]], method = "REML"))[3, ]
  }
}
### END CODE FOR RANDOM-EFFECTs & FIXED-EFFECT META-ANALYSIS



### START CODE FOR AIC: 'BEST' CORRELATION COEFFICIENT
## A matrix of AICs: extractor in column and correlation in row
table.aic <- matrix(unlist(AIC), nrow = length(corr), ncol = 6, byrow = T)
rownames(table.aic) <- c(as.character(corr))



## Find correlation with the smallest AIC (expected, actually)
which(table.aic == apply(table.aic, 2, min), arr.ind = TRUE)
### END CODE FOR AIC: 'BEST' CORRELATION COEFFICIENT



## Create a dataframe to visualise
(dataset <- data.frame(unlist(est.RE[[20]]), unlist(std.err.RE[[20]]), 1:6))
colnames(dataset) <- c("MD", "SE", "extractor"); dataset

fig <- ggplot(data = dataset, aes(x = as.factor(extractor), y = MD, ymin = MD - 1.93*SE, ymax = MD + 1.93*SE) ) +  
  geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, lty = 2, col = "grey") +
  geom_point(size = 1.5,  colour = "white", stroke = 0.3, position = position_dodge(width = 0.5)) + 
  geom_text(aes(x = extractor, y = round(MD, 2), label = round(MD, 2)), color = "blue", hjust = -0.3, vjust = -0.4, size = 4, 
            check_overlap = F, parse = F, position = position_dodge(width = 0.8), inherit.aes = T) +
  labs(x = "Extractor", y = "Mean difference") +
  coord_flip() +
  ggtitle("Mean difference between original and extracted data") +
  theme_classic() + 
  theme(axis.title.y = element_text(color = "black", size = 14), axis.title.x = element_text(color = "black", size = 14),
        axis.text.x = element_text(color = "black", size = 14), axis.text.y = element_text(color = "black", size = 14),
        plot.title = element_text(size = 14, hjust = 0.5))  
fig

ggsave("MD_orig_extrac.tiff",
       dpi = 600)


### Between-trial variance is estimated to be zero in all extractors and the forestplot ensures that! Check the code below:
for (i in 1:6) {
  here <- rma(yi = trial.MD[[i]], vi = trial.var[[20]][[i]], method = "SJ")
  confint(here)
  
  tiff(paste("./forest", i, ".tiff", sep = ""), height = 18, width = 16, units = "cm", compression = "lzw", res = 600)
  forest(here, 
         order = "prec",
         header= c("Study", "Mean difference [95% CI]"),
         top = 1,
         slab = data_wide$study)
  dev.off()
}


## Create a dataframe to visualise the extractor-specific pooled MD per correlation coefficient
(dataset2 <- data.frame(unlist(est.RE), unlist(std.err.RE), rep(paste("Extractor",1:6), 20), rep(corr, each = 6)))
colnames(dataset2) <- c("MD", "se.MD", "Extractor", "rho")
quantile(unlist(est.RE), c(0.25, 0.50, 0.75))

tiff("./Extractor-specific MD per pho.tiff", height = 20, width = 25, units = "cm", compression = "lzw", res = 600)
ggplot(data = dataset2, aes(x = as.factor(rho), y = MD, ymin = MD - 1.93*se.MD, ymax = MD + 1.93*se.MD)) +  
  geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, lty = 2, col = "grey") +
  geom_point(size = 1.8,  colour = "white", stroke = 0.3, position = position_dodge(width = 0.5)) + 
  labs(x = "Correlation coefficient", y = "Mean difference") +
  facet_wrap(Extractor ~.)+
  coord_flip() +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12, face = "bold"), strip.text = element_text(size = 12, face = "bold"))
dev.off()


## Create a dataframe to visualise the FE- vs RE-MA for each extractor and a randomly selected rho (0.95)
(dataset3 <- data.frame(c(unlist(est.RE[[1]]), unlist(est.RE[[11]]), unlist(est.RE[[20]]), 
                          unlist(est.FE[[1]]), unlist(est.FE[[11]]), unlist(est.FE[[20]])), 
                        c(unlist(std.err.RE[[1]]), unlist(std.err.RE[[11]]), unlist(std.err.RE[[20]]), 
                          unlist(std.err.FE[[1]]), unlist(std.err.FE[[11]]), unlist(std.err.FE[[20]])), rep(1:6, 3*2), 
                        rep(c("Random-effects", "Fixed-effect"), each = 6*3), rep(rep(c(paste("Correlation coefficient =", 0), paste("Correlation coefficient =", 0.50), paste("Correlation coefficient =", 0.95)), each = 6), 2) ))  
colnames(dataset3) <- c("MD", "se.MD", "Extractor", "Model", "Rho")

tiff("./RE vs FE per extractor.tiff", height = 20, width = 25, units = "cm", compression = "lzw", res = 600)
ggplot(data = dataset3, aes(x = as.factor(Extractor), y = MD, ymin = MD - 1.93*se.MD, ymax = MD + 1.93*se.MD, colour = Model, group = Model)) +  
  geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
  geom_point(colour = "white", size = 1.8, stroke = 0.3, position = position_dodge(width = 0.5)) + 
  geom_hline(yintercept = 0, lty = 2, col = "grey") +
  labs(x = "Extractor", y = "Mean difference") +
  facet_wrap(Rho ~.)+
  coord_flip() +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom",
        legend.text = element_text(size = 12), legend.title = element_text(size = 12, face = "bold"), strip.text = element_text(size = 12, face = "bold"))
dev.off()

##end##