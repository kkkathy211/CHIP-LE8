############################################################
############ clean version - used for manuscript############
############################################################

# 0) Load required libraries
library(data.table)
library(dplyr)
library(survival)
library(broom)
library(tidyr)
library(ggplot2)
library(stringr)

setwd("/medpop/esp2")

# =======================================
# 1) Data Input & Preparation
# =======================================

# Load LE8 composite score data
le8_data <- fread("/medpop/esp2/yxliu/LE8score_calculation/updated_LE8score.txt")

# Load CHIP calls
chip_data <- fread("/medpop/esp2/tnakao/data/UKBB/phenotype/CHIP/450k/2022Oct/CHIP_calls_Oct16_2022.txt")

# Load the "mask" files for filtering
id_450k <- fread("/medpop/esp2/mesbah/datasets/CHIP/UKBB/450k/ukb450k.eid.list") %>% rename(id = V1)
no_consent <- fread("/medpop/esp2/projects/UK_Biobank/withdrawn_samples/w7089_20241217.csv") %>% rename(id = V1)
relate_remove <- fread("/medpop/esp2/zyu/chip_protemoics/listofremovefor450K.txt") %>% rename(id = IID)

# This file also needs to contain the baseline variables (age, sex, PCs, etc.)
pheno_master <- fread("/medpop/esp2/zyu/UKB_bigquery_2023pheno/20250711_pheno_noGP_clean.tsv.gz") %>%
  as.data.frame()

ori_pheno <- fread("/medpop/esp2/mzekavat/UKBB/ukbb_PhenoFile.ALL_500k.UpdatedIncdPhenos_202020.REAL.txt") %>%
  as.data.frame() %>%
  filter(Submitted_Gender == Inferred_Gender,
         Non_Consented     == 0,
         Prev_Myeloid_leukemia  == 0,
         Prev_Lymphoid_leukemia == 0
  ) %>%
  select(
    id, age, Sex_numeric, in_white_British_ancestry_subset,
    paste0("PC",1:10)
  )

pheno_master2 <- merge(
  pheno_master, ori_pheno,
  by = "id",
  all = TRUE
)

# Subset pheno_master to just the 450K participants with extensive QC
# and rename outcome columns to match existing code's expected pattern
pheno_filtered_renamed <- pheno_master2 %>%
  # keep only those in the 450K ID list
  filter(id %in% id_450k$id) %>%
  # remove withdrawn consent
  filter(!id %in% no_consent$id) %>%
  # only unrelated individuals, drop related ones
  filter(!id %in% relate_remove$id) %>%
  # apply existing QC steps (assuming these fields exist in pheno_master)
  # Select baseline columns and all relevant outcome columns
  select(
    id,
    age, 
    Sex_numeric, 
    in_white_British_ancestry_subset, 
    paste0("PC", 1:10), 
    # Select and rename prevalence flags
    Prev_composite_mi_cad_stroke_death = prevalent_disease_Coronary_Artery_Disease_SOFT,
    Prev_Coronary_Artery_Disease_INTERMEDIATE = prevalent_disease_Coronary_Artery_Disease_INTERMEDIATE,
    Prev_Coronary_Artery_Disease_Death = prevalent_disease_CAD_Death,
    # Select and rename follow-up times
    FollowUp_composite_mi_cad_stroke_death = followup_Coronary_Artery_Disease_SOFT,
    FollowUp_Coronary_Artery_Disease_INTERMEDIATE = followup_Coronary_Artery_Disease_INTERMEDIATE,
    FollowUp_Coronary_Artery_Disease_Death = followup_CAD_Death,
    # Select and rename event indicators
    Incd_composite_mi_cad_stroke_death = incident_disease_Coronary_Artery_Disease_SOFT,
    Incd_Coronary_Artery_Disease_INTERMEDIATE = incident_disease_Coronary_Artery_Disease_INTERMEDIATE,
    Incd_Coronary_Artery_Disease_Death = incident_disease_CAD_Death
  )


# Merge LE8 and the fully processed phenotype data (inner-join so only complete cases remain)
merged_data <- merge(
  le8_data, pheno_filtered_renamed,
  by = "id",
  all = FALSE
)


# --- CHECK SAMPLE SIZE AND MEDIAN FOLLOW-UP  ---

# 1. Exclude prevalent CAD 
cad0 <- merged_data %>%
  filter(Prev_Coronary_Artery_Disease_INTERMEDIATE == 0)

# 2. Compute sample size
n_no_cad <- nrow(cad0)
print(n_no_cad)
# [1] 410969


# 3. Compute median follow-up and IQR 
followup_times_int <- cad0$FollowUp_Coronary_Artery_Disease_INTERMEDIATE

median_fu_int <- median(followup_times_int, na.rm = TRUE)
iqr_fu_int <- quantile(followup_times_int, probs = c(0.25, 0.75), na.rm = TRUE)

cat("Median follow-up (Int CAD):", round(median_fu_int, 1),
    "years (IQR ",
    round(iqr_fu_int[1], 1), "–", round(iqr_fu_int[2], 1), " years)\n")

# Median follow-up (Soft CAD): 13.6 years (IQR  12.7 – 14.3  years)




# =======================================
# 2) Process CHIP Data (5 subtypes)
# =======================================
# now: AF >= 0.02 (2%)
CHIP_ids <- unique(chip_data[AF >= 0.02,]$Broad_ID)
ASXL1_ids <- unique(chip_data[Gene.refGene == "ASXL1" & AF >= 0.02,]$Broad_ID)
TET2_ids <- unique(chip_data[Gene.refGene == "TET2" & AF >= 0.02,]$Broad_ID)
DNMT3A_ids <- unique(chip_data[Gene.refGene == "DNMT3A" & AF >= 0.02,]$Broad_ID)
nonDNMT3A_ids <- unique(chip_data[AF >= 0.02 & Gene.refGene != "DNMT3A",]$Broad_ID)

merged_data <- merged_data %>%
  mutate(
    CHIP = as.integer(id %in% CHIP_ids),
    ASXL1 = as.integer(id %in% ASXL1_ids),
    TET2 = as.integer(id %in% TET2_ids),
    DNMT3A = as.integer(id %in% DNMT3A_ids),
    non_DNMT3A = as.integer(id %in% nonDNMT3A_ids)
  )

# =======================================
# 3) Process LE8 Composite Exposure
# =======================================
merged_data <- merged_data %>%
  mutate(
    LE8_composite = (diet_score + BP_Score + PA_score + final_smoking_score +
                       sleep_health_points + bmi_points + Blood_lipid_points + HbA1c_score) / 8,
    LE8_composite_scaled = LE8_composite / 10,
    LE8_cat = case_when(
      LE8_composite < 50 ~ "<50",
      LE8_composite < 80 ~ "50–<80",
      TRUE ~ "≥80"
    )
  ) %>%
  mutate(LE8_cat = factor(LE8_cat, levels = c("<50", "50–<80", "≥80")))

# =======================================
# 4) Stratified Cox Regression
# =======================================
outcome_list <- c(
  "composite_mi_cad_stroke_death",
  "Coronary_Artery_Disease_INTERMEDIATE",
  "Coronary_Artery_Disease_Death"
)
chip_var <- c("CHIP", "ASXL1", "TET2", "DNMT3A", "non_DNMT3A")
le8_variables <- c("LE8_composite_scaled")

results_strat <- data.frame(
  exposure = character(),
  chip = character(),
  outcome = character(),
  n_nochip = numeric(),
  hr_nochip = numeric(),
  se_nochip = numeric(),
  z_nochip = numeric(),
  pval_nochip = numeric(),
  lci_nochip = numeric(),
  uci_nochip = numeric(),
  n_chip = numeric(),
  hr_chip = numeric(),
  se_chip = numeric(),
  z_chip = numeric(),
  pval_chip = numeric(),
  lci_chip = numeric(),
  uci_chip = numeric(),
  beta_inter = numeric(),
  se_inter = numeric(),
  z_inter = numeric(),
  pval_inter = numeric(),
  stringsAsFactors = FALSE
)

for (outcome in outcome_list) {
  FollowUp_var <- paste0("FollowUp_", outcome)
  Incd_var <- paste0("Incd_", outcome)
  Prev_var <- paste0("Prev_", outcome)
  
  # exclude prevalent cases
  temp_all0 <- merged_data %>%
    filter(.data[[Prev_var]] == 0)
  
  for (current_chip in chip_var) {
    temp_nochip <- filter(temp_all0, .data[[current_chip]] == 0)
    temp_chip <- filter(temp_all0, .data[[current_chip]] == 1)
    
    for (current_exposure in le8_variables) {
      base_covs <- paste(
        "in_white_British_ancestry_subset + age.x + Sex_numeric +",
        paste0("PC", 1:10),
        collapse = " + "
      )
      fmla_no_str <- sprintf(
        "Surv(%s, %s) ~ %s + %s",
        FollowUp_var, Incd_var, current_exposure, base_covs
      )
      fmla_inter_str <- sprintf(
        "Surv(%s, %s) ~ %s * %s + %s",
        FollowUp_var, Incd_var, current_exposure, current_chip, base_covs
      )
      
      m_no <- coxph(as.formula(fmla_no_str), data = temp_nochip)
      m_chip <- coxph(as.formula(fmla_no_str), data = temp_chip)
      m_int <- coxph(as.formula(fmla_inter_str), data = temp_all0)
      
      s_no <- summary(m_no)
      s_chip <- summary(m_chip)
      coef_int <- summary(m_int)$coef
      idx_int <- grep(
        paste0(current_exposure, ":", current_chip),
        rownames(coef_int)
      )
      
      res_row <- data.frame(
        exposure = current_exposure,
        chip = current_chip,
        outcome = outcome,
        n_nochip = s_no$n,
        hr_nochip = s_no$coef[1, "exp(coef)"],
        se_nochip = s_no$coef[1, "se(coef)"],
        z_nochip = s_no$coef[1, "z"],
        pval_nochip = s_no$coef[1, "Pr(>|z|)"],
        lci_nochip = exp(confint(m_no)[1, 1]),
        uci_nochip = exp(confint(m_no)[1, 2]),
        n_chip = s_chip$n,
        hr_chip = s_chip$coef[1, "exp(coef)"],
        se_chip = s_chip$coef[1, "se(coef)"],
        z_chip = s_chip$coef[1, "z"],
        pval_chip = s_chip$coef[1, "Pr(>|z|)"],
        lci_chip = exp(confint(m_chip)[1, 1]),
        uci_chip = exp(confint(m_chip)[1, 2]),
        beta_inter = if (length(idx_int) == 1) coef_int[idx_int, "coef"] else NA,
        se_inter = if (length(idx_int) == 1) coef_int[idx_int, "se(coef)"] else NA,
        z_inter = if (length(idx_int) == 1) coef_int[idx_int, "z"] else NA,
        pval_inter = if (length(idx_int) == 1) coef_int[idx_int, "Pr(>|z|)"] else NA,
        stringsAsFactors = FALSE
      )
      
      results_strat <- bind_rows(results_strat, res_row)
    }
  }
}


# Write stratified Cox results
write.table(
  results_strat,
  file = "/medpop/esp2/yxliu/CHIP_LE8/Stratified_analysis/Stratified_incidence/stratified_cox_results_VAF2.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# results_stratified <- fread("/medpop/esp2/yxliu/updated_stratified_cox_results_VAF2.txt")

# =======================================
# 5) Absolute Incidence Rates
# =======================================

# Build the incidence_all table
incidence_list <- list()
for (outcome in outcome_list) {
  FU <- paste0("FollowUp_", outcome)
  EV <- paste0("Incd_", outcome)
  PV <- paste0("Prev_", outcome) # Prevalence variable
  
  for (ct in chip_var) {
    df_tmp <- merged_data %>%
      filter(.data[[PV]] == 0) %>% # Exclude prevalent cases based on the specific outcome
      mutate(CHIP_status = ifelse(.data[[ct]] == 1, ct, paste0("No_", ct))) %>%
      group_by(LE8_cat, CHIP_status) %>%
      summarise(
        events = sum(.data[[EV]], na.rm = TRUE),
        py = sum(.data[[FU]], na.rm = TRUE),
        rate = events / py, # per 1 PY
        .groups = "drop"
      ) %>%
      mutate(chip_type = ct, outcome = outcome)
    
    incidence_list[[paste(outcome, ct, sep = "_")]] <- df_tmp
  }
}
incidence_all <- bind_rows(incidence_list) %>%
  # convert LE8_cat to low/intermediate/high
  mutate(
    LE8_cat = factor(
      LE8_cat,
      levels = c("<50", "50–<80", "≥80"),
      labels = c("low", "intermediate", "high")
    )
  )

# Re-scale to per 1000 PY and enforce chip_type order
comp_order <- c("CHIP", "DNMT3A", "non_DNMT3A", "TET2", "ASXL1")
incidence_all <- incidence_all %>%
  mutate(
    rate = rate * 1000, # per 1000 PY
    chip_type = factor(chip_type, levels = comp_order)
  )

incidence_all <- incidence_all %>%
  mutate(
    Status = ifelse(CHIP_status == chip_type, "Carrier", "Non-carrier"),
    Status = factor(Status, levels = c("Carrier", "Non-carrier"))
  )

# =========================================================================
# Calculate CI/P-values for Incidence Rates and Export Table
# =========================================================================

# We use rowwise() to apply poisson.test to each row of the dataframe.
incidence_with_stats <- incidence_all %>%
  rowwise() %>%
  mutate(
    # Run the poisson test for the number of events given the person-years
    poisson_results = list(poisson.test(x = events, T = py)),
    # Extract the p-value
    p_value = poisson_results$p.value,
    # Extract the lower and upper bounds of the 95% CI for the rate (per 1 PY)
    lower_ci_raw = poisson_results$conf.int[1],
    upper_ci_raw = poisson_results$conf.int[2]
  ) %>%
  ungroup() %>% # It's important to ungroup after a rowwise operation
  # Scale the CIs to be per 100 Person-Years, matching the 'rate' column
  mutate(
    lower_ci_1000py = lower_ci_raw * 1000,
    upper_ci_1000py = upper_ci_raw * 1000,
    # Create a formatted string for the 95% CI for easy reading
    ci_95_str = paste0("(", sprintf("%.2f", lower_ci_1000py), " - ", sprintf("%.2f", upper_ci_1000py), ")")
  )

# Now, create the final output table by selecting and renaming columns as requested
# and filtering for only the three outcomes of interest.
final_summary_table <- incidence_with_stats %>%
  filter(
    outcome %in% c(
      "Coronary_Artery_Disease_INTERMEDIATE",
      "composite_mi_cad_stroke_death"
    )
  ) %>%
  select(
    exposure = LE8_cat,
    chip_type,
    chip_status = Status,
    outcome,
    incidence_rate_1000py = rate,
    p_value,
    ci_95 = ci_95_str
  ) %>%
  # Arrange the table for better readability
  arrange(outcome, chip_type, exposure, chip_status)


# Write the final summary table to a tab-separated file
write.table(
  final_summary_table,
  file = "/medpop/esp2/yxliu/CHIP_LE8/Stratified_analysis/Stratified_incidence/incidence_rates_summary.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

#print("Successfully generated the incidence rate summary table.")
#print(head(final_summary_table))

#Incidence <- fread("/medpop/esp2/yxliu/incidence_rates_summary_with_stats.txt")

# =======================================
# 6) Combined incidence bar-plots by LE8_cat & CHIP subtype
#    for three outcomes:
#      • Intermediate CAD
#      • Composite CVD (MI + CAD + stroke + death)
#      • CAD death
# =======================================

# 6.1) Standardize Status labels
incidence_all2 <- incidence_all %>%
  mutate(
    Status = ifelse(CHIP_status == chip_type, "Carrier", "Non-carrier"),
    Status = factor(Status, levels = c("Carrier", "Non-carrier"))
  )


# 6.2) Create output folder if needed
output_path <- "/medpop/esp2/yxliu/CHIP_LE8/Stratified_analysis/Stratified_incidence/Plots"
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

# 6.3) Helper to plot & save each outcome
plot_incidence <- function(df, code, title, fname) {
  # 1. filter for this outcome
  dfp <- df %>%
    filter(outcome == code)
  if (nrow(dfp) == 0) {
    warning("No data for outcome: ", code)
    return(NULL)
  }
  
  # 2. Use the provided title directly
  title_expr <- title
  
  # 3. make the plot
  p <- ggplot(dfp, aes(x = LE8_cat, y = rate, fill = Status)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_text(aes(label = sprintf("%.1f", rate)),
              position = position_dodge(width = 0.7),
              vjust = -0.3, size = 3
    ) +
    facet_wrap(~chip_type, nrow = 1) +
    scale_fill_manual(
      values = c(
        "Carrier" = "#1b9e77",
        "Non-carrier" = "#d95f02"
      ),
      name = "",
      labels = c("Carrier", "Non-carrier")
    ) +
    labs(
      title = title_expr,
      x = "LE8 Category",
      y = "Incidence Rate (per 100 PY)"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      axis.text.x = element_text(angle = 30, hjust = 1),
      strip.text = element_text(face = "italic")
    )
  
  # 4. save to disk
  ggsave(
    filename = file.path(output_path, fname),
    plot = p,
    width = 12,
    height = 4,
    dpi = 300
  )
}

# 6.4) Draw & save all three plots

# Intermediate CAD
plot_incidence(
  incidence_all2,
  "Coronary_Artery_Disease_INTERMEDIATE",
  "Incidence of Intermediate CAD by LE8 Category & CHIP Subtype",
  "incidence_CAD_intermediate_by_LE8_and_chip.png"
)

# Composite CVD (MI + CAD + stroke + death)
plot_incidence(
  incidence_all2,
  "composite_mi_cad_stroke_death",
  "Incidence of Composite CVD Outcome by LE8 Category & CHIP Subtype",
  "incidence_composite_CVD_by_LE8_and_chip.png"
)


#######################################################
############## Make two plots together ################# 
#######################################################

library(patchwork)

output_path <- "/medpop/esp2/yxliu/CHIP_LE8/Stratified_analysis/Stratified_incidence/Plots"
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

plot_incidence_single <- function(df, code, title, subplot_label) {
  # 1. filter for this outcome
  dfp <- df %>%
    filter(outcome == code)
  if (nrow(dfp) == 0) {
    warning("No data for outcome: ", code)
    return(NULL)
  }
  
  dfp <- dfp %>%
    mutate(
      chip_type_label = case_when(
        chip_type == "CHIP" ~ "italic('CHIP')",
        chip_type == "ASXL1" ~ "italic('ASXL1')",
        chip_type == "TET2" ~ "italic('TET2')",
        chip_type == "DNMT3A" ~ "italic('DNMT3A')",
        chip_type == "non_DNMT3A" ~ "plain('non-')*italic('DNMT3A')", # More complex for non-DNMT3A
        TRUE ~ chip_type
      ),
      chip_type_label = factor(chip_type_label, levels = c("italic('CHIP')", "italic('DNMT3A')", "plain('non-')*italic('DNMT3A')", "italic('TET2')", "italic('ASXL1')")) # Ensure order
    )
  
  # 3. make the plot
  p <- ggplot(dfp, aes(x = LE8_cat, y = rate, fill = Status)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_text(aes(label = sprintf("%.1f", rate)),
              position = position_dodge(width = 0.7),
              vjust = -0.3, size = 3.5, color = "black" 
    ) +
    facet_wrap(~chip_type_label, nrow = 1, labeller = label_parsed) + 
    scale_fill_manual(
      values = c(
        "Carrier" = "#1b9e77", 
        "Non-carrier" = "#d95f02" 
      ),
      name = "CHIP Status", 
      labels = c("Carrier", "Non-carrier")
    ) +
    labs(
      title = paste0(subplot_label, " ", title), 
      x = "LE8 Category",
      y = "Incidence Rate (per 1000 PY)"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + 
    theme_minimal(base_size = 14) + 
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 30, hjust = 1, size = 11),
      axis.text.y = element_text(size = 11),
      axis.title = element_text(size = 13, face = "bold"),
      strip.text = element_text(face = "bold.italic", size = 12, color = "darkblue"), 
      panel.spacing.x = unit(1, "lines") #
    )
  
  return(p) 
}


# Intermediate CAD
p_intermediate <- plot_incidence_single(
  incidence_all2,
  "Coronary_Artery_Disease_INTERMEDIATE",
  "Incidence of CAD by LE8 Category & CHIP Subtype",
  "(A)"
)

# Composite CVD (MI + CAD + stroke + death)
p_composite_cvd <- plot_incidence_single(
  incidence_all2,
  "composite_mi_cad_stroke_death",
  "Incidence of Composite CVD Outcome by LE8 Category & CHIP Subtype",
  "(B)"
)



# Stack plots vertically
combined_plot <- (p_intermediate / p_composite_cvd) +
  plot_layout(guides = "collect") + 
  plot_annotation(
    title = 'Incidence Rates by LE8 Category and CHIP Subtype for Cardiovascular Outcomes',
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  )

ggsave(
  filename = file.path(output_path, "combined_incidence_plots_aesthetic.png"),
  plot = combined_plot,
  width = 16, # Adjust width for better readability of multiple facets
  height = 10, # Adjust height for three stacked plots
  dpi = 300
)

message("Combined incidence plot saved to: ", file.path(output_path, "combined_incidence_plots_aesthetic.png"))





################# LE8 Individual Component #################

# =======================================
# 1) Data Input & Preparation
# =======================================


# 0) Load required libraries
library(data.table)
library(dplyr)
library(survival)
library(broom)
library(tidyr)
library(ggplot2)
library(stringr)

setwd("/medpop/esp2")

# =======================================
# 1) Data Input & Preparation
# =======================================

# Load LE8 composite score data
le8_data <- fread("/medpop/esp2/yxliu/LE8score_calculation/updated_LE8score.txt")

# Load CHIP calls
chip_data <- fread("/medpop/esp2/tnakao/data/UKBB/phenotype/CHIP/450k/2022Oct/CHIP_calls_Oct16_2022.txt")

# Load the "mask" files for filtering
id_450k <- fread("/medpop/esp2/mesbah/datasets/CHIP/UKBB/450k/ukb450k.eid.list") %>% rename(id = V1)
no_consent <- fread("/medpop/esp2/projects/UK_Biobank/withdrawn_samples/w7089_20241217.csv") %>% rename(id = V1)
relate_remove <- fread("/medpop/esp2/zyu/chip_protemoics/listofremovefor450K.txt") %>% rename(id = IID)

# This file also needs to contain the baseline variables (age, sex, PCs, etc.)
pheno_master <- fread("/medpop/esp2/zyu/UKB_bigquery_2023pheno/20250711_pheno_noGP_clean.tsv.gz") %>%
  as.data.frame()

ori_pheno <- fread("/medpop/esp2/mzekavat/UKBB/ukbb_PhenoFile.ALL_500k.UpdatedIncdPhenos_202020.REAL.txt") %>%
  as.data.frame() %>%
  filter(Submitted_Gender == Inferred_Gender,
         Non_Consented     == 0,
         Prev_Myeloid_leukemia  == 0,
         Prev_Lymphoid_leukemia == 0
  ) %>%
  select(
    id, age, Sex_numeric, in_white_British_ancestry_subset,
    paste0("PC",1:10)
  )

pheno_master2 <- merge(
  pheno_master, ori_pheno,
  by = "id",
  all = TRUE
)

# Subset pheno_master to just the 450K participants with extensive QC
# and rename outcome columns to match existing code's expected pattern
pheno_filtered_renamed <- pheno_master2 %>%
  # keep only those in the 450K ID list
  filter(id %in% id_450k$id) %>%
  # remove withdrawn consent
  filter(!id %in% no_consent$id) %>%
  # only unrelated individuals, drop related ones
  filter(!id %in% relate_remove$id) %>%
  # apply existing QC steps (assuming these fields exist in pheno_master)
  # Select baseline columns and all relevant outcome columns
  select(
    id,
    age, 
    Sex_numeric, 
    in_white_British_ancestry_subset, 
    paste0("PC", 1:10), 
    # Select and rename prevalence flags
    Prev_composite_mi_cad_stroke_death = prevalent_disease_Coronary_Artery_Disease_SOFT,
    Prev_Coronary_Artery_Disease_INTERMEDIATE = prevalent_disease_Coronary_Artery_Disease_INTERMEDIATE,
    Prev_Coronary_Artery_Disease_Death = prevalent_disease_CAD_Death,
    # Select and rename follow-up times
    FollowUp_composite_mi_cad_stroke_death = followup_Coronary_Artery_Disease_SOFT,
    FollowUp_Coronary_Artery_Disease_INTERMEDIATE = followup_Coronary_Artery_Disease_INTERMEDIATE,
    FollowUp_Coronary_Artery_Disease_Death = followup_CAD_Death,
    # Select and rename event indicators
    Incd_composite_mi_cad_stroke_death = incident_disease_Coronary_Artery_Disease_SOFT,
    Incd_Coronary_Artery_Disease_INTERMEDIATE = incident_disease_Coronary_Artery_Disease_INTERMEDIATE,
    Incd_Coronary_Artery_Disease_Death = incident_disease_CAD_Death
  )


# Merge LE8 and the fully processed phenotype data (inner-join so only complete cases remain)
merged_data <- merge(
  le8_data, pheno_filtered_renamed,
  by = "id",
  all = FALSE
)


# =======================================
# 2) Process CHIP Data (5 subtypes)
# =======================================
CHIP_ids      <- unique(chip_data[AF >= 0.02,]$Broad_ID)
ASXL1_ids     <- unique(chip_data[Gene.refGene=="ASXL1"  & AF >= 0.02,]$Broad_ID)
TET2_ids      <- unique(chip_data[Gene.refGene=="TET2"   & AF >= 0.02,]$Broad_ID)
DNMT3A_ids    <- unique(chip_data[Gene.refGene=="DNMT3A" & AF >= 0.02,]$Broad_ID)
nonDNMT3A_ids <- unique(chip_data[AF >= 0.02 & Gene.refGene!="DNMT3A",]$Broad_ID)

merged_data <- merged_data %>%
  mutate(
    CHIP       = as.integer(id %in% CHIP_ids),
    ASXL1      = as.integer(id %in% ASXL1_ids),
    TET2       = as.integer(id %in% TET2_ids),
    DNMT3A     = as.integer(id %in% DNMT3A_ids),
    non_DNMT3A = as.integer(id %in% nonDNMT3A_ids)
  )

# =======================================
# 3) Process INDIVIDUAL LE8 Component Exposures
# =======================================

le8_component_vars <- c(
  "diet_score", "BP_Score", "PA_score", "final_smoking_score",
  "sleep_health_points", "bmi_points", "Blood_lipid_points", "HbA1c_score"
)

merged_data <- merged_data %>%
  mutate(
    across(
      all_of(le8_component_vars),
      ~ . / 10,
      .names = "{.col}_scaled" 
    ),
    across(
      all_of(le8_component_vars),
      ~ case_when(
        . < 50 ~ "Low",
        . < 80 ~ "Intermediate",
        TRUE   ~ "High"
      ) %>% factor(levels = c("Low", "Intermediate", "High")),
      .names = "{.col}_cat" # Creates diet_score_cat, BP_Score_cat, etc.
    )
  )


# ======================================================
# 4) Stratified Cox Regression (for each LE8 component)
# ======================================================
outcome_list <- c(
  "composite_mi_cad_stroke_death",
  "Coronary_Artery_Disease_INTERMEDIATE"
)

chip_var        <- c("CHIP","ASXL1","TET2","DNMT3A","non_DNMT3A")

le8_variables_scaled <- paste0(le8_component_vars, "_scaled")

results_strat <- data.frame() # Initialize empty data frame

for (outcome in outcome_list) {
  FollowUp_var <- paste0("FollowUp_", outcome)
  Incd_var     <- paste0("Incd_",     outcome)
  Prev_var     <- paste0("Prev_",     outcome)
  
  temp_all0 <- merged_data %>% filter(.data[[Prev_var]] == 0)
  
  for (current_chip in chip_var) {
    temp_nochip <- filter(temp_all0, .data[[current_chip]] == 0)
    temp_chip   <- filter(temp_all0, .data[[current_chip]] == 1)
    
    for (current_exposure in le8_variables_scaled) {
      base_covs <- paste(
        "in_white_British_ancestry_subset + age.x + Sex_numeric +",
        paste0("PC", 1:10, collapse=" + ")
      )
      fmla_no_str <- sprintf("Surv(%s, %s) ~ %s + %s", FollowUp_var, Incd_var, current_exposure, base_covs)
      fmla_inter_str <- sprintf("Surv(%s, %s) ~ %s * %s + %s", FollowUp_var, Incd_var, current_exposure, current_chip, base_covs)
      
      m_no   <- coxph(as.formula(fmla_no_str), data = temp_nochip)
      m_chip <- coxph(as.formula(fmla_no_str), data = temp_chip)
      m_int  <- coxph(as.formula(fmla_inter_str), data = temp_all0)
      
      s_no    <- summary(m_no)
      s_chip  <- summary(m_chip)
      coef_int <- summary(m_int)$coef
      idx_int  <- grep(paste0(current_exposure, ":", current_chip), rownames(coef_int))
      
      res_row <- data.frame(
        exposure    = str_remove(current_exposure, "_scaled"), 
        chip        = current_chip,
        outcome     = outcome,
        n_nochip    = s_no$n,
        hr_nochip   = s_no$coef[1,"exp(coef)"], se_nochip   = s_no$coef[1,"se(coef)"],
        pval_nochip = s_no$coef[1,"Pr(>|z|)"], lci_nochip  = exp(confint(m_no)[1,1]), uci_nochip  = exp(confint(m_no)[1,2]),
        n_chip      = s_chip$n,
        hr_chip     = s_chip$coef[1,"exp(coef)"], se_chip     = s_chip$coef[1,"se(coef)"],
        pval_chip   = s_chip$coef[1,"Pr(>|z|)"], lci_chip    = exp(confint(m_chip)[1,1]), uci_chip    = exp(confint(m_chip)[1,2]),
        pval_inter  = if(length(idx_int)==1) coef_int[idx_int,"Pr(>|z|)"] else NA
      )
      
      results_strat <- bind_rows(results_strat, res_row)
    }
  }
}

# Write stratified Cox results for all components
write.table(
  results_strat,
  file      = "/medpop/esp2/yxliu/CHIP_LE8/Stratified_analysis/Stratified_incidence/stratified_cox_results_component.txt",
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)


# =======================================
# Load Pre-computed Stratified Cox Regression Results (SKIPPING THE LONG LOOP)
# =======================================
# cat("Skipping the long Cox regression loop and loading results from file...\n")
# results_strat <- fread("/medpop/esp2/yxliu/stratified_cox_results_by_component.txt")
# cat("Cox regression results loaded successfully.\n")



# =======================================
# 4.5-) Absolute Incidence Rates (for each component) - MEMORY EFFICIENT VERSION
# =======================================

cat("Calculating incidence rates using the memory-efficient method...\n")

# Define the two main outcomes for which plots will be generated
outcome_list <- c(
  "composite_mi_cad_stroke_death", # Represents Composite CVD outcome
  "Coronary_Artery_Disease_INTERMEDIATE" # Represents CAD outcome
)

# Define the five CHIP mutation types for subplots
chip_var <- c("CHIP", "ASXL1", "TET2", "DNMT3A", "non_DNMT3A")

# Define the eight LE8 components, each of which will have its own plot
le8_component_vars <- c(
  "diet_score", "BP_Score", "PA_score", "final_smoking_score",
  "sleep_health_points", "bmi_points", "Blood_lipid_points", "HbA1c_score"
)

incidence_list <- list()

# Loop through each outcome (e.g., CAD, Composite CVD)
for (outcome in outcome_list) {
  df_outcome_filtered <- merged_data %>%
    filter(.data[[paste0("Prev_", outcome)]] == 0)
  
  for (component in le8_component_vars) {
    
    component_cat_var <- paste0(component, "_cat")
    
    for (ct in chip_var) {
      
      df_tmp <- df_outcome_filtered %>%
        mutate(CHIP_status = ifelse(.data[[ct]] == 1, ct, paste0("No_", ct))) %>%
        group_by(le8_category = .data[[component_cat_var]], CHIP_status) %>%
        summarise(
          events = sum(.data[[paste0("Incd_", outcome)]], na.rm = TRUE), # Sum of incident events
          py     = sum(.data[[paste0("FollowUp_", outcome)]], na.rm = TRUE), # Sum of person-years
          .groups = "drop" # Drop grouping to avoid subsequent issues
        ) %>%
        # Add identifiers to track which component, chip type, and outcome these results belong to
        mutate(
          le8_component = component, # Name of the current LE8 component
          chip_type = ct, # Name of the current CHIP mutation type
          outcome = outcome # Name of the current outcome
        )
      
      # Store the results in a list, using a unique key for each combination
      key <- paste(outcome, component, ct, sep = "_")
      incidence_list[[key]] <- df_tmp
    }
  }
}

# Combine all individual results into a single data frame for plotting and summary
incidence_all2 <- bind_rows(incidence_list) %>%
  mutate(
    rate      = (events / py) * 1000, # Calculate incidence rate per 1000 Person-Years
    Status    = factor(ifelse(CHIP_status == chip_type, "Carrier", "Non-carrier"), levels = c("Carrier", "Non-carrier")),
    # Ensure chip_type is a factor with the desired order for consistent subplot arrangement
    chip_type = factor(chip_type, levels = c("CHIP", "DNMT3A", "non_DNMT3A", "TET2", "ASXL1")),
    # Ensure LE8 category is a factor with the desired order for x-axis
    le8_category = factor(le8_category, levels = c("Low", "Intermediate", "High"))
  )


# =========================================================================
# 5-) Generate Final Summary Table with Interaction P-values
# =========================================================================
# This section will:
# 1. Calculate 95% CIs and p-values for the incidence rates using a Poisson model.
# 2. Extract the interaction p-values from the Cox model results (assuming results_strat is loaded).
# 3. Join these two pieces of information into a single summary table.
# 4. Save the final table to a file.

cat("Generating final summary table with incidence rates and interaction p-values...\n")

# --- Step 1: Calculate CIs and p-values for the incidence rates ---
incidence_with_stats <- incidence_all2 %>%
  # Filter out groups with zero person-years to prevent 0/0 -> NaN errors
  filter(py > 0) %>%
  rowwise() %>%
  mutate(
    # Run the poisson test for the number of events given the person-years
    poisson_results = list(poisson.test(x = events, T = py)),
    # Extract the p-value for the rate
    p_value_rate = poisson_results$p.value,
    # Extract the lower and upper bounds of the 95% CI (per 1 PY)
    lower_ci_raw = poisson_results$conf.int[1],
    upper_ci_raw = poisson_results$conf.int[2]
  ) %>%
  ungroup() %>%
  # Scale CIs to be per 1000 Person-Years and create a formatted string
  mutate(
    lower_ci_1000py = lower_ci_raw * 1000,
    upper_ci_1000py = upper_ci_raw * 1000,
    ci_95_str = paste0("(", sprintf("%.2f", lower_ci_1000py), " - ", sprintf("%.2f", upper_ci_1000py), ")")
  )

# --- Step 2: Prepare the interaction p-values for joining ---

if (exists("results_strat")) { 
  interaction_pvals <- results_strat %>%
    select(
      le8_component = exposure, 
      chip_type = chip,          
      outcome,
      pval_interaction = pval_inter
    )
} else {
  warning("results_strat not found. Interaction p-values will not be included in the summary table.")
  interaction_pvals <- data.frame(le8_component = character(), chip_type = character(), outcome = character(), pval_interaction = numeric())
}


# --- Step 3: Join the incidence stats with the interaction p-values ---

final_summary_table <- left_join(
  incidence_with_stats,
  interaction_pvals,
  by = c("le8_component", "chip_type", "outcome")
)

# --- Step 4: Finalize the table and save to disk ---
# Filter for the two requested outcomes, select and rename columns, and arrange for clarity.
final_output <- final_summary_table %>%
  filter(
    outcome %in% c(
      "Coronary_Artery_Disease_INTERMEDIATE",
      "composite_mi_cad_stroke_death"
    )
  ) %>%
  select(
    exposure_component = le8_component,
    exposure_category = le8_category,
    chip_type,
    chip_status = Status,
    outcome,
    incidence_rate_1000py = rate, # Renamed for clarity to match rate unit
    ci_95 = ci_95_str,
    p_value_rate,
    pval_interaction
  ) %>%
  # Arrange the table for better readability
  arrange(outcome, exposure_component, chip_type, exposure_category, chip_status)

# Write the final summary table to a tab-separated file
output_file_path <- "/medpop/esp2/yxliu/CHIP_LE8/Stratified_analysis/Stratified_incidence/component_incidence_interaction.txt"
write.table(
  final_output,
  file = output_file_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Successfully generated and saved the summary table to:", output_file_path, "\n")
# print(head(final_output))


# =======================================
# 6) Combined incidence bar-plots (16 plots total)
# =======================================

output_path_components <- "/medpop/esp2/yxliu/CHIP_LE8/Stratified_analysis/Stratified_incidence/Plots_by_Component"
if (!dir.exists(output_path_components)) dir.create(output_path_components, recursive = TRUE)

# A more flexible plotting function
plot_incidence_by_component <- function(df, plot_title, file_name) {
  # Check if data exists for plotting
  if (nrow(df) == 0) {
    warning("No data for this plot: ", file_name)
    return(NULL)
  }
  
  p <- ggplot(df, aes(x = le8_category, y = rate, fill = Status)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) + 
    geom_text(aes(label = sprintf("%.1f", rate)), 
              position = position_dodge(width = 0.7),
              vjust = -0.3, size = 3) +
    facet_wrap(~ chip_type, nrow = 1) + 
    scale_fill_manual(
      values = c("Carrier" = "#1b9e77", "Non-carrier" = "#d95f02"),
      name = "" # No legend title
    ) +
    labs(
      title = plot_title, # Dynamic plot title
      x     = "LE8 Category", # X-axis label
      y     = "Incidence Rate (per 1000 PY)" # Y-axis label
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + # Expand y-axis to make space for labels
    theme_minimal(base_size = 13) + # Minimal theme with base font size
    theme(
      plot.title      = element_text(face = "bold", hjust = 0.5, size = 16), # Title styling
      panel.grid.major  = element_blank(), # Remove major grid lines
      panel.grid.minor  = element_blank(), # Remove minor grid lines
      legend.position   = "top", # Place legend at the top
      axis.text.x     = element_text(angle = 30, hjust = 1), # Rotate x-axis labels for readability
      strip.text      = element_text(face = "italic", size = 12) # Styling for facet labels (chip types)
    )
  
  # Save the plot to the specified output folder
  ggsave(
    filename = file.path(output_path_components, file_name),
    plot     = p,
    width    = 12, # Plot width
    height   = 5, # Plot height
    dpi      = 300 # Resolution
  )
  
  return(p) 
}


# --- Plotting Loop: Generate all 16 plots ---

outcome_titles <- c(
  "Coronary_Artery_Disease_INTERMEDIATE" = "CAD",
  "composite_mi_cad_stroke_death"        = "Composite CVD"
)

component_titles <- c(
  "diet_score"          = "Diet Score",
  "BP_Score"            = "Blood Pressure Score",
  "PA_score"            = "Physical Activity Score",
  "final_smoking_score" = "Smoking Score",
  "sleep_health_points" = "Sleep Health Score",
  "bmi_points"          = "BMI Score",
  "Blood_lipid_points"  = "Blood Lipids Score",
  "HbA1c_score"         = "Glucose Score"
)

for (outcome_code in names(outcome_titles)) {
  for (component_code in names(component_titles)) {
    
    plot_data <- incidence_all2 %>%
      filter(outcome == outcome_code, le8_component == component_code)
    
    current_plot_title <- sprintf("Incidence of %s by %s & CHIP Subtype",
                                  outcome_titles[outcome_code],
                                  component_titles[component_code])
    
    current_file_name <- sprintf("incidence_%s_by_%s.png",
                                 outcome_code,
                                 component_code)
    
    cat("Generating plot:", current_file_name, "\n")
    plot_incidence_by_component(
      df = plot_data,
      plot_title = current_plot_title,
      file_name = current_file_name
    )
  }
}
