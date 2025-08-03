############################################################
############ clean version - used for manuscript############
############################################################

# =======================================
# Full script: Single-reference “No-CHIP & High” forest plots for 4 outcomes
# ======================================

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
# Full script: Single-reference “No-CHIP & High” forest plots for 4 outcomes
# =======================================


# 1.1 Define CHIP & subtypes (binary indicators)
df <- merged_data %>%
  mutate(
    CHIP       = as.integer(id %in% chip_data[AF >= 0.02, Broad_ID]),
    ASXL1      = as.integer(id %in% chip_data[Gene.refGene=="ASXL1"  & AF >= 0.02, Broad_ID]),
    TET2       = as.integer(id %in% chip_data[Gene.refGene=="TET2"   & AF >= 0.02, Broad_ID]),
    DNMT3A     = as.integer(id %in% chip_data[Gene.refGene=="DNMT3A" & AF >= 0.02, Broad_ID]),
    non_DNMT3A = as.integer(id %in% chip_data[AF >= 0.02 & Gene.refGene!="DNMT3A", Broad_ID])
  )

# 1.2 Define LE8 category: High (≥80), Intermediate (50–<80), Low (<50)
df <- df %>%
  mutate(
    LE8_composite = (
      diet_score + BP_Score + PA_score + final_smoking_score +
        sleep_health_points + bmi_points + Blood_lipid_points + HbA1c_score
    ) / 8,
    LE8_cat = case_when(
      LE8_composite >= 80 ~ "High",
      LE8_composite >= 50 ~ "Intermediate",
      TRUE                ~ "Low"
    ) %>%
      factor(levels = c("High","Intermediate","Low"))
  )

# 1.3 Define single ‘subgroup’ variable —
df <- df %>%
  mutate(
    subgroup = case_when(
      CHIP == 0 & ASXL1==0 & TET2==0 & DNMT3A==0            ~ "No-CHIP",
      ASXL1 == 1                                           ~ "ASXL1",
      TET2  == 1                                           ~ "TET2",
      DNMT3A == 1                                          ~ "DNMT3A",
      non_DNMT3A == 1                                      ~ "non-DNMT3A",
      TRUE                                                 ~ NA_character_
    ) %>%
      factor(levels = c("No-CHIP","DNMT3A","non-DNMT3A","TET2","ASXL1"))
  )

# 1.4 Build combined factor ‘combo’ = subgroup + LE8_cat
combo_levels <- c(
  "No-CHIP_High","No-CHIP_Intermediate","No-CHIP_Low",
  "DNMT3A_Intermediate","DNMT3A_Low",
  "non-DNMT3A_Intermediate","non-DNMT3A_Low",
  "TET2_Intermediate","TET2_Low",
  "ASXL1_Intermediate","ASXL1_Low"
)

df <- df %>%
  mutate(
    combo = paste(subgroup, LE8_cat, sep = "_") %>%
      factor(levels = combo_levels)
  )


# 2) Outcomes & covariates --------------------------------------------------

outcomes <- list(
  "Intermediate CAD"  = c("FollowUp_Coronary_Artery_Disease_INTERMEDIATE",
                          "Incd_Coronary_Artery_Disease_INTERMEDIATE"),
  "Composite CVD"     = c("FollowUp_composite_mi_cad_stroke_death",
                          "Incd_composite_mi_cad_stroke_death")
)

# Base covariates
covariate_string <- paste(
  "in_white_British_ancestry_subset",
  "age.x",
  "Sex_numeric",
  paste0("PC", 1:10),
  sep = " + "
)

# 3) Fit Cox models & extract HRs ------------------------------------------

forest_list <- list()

for (nm in names(outcomes)) {
  FU <- outcomes[[nm]][1]
  EV <- outcomes[[nm]][2]
  
  df_fit <- df %>% filter(!is.na(combo))
  
  fit <- coxph(
    as.formula(paste0("Surv(", FU, ", ", EV, ") ~ combo + ", covariate_string)),
    data = df_fit
  )
  
  td <- tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(grepl("^combo", term)) %>%
    mutate(
      combo   = sub("^combo", "", term),
      outcome = nm
    )
  forest_list[[nm]] <- td
}

forest_df <- bind_rows(forest_list) %>%
  mutate(
    combo    = factor(combo, levels = combo_levels),
    subgroup = factor(sub("_.*", "", combo),
                      levels = c("No-CHIP","CHIP","DNMT3A","non-DNMT3A","TET2","ASXL1")),
    LE8_cat  = factor(sub(".*_", "", combo),
                      levels = c("High","Intermediate","Low"))
  )

write.table(
  forest_df,
  file      = "/medpop/esp2/yxliu/CHIP_LE8/Stratified_analysis/Stratified_HR/Results_HR_stratified_singleRef.txt",
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

#forest_df <- fread("/medpop/esp2/yxliu/Results_HR_stratified.txt")

# 4) Forest-plot function ---------------------------------------------------

# 1) 带缩进的 y 轴标签（包括组名 Header）
axis_levels <- c(
  "No-CHIP_Header",
  "No-CHIP_High","No-CHIP_Intermediate","No-CHIP_Low",
  "DNMT3A_Header",
  "DNMT3A_Intermediate","DNMT3A_Low",
  "non-DNMT3A_Header",
  "non-DNMT3A_Intermediate","non-DNMT3A_Low",
  "TET2_Header",
  "TET2_Intermediate","TET2_Low",
  "ASXL1_Header",
  "ASXL1_Intermediate","ASXL1_Low"
)

labels_combo <- setNames(
  c(
    expression(bold(plain("No-"))*italic("CHIP")),              # No-CHIP_Header
    "   High",
    "   Intermediate",
    "   Low",
    expression(bold(italic("DNMT3A"))),                            # DNMT3A_Header
    "   Intermediate",
    "   Low",
    expression(bold(plain("non-"))*italic("DNMT3A")),               # non-DNMT3A_Header
    "   Intermediate",
    "   Low",
    expression(bold(italic("TET2"))),                              # TET2_Header
    "   Intermediate",
    "   Low",
    expression(bold(italic("ASXL1"))),                             # ASXL1_Header
    "   Intermediate",
    "   Low"
  ),
  axis_levels
)

plot_forest <- function(df, title_text) {
  header_df <- data.frame(
    combo     = axis_levels[grep("_Header$", axis_levels)],
    estimate  = NA_real_,
    conf.low  = NA_real_,
    conf.high = NA_real_,
    LE8_cat   = factor(NA, levels=c("High","Intermediate","Low")),
    stringsAsFactors = FALSE
  )
  
  full_df <- bind_rows(
    header_df,
    df %>% select(combo, estimate, conf.low, conf.high, LE8_cat)
  ) %>%
    mutate(combo = factor(combo, levels = axis_levels))
  
  ggplot() +
    geom_point(
      data = filter(full_df, !is.na(estimate)),
      aes(x = estimate, y = combo, color = LE8_cat),
      size = 3, position = position_dodge(width = 0.7)
    ) +
    geom_errorbarh(
      data = filter(full_df, !is.na(estimate)),
      aes(xmin = conf.low, xmax = conf.high, y = combo, color = LE8_cat),  # <--- 映射 color
      height = 0.2, position = position_dodge(width = 0.7)
    ) +
    # 虚线 x=1
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
    # y 轴顺序 & 标签
    scale_y_discrete(
      limits = rev(axis_levels),
      labels = rev(labels_combo[axis_levels])
    ) +
    scale_x_log10() +
    scale_color_manual(values = c("Intermediate" = "#E69F00", "Low" = "#56B4E9")) +
    labs(
      title = title_text,
      x     = "Hazard Ratio (log scale)",
      y     = NULL,
      color = "LE8 Category"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title         = element_text(face = "bold", size = 16, hjust = 0.5),
      axis.text.y        = element_text(hjust = 0),
      axis.text.x        = element_text(size = 12),
      panel.grid.major.x = element_line(color = "gray80"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.position    = "bottom",
      legend.key.width   = grid::unit(1.5, "lines")
    )
}


library(patchwork)

p1 <- plot_forest(
  forest_df %>% filter(outcome == "Intermediate CAD"),
  "Hazard Ratios for CAD"
)
p3 <- plot_forest(
  forest_df %>% filter(outcome == "CAD Death (HARD)"),
  "Hazard Ratios for CAD Death"
)
p2 <- plot_forest(
  forest_df %>% filter(outcome == "Composite CVD"),
  "Hazard Ratios for CVD"
)

# ====== Combine the plots and add labels, then collect legend ======

p1_no_legend <- p1 + theme(legend.position = "none")
p3_no_legend <- p3 + theme(legend.position = "none")
p2_no_legend <- p2 + theme(legend.position = "none")

# Combine the plots without individual legends
combined_plot_no_legend <- p1_no_legend + p3_no_legend + p2_no_legend +
  plot_layout(ncol = 1) + 
  plot_annotation(tag_levels = '(A)') & 
  theme(plot.tag = element_text(face = 'bold')) 

combined_plot <- (p1 + p3 + p2) + 
  plot_layout(
    ncol = 1,
    guides = 'collect' 
  ) +
  plot_annotation(
    tag_levels = 'A',        
    tag_prefix = '(',      
    tag_suffix = ')'         
  ) &
  theme(
    plot.tag = element_text(face = 'bold'), 
    legend.position = 'bottom',             
    legend.title = element_text(face = "bold"), 
    legend.text = element_text(size = 9)    
  )

print(combined_plot)

# ====== Save the combined plot to disk ======
ggsave(
  "/medpop/esp2/yxliu/CHIP_LE8/Stratified_analysis/Stratified_HR/Plots_SingleRef/2HR_Combined_Outcomes_SingleLegend2.png",
  combined_plot,
  width = 8,
  height = 15 
)



###################################################################
###################################################################
##### Unique ref = No-CHIP_High + within-subgroup comparisons #####
###################################################################
###################################################################


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
le8_data <- fread("/medpop/esp2/yxliu/updated_LE8score.txt")

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


# 1.1 Define CHIP & subtypes (binary indicators)
df <- merged_data %>%
  mutate(
    CHIP       = as.integer(id %in% chip_data[AF >= 0.02, Broad_ID]),
    ASXL1      = as.integer(id %in% chip_data[Gene.refGene=="ASXL1"  & AF >= 0.02, Broad_ID]),
    TET2       = as.integer(id %in% chip_data[Gene.refGene=="TET2"   & AF >= 0.02, Broad_ID]),
    DNMT3A     = as.integer(id %in% chip_data[Gene.refGene=="DNMT3A" & AF >= 0.02, Broad_ID]),
    non_DNMT3A = as.integer(id %in% chip_data[AF >= 0.02 & Gene.refGene!="DNMT3A", Broad_ID])
  )

# 1.2 Define LE8 category: High (≥80), Intermediate (50–<80), Low (<50)
df <- df %>%
  mutate(
    LE8_composite = (
      diet_score + BP_Score + PA_score + final_smoking_score +
        sleep_health_points + bmi_points + Blood_lipid_points + HbA1c_score
    ) / 8,
    LE8_cat = case_when(
      LE8_composite >= 80 ~ "High",
      LE8_composite >= 50 ~ "Intermediate",
      TRUE                ~ "Low"
    ) %>%
      factor(levels = c("High","Intermediate","Low"))
  )


# 1.3 Define single ‘subgroup’ variable —
df <- df %>%
  mutate(
    subgroup = case_when(
      CHIP == 0 & ASXL1==0 & TET2==0 & DNMT3A==0            ~ "No-CHIP",
      ASXL1 == 1                                           ~ "ASXL1",
      TET2  == 1                                           ~ "TET2",
      DNMT3A == 1                                          ~ "DNMT3A",
      non_DNMT3A == 1                                      ~ "non-DNMT3A",
      TRUE                                                 ~ NA_character_
    ) %>%
      factor(levels = c("No-CHIP","DNMT3A","non-DNMT3A","TET2","ASXL1"))
  )

# 1.4 Build combined factor ‘combo’ = subgroup + LE8_cat
# MODIFICATION: Define levels to include ALL combinations that will appear on the plot
combo_levels <- c(
  "No-CHIP_High", "No-CHIP_Intermediate", "No-CHIP_Low",
  "DNMT3A_High", "DNMT3A_Intermediate", "DNMT3A_Low",
  "non-DNMT3A_High", "non-DNMT3A_Intermediate", "non-DNMT3A_Low",
  "TET2_High", "TET2_Intermediate", "TET2_Low",
  "ASXL1_High", "ASXL1_Intermediate", "ASXL1_Low"
)

df <- df %>%
  mutate(
    combo = paste(subgroup, LE8_cat, sep = "_") %>%
      factor(levels = combo_levels)
  )


# 2) Outcomes & covariates --------------------------------------------------

outcomes <- list(
  "Intermediate CAD"  = c("FollowUp_Coronary_Artery_Disease_INTERMEDIATE",
                          "Incd_Coronary_Artery_Disease_INTERMEDIATE"),
  "Composite CVD"     = c("FollowUp_composite_mi_cad_stroke_death",
                          "Incd_composite_mi_cad_stroke_death"),
  "CAD Death (HARD)"  = c("FollowUp_Coronary_Artery_Disease_Death",
                          "Incd_Coronary_Artery_Disease_Death")
)

# Base covariates
covariate_string <- paste(
  "in_white_British_ancestry_subset",
  "age.x", 
  paste0("PC", 1:10, collapse = " + "),
  sep = " + "
)

# 3) Fit Cox models & extract HRs ------------------------------------------

all_results_list <- list()

# Get list of CHIP subtypes to loop through for specific comparisons
chip_subgroups <- c("DNMT3A", "non-DNMT3A", "TET2", "ASXL1")

for (nm in names(outcomes)) {
  FU <- outcomes[[nm]][1]
  EV <- outcomes[[nm]][2]
  
  df_fit <- df %>% filter(!is.na(combo))
  
  # --- Step 3.1: Main model with No-CHIP_High as the reference ---
  # This gives us HRs for:
  #   - No-CHIP Intermediate/Low vs. No-CHIP High
  #   - Each CHIP-subtype High vs. No-CHIP High
  
  fit_main <- coxph(
    as.formula(paste0("Surv(", FU, ", ", EV, ") ~ combo + ", covariate_string)),
    data = df_fit
  )
  
  # Extract results from the main model
  results_main <- tidy(fit_main, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(grepl("^combo", term)) %>%
    mutate(term = sub("^combo", "", term)) %>%
    # Keep only the comparisons relative to the main reference we need
    filter(
      grepl("^No-CHIP_", term) | grepl("_High$", term)
    )
  
  # Add the reference group manually
  ref_row <- tibble(
    term = "No-CHIP_High", estimate = 1, std.error = NA,
    statistic = NA, p.value = NA, conf.low = 1, conf.high = 1
  )
  
  results_main <- bind_rows(ref_row, results_main)
  
  # --- Step 3.2: Loop through CHIP subtypes to get within-subgroup comparisons ---
  # This gives us HRs for:
  #   - DNMT3A Intermediate/Low vs. DNMT3A High
  #   - TET2 Intermediate/Low vs. TET2 High
  #   - etc.
  
  results_within_subgroups <- list()
  for (subgroup_name in chip_subgroups) {
    
    # Define the reference level for this specific subgroup model
    ref_level <- paste0(subgroup_name, "_High")
    
    # Create a temporary factor with the correct reference level
    df_fit_releveled <- df_fit %>%
      mutate(combo_releveled = relevel(combo, ref = ref_level))
    
    # Fit the model with the new reference level
    fit_subgroup <- coxph(
      as.formula(paste0("Surv(", FU, ", ", EV, ") ~ combo_releveled + ", covariate_string)),
      data = df_fit_releveled
    )
    
    # Tidy results and keep only the Intermediate/Low for the current subgroup
    results_sub <- tidy(fit_subgroup, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(grepl("^combo_releveled", term)) %>%
      mutate(term = sub("^combo_releveled", "", term)) %>%
      filter(
        grepl(paste0("^", subgroup_name, "_Intermediate$"), term) |
          grepl(paste0("^", subgroup_name, "_Low$"), term)
      )
    
    results_within_subgroups[[subgroup_name]] <- results_sub
  }
  
  # --- Step 3.3: Combine all results for this outcome ---
  
  final_df_for_outcome <- bind_rows(
    results_main,
    bind_rows(results_within_subgroups)
  ) %>%
    mutate(
      combo   = term,
      outcome = nm
    ) %>%
    select(outcome, combo, estimate, conf.low, conf.high)
  
  all_results_list[[nm]] <- final_df_for_outcome
}

# Bind all results from all outcomes into one final data frame
forest_df <- bind_rows(all_results_list) %>%
  mutate(
    combo    = factor(combo, levels = combo_levels),
    subgroup = factor(sub("_.*", "", combo),
                      levels = c("No-CHIP", "DNMT3A", "non-DNMT3A", "TET2", "ASXL1")),
    LE8_cat  = factor(sub(".*_", "", combo),
                      levels = c("High", "Intermediate", "Low"))
  )

# Write results to file
write.table(
  forest_df,
  file      = "/medpop/esp2/yxliu/CHIP_LE8/Stratified_analysis/Stratified_HR/Results_HR_stratified_uniRef.txt",
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

# Optional: read back in for plotting if running script in sections
# forest_df <- fread("/medpop/esp2/yxliu/Results_HR_stratified_v2.txt")


# 4) Forest-plot function ---------------------------------------------------

# 4.1 Define y-axis levels and labels, now including all groups
axis_levels <- c(
  "No-CHIP_Header",
  "No-CHIP_High", "No-CHIP_Intermediate", "No-CHIP_Low",
  "DNMT3A_Header",
  "DNMT3A_High", "DNMT3A_Intermediate", "DNMT3A_Low",
  "non-DNMT3A_Header",
  "non-DNMT3A_High", "non-DNMT3A_Intermediate", "non-DNMT3A_Low",
  "TET2_Header",
  "TET2_High", "TET2_Intermediate", "TET2_Low",
  "ASXL1_Header",
  "ASXL1_High", "ASXL1_Intermediate", "ASXL1_Low"
)

labels_combo <- setNames(
  c(
    expression(bold(plain("No-"))*italic("CHIP")),            # No-CHIP_Header
    "   High", "   Intermediate", "   Low",
    expression(bold(italic("DNMT3A"))),                          # DNMT3A_Header
    "   High", "   Intermediate", "   Low",
    expression(bold(plain("non-"))*italic("DNMT3A")),             # non-DNMT3A_Header
    "   High", "   Intermediate", "   Low",
    expression(bold(italic("TET2"))),                            # TET2_Header
    "   High", "   Intermediate", "   Low",
    expression(bold(italic("ASXL1"))),                           # ASXL1_Header
    "   High", "   Intermediate", "   Low"
  ),
  axis_levels
)

# 4.2 Define colors for all three LE8 categories
le8_colors <- c("High" = "#009E73", "Intermediate" = "#E69F00", "Low" = "#56B4E9")

# 4.3 Modified plotting function
plot_forest <- function(df, title_text) {
  # Create header rows (no data points, for labels only)
  header_df <- data.frame(
    combo     = axis_levels[grep("_Header$", axis_levels)],
    estimate  = NA_real_,
    conf.low  = NA_real_,
    conf.high = NA_real_,
    LE8_cat   = factor(NA, levels=c("High","Intermediate","Low")),
    stringsAsFactors = FALSE
  )
  
  # Combine header rows with the actual data
  full_df <- bind_rows(
    header_df,
    df %>% select(combo, estimate, conf.low, conf.high, LE8_cat)
  ) %>%
    mutate(combo = factor(combo, levels = axis_levels))
  
  ggplot(full_df, aes(x = estimate, y = combo, color = LE8_cat)) +
    # Draw points and error bars only for rows with data
    geom_point(
      data = . %>% filter(!is.na(estimate)),
      size = 3, position = position_dodge(width = 0.7)
    ) +
    geom_errorbarh(
      data = . %>% filter(!is.na(estimate)),
      aes(xmin = conf.low, xmax = conf.high),
      height = 0.2, position = position_dodge(width = 0.7)
    ) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
    # Use the expanded y-axis definitions
    scale_y_discrete(
      limits = rev(axis_levels),
      labels = rev(labels_combo[axis_levels])
    ) +
    scale_x_log10() +
    # Use the new color palette for 3 levels
    scale_color_manual(values = le8_colors, na.value="transparent") +
    labs(
      title = title_text,
      x     = "Hazard Ratio (log scale)",
      y     = NULL,
      color = "LE8 Category"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title         = element_text(face = "bold", size = 16, hjust = 0.5),
      axis.text.y        = element_text(hjust = 0),
      axis.text.x        = element_text(size = 12),
      panel.grid.major.x = element_line(color = "gray80"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.position    = "bottom",
      legend.key.width   = grid::unit(1.5, "lines")
    )
}

# ====== 4.4 Call plotting function & save ======

p1 <- plot_forest(
  forest_df %>% filter(outcome == "Intermediate CAD"),
  "Hazard Ratios for CAD"
)
p2 <- plot_forest(
  forest_df %>% filter(outcome == "Composite CVD"),
  "Hazard Ratios for CVD"
)
p3 <- plot_forest(
  forest_df %>% filter(outcome == "CAD Death (HARD)"),
  "Hazard Ratios for CAD Death"
)

p1_no_legend <- p1 + theme(legend.position = "none")
p3_no_legend <- p3 + theme(legend.position = "none")
p2_no_legend <- p2 + theme(legend.position = "none")

# Combine the plots without individual legends
combined_plot_no_legend <- p1_no_legend + p3_no_legend + p2_no_legend +
  plot_layout(ncol = 1) + 
  plot_annotation(tag_levels = '(A)') & 
  theme(plot.tag = element_text(face = 'bold')) 

combined_plot <- (p1 + p3 + p2) + 
  plot_layout(
    ncol = 1,
    guides = 'collect' 
  ) +
  plot_annotation(
    tag_levels = 'A',        
    tag_prefix = '(',        
    tag_suffix = ')'       
  ) &
  theme(
    plot.tag = element_text(face = 'bold'), 
    legend.position = 'bottom',             
    legend.title = element_text(face = "bold"), 
    legend.text = element_text(size = 9)    
  )


print(combined_plot)

# ====== Save the combined plot to disk ======
ggsave(
  "/medpop/esp2/yxliu/CHIP_LE8/Stratified_analysis/Stratified_HR/Plot_uniRef/2HR_Combined_Outcomes.png",
  combined_plot,
  width = 8,
  height = 15 
)




