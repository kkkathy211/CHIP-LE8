############################################################
############ clean version - used for manuscript############
############################################################


# ================= Step 1: Load Data ================= #
library(data.table)
library(dplyr)
library(broom)
library(ggplot2)
library(tableone)

# Load LE8/CVH dataset
le8_data <- fread("/medpop/esp2/yxliu/LE8score_calculation/updated_LE8score.txt")

# Load raw tables
chip_data   <- fread("/medpop/esp2/tnakao/data/UKBB/phenotype/CHIP/450k/2022Oct/CHIP_calls_Oct16_2022.txt")
pheno       <- fread("/medpop/esp2/mzekavat/UKBB/ukbb_PhenoFile.ALL_500k.UpdatedIncdPhenos_202020.REAL.txt") %>% as.data.frame()

# Load the “mask” files
id_450k      <- fread("/medpop/esp2/mesbah/datasets/CHIP/UKBB/450k/ukb450k.eid.list") %>% rename(id = V1)
no_consent   <- fread("/medpop/esp2/projects/UK_Biobank/withdrawn_samples/w7089_20241217.csv") %>% rename(id = V1)
relate_remove<- fread("/medpop/esp2/zyu/chip_protemoics/listofremovefor450K.txt") %>% rename(id = IID)

# Subset pheno to just the 450K participants…
pheno_filt <- pheno %>%
  # keep only those in the 450K ID list
  filter(id %in% id_450k$id) %>%
  # remove withdrawn consent
  filter(!id %in% no_consent$id) %>%
  # only unrelated individuals, drop related ones
  filter(!id %in% relate_remove$id) %>%
  # apply existing QC steps:
  filter(Submitted_Gender == Inferred_Gender,
         Non_Consented     == 0,
         Prev_Myeloid_leukemia  == 0,
         Prev_Lymphoid_leukemia == 0
  ) %>%
  select(id, age, Sex, Sex_numeric, Race, paste0("PC", 1:10))




# ================= Step 2: Process CHIP Data ================= #
# Identify *small* CHIP clones (AF >= 0.02) — **no minAD filter**
CHIP        <- unique(chip_data[chip_data$AF >= 0.02,]$Broad_ID)
DNMT3A      <- unique(chip_data[chip_data$Gene.refGene == "DNMT3A" & chip_data$AF >= 0.02,]$Broad_ID)
nonDNMT3A   <- unique(chip_data[chip_data$AF >= 0.02 & chip_data$Gene.refGene != "DNMT3A",]$Broad_ID)
TET2        <- unique(chip_data[chip_data$Gene.refGene == "TET2"   & chip_data$AF >= 0.02,]$Broad_ID)
ASXL1       <- unique(chip_data[chip_data$Gene.refGene == "ASXL1"  & chip_data$AF >= 0.02,]$Broad_ID)

chip_genes <- c(
  "CHIP", "DNMT3A", "nonDNMT3A", "TET2","ASXL1")

# Merge LE8 dataset with phenotype data
merged_data <- merge(le8_data, pheno_filt, by = "id")

# Add binary indicators for each subtype
for (gene in chip_genes) {
  merged_data[[gene]] <- ifelse(merged_data$id %in% get(gene), 1, 0)
}

# Rename age.x to age if needed
if ("age.x" %in% colnames(merged_data)) {
  merged_data <- merged_data %>% rename(age = age.x)
}




# ================= Step 2.5: Create Composite & Component Scores ================= #
merged_data <- merged_data %>%
  mutate(
    LE8_composite        = (diet_score + BP_Score + PA_score + final_smoking_score +
                              sleep_health_points + bmi_points + Blood_lipid_points + HbA1c_score) / 8,
    LE8_composite_scaled = LE8_composite / 10,
    LE8_composite_std    = as.numeric(scale(LE8_composite)),
    LE8_cat = case_when(
      LE8_composite <  50 ~ "low",
      LE8_composite <  80 ~ "intermediate",
      TRUE                ~ "high"
    )
  ) %>%
  mutate(
    LE8_cat = factor(LE8_cat, levels = c("low", "intermediate", "high")),
    LE8_cat = relevel(LE8_cat, ref = "high")
  ) %>%
  mutate(
    # 10-point scaled
    diet_score_scaled          = diet_score / 10,
    BP_Score_scaled            = BP_Score / 10,
    PA_score_scaled            = PA_score / 10,
    final_smoking_score_scaled = final_smoking_score / 10,
    sleep_health_points_scaled = sleep_health_points / 10,
    bmi_points_scaled          = bmi_points / 10,
    Blood_lipid_points_scaled  = Blood_lipid_points / 10,
    HbA1c_score_scaled         = HbA1c_score / 10,
    # 1-SD standardized
    diet_score_std             = as.numeric(scale(diet_score)),
    BP_Score_std               = as.numeric(scale(BP_Score)),
    PA_score_std               = as.numeric(scale(PA_score)),
    final_smoking_score_std    = as.numeric(scale(final_smoking_score)),
    sleep_health_points_std    = as.numeric(scale(sleep_health_points)),
    bmi_points_std             = as.numeric(scale(bmi_points)),
    Blood_lipid_points_std     = as.numeric(scale(Blood_lipid_points)),
    HbA1c_score_std            = as.numeric(scale(HbA1c_score))
  )


#================================== Basic Info ===================================#
# ================= Generate Cohort Statistics ================= #

# Calculate total number of participants
n_cohort <- nrow(merged_data)

# Calculate mean and SD of age
mean_age <- mean(merged_data$age, na.rm = TRUE)
sd_age <- sd(merged_data$age, na.rm = TRUE)

# Calculate percentage of female participants
n_female <- sum(merged_data$Sex_numeric == 1, na.rm = TRUE)
percent_female <- (n_female / n_cohort) * 100

# Calculate median and IQR of LE8 composite score
median_le8 <- median(merged_data$LE8_composite, na.rm = TRUE)
iqr_le8 <- quantile(merged_data$LE8_composite, probs = c(0.25, 0.75), na.rm = TRUE)

# Print the statistics in the requested format
cat(sprintf(
  "The analytic cohort comprised %d UK Biobank participants (mean [SD] age: %.1f [%.1f] years; %.1f%% female), median (IQR) LE8 score was %.1f (%.1f-%.1f).\n",
  n_cohort,
  mean_age,
  sd_age,
  percent_female,
  median_le8,
  iqr_le8[1],
  iqr_le8[2]
))


racial_group_counts <- table(merged_data$Race)
racial_group_proportions <- prop.table(racial_group_counts) * 100

cat("\nBreakdown by Racial Group:\n")
print(racial_group_counts)
cat("\nPercentage Breakdown by Racial Group:\n")
print(racial_group_proportions)

#============================== Table1 Construction ==============================#

# ── Study Population Table ──
library(tableone)
library(Hmisc)

# tag everyone as “UKBB”
merged_data$dataset <- "UKBB"

# if your data has Sex.x, rename it to Sex
if ("Sex.x" %in% names(merged_data)) {
  merged_data <- merged_data %>% rename(Sex = Sex.x)
}

# ensure age is numeric (continuous)
merged_data$age <- as.numeric(merged_data$age)

# 1) Specify variables and which are factors
tab1_vars       <- c("age", "Sex", "Race", "LE8_cat", "CHIP",
                     "DNMT3A", "nonDNMT3A", "ASXL1", "TET2")
tab1_factorVars <- c("Sex", "Race", "LE8_cat", "CHIP",
                     "DNMT3A", "nonDNMT3A", "ASXL1", "TET2")

# 2) Attach pretty labels via Hmisc
label(merged_data$age)       <- "Age (yr)"
label(merged_data$Sex)       <- "Sex"
label(merged_data$Race)      <- "Race"
label(merged_data$LE8_cat)   <- "LE8 Score"
label(merged_data$CHIP)      <- "CHIP"
label(merged_data$DNMT3A)    <- "DNMT3A"
label(merged_data$nonDNMT3A) <- "non-DNMT3A"
label(merged_data$TET2)      <- "TET2"
label(merged_data$ASXL1)     <- "ASXL1"

# 3) Build the TableOne object (no strata → one column)
table1_obj <- CreateTableOne(
  vars       = tab1_vars,
  data       = merged_data,
  factorVars = tab1_factorVars,
  includeNA  = FALSE
)

# 4) Print to matrix, continuous=mean (SD), categorical=n (%)
table1_mat <- print(
  table1_obj,
  showAllLevels = TRUE,
  contDigits    = 1,
  catDigits     = 1,
  varLabels     = TRUE,
  printToggle   = FALSE
)

# 5) Drop the “n” row so sample size appears only in header
if ("n" %in% rownames(table1_mat)) {
  table1_mat <- table1_mat[rownames(table1_mat) != "n", , drop = FALSE]
}

# 6) Rename the single column to include UKBB sample size
N <- nrow(merged_data)
colnames(table1_mat)[2] <- paste0("UKBB (N=", nrow(merged_data), ")")

# 7) Render with kableExtra, add a footnote and proper caption
library(knitr)
library(kableExtra)
library(tibble)
knitr::kable(
  table1_mat,
  row.names  = TRUE,                    
  col.names  = c("Metrics¹",            
                 paste0("UKBB (N=",nrow(merged_data),")")),
  caption    = "Table 1. Characteristics of the study population in UKBB",
  booktabs   = TRUE
) %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped","hover")) %>%
  footnote(
    general       = "Metrics are presented as mean (SD) for continuous variables and n (%) for categorical variables.",
    general_title = ""
  )


### Table 1 Other information listed here ######

merged_data2 <- merged_data %>%
  mutate(
    CHIP_any = if_else(
      DNMT3A == 1 | nonDNMT3A == 1 | TET2 == 1 | ASXL1 == 1,
      1, 0
    ),
    CHIP_group = factor(
      CHIP_any,
      levels = c(0, 1),
      labels = c("Non_CHIP", "CHIP")
    )
  )

merged_data2 %>%
  count(CHIP_group) %>%
  mutate(pct = n / sum(n) * 100) %>%
  print()
# CHIP_group      n       pct
# 1:   Non_CHIP 399596 93.214582
# 2:       CHIP  29088  6.785418


merged_data2 %>%
  group_by(CHIP_group) %>%
  summarise(
    n        = n(),
    mean_age = mean(age, na.rm = TRUE),
    sd_age   = sd(age,   na.rm = TRUE)
  ) %>%
  print()

#CHIP_group      n mean_age sd_age
#<fct>       <int>    <dbl>  <dbl>
# 1 Non_CHIP   399596     56.3   8.08
# 2 CHIP        29088     59.7   7.18


merged_data2 %>%
  group_by(CHIP_group, Sex) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(CHIP_group) %>%
  mutate(pct = n / sum(n) * 100) %>%
  print()
# A tibble: 4 x 4
# Groups:   CHIP_group [2]
# CHIP_group Sex             n   pct
#<fct>      <labelled>  <int> <dbl>
#  1 Non_CHIP   Female     216485  54.2
#  2 Non_CHIP   Male       183111  45.8
#  3 CHIP       Female      15614  53.7
#  4 CHIP       Male        13474  46.3


merged_data2 %>%
  group_by(CHIP_group, Race) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(CHIP_group) %>%
  mutate(pct = n / sum(n) * 100) %>%
  print()
# A tibble: 10 x 4
# Groups:   CHIP_group [2]
#CHIP_group Race             n    pct
#<fct>      <labelled>   <int>  <dbl>
#1 Non_CHIP   Black         6185  1.55 
#2 Non_CHIP   Chinese       1246  0.312
#3 Non_CHIP   Other         7170  1.79 
#4 Non_CHIP   South Asian   7752  1.94 
#5 Non_CHIP   White       377243 94.4  
#6 CHIP       Black          329  1.13 
#7 CHIP       Chinese         83  0.285
#8 CHIP       Other          462  1.59 
#9 CHIP       South Asian    456  1.57 
#10 CHIP       White        27758 95.4  


merged_data2 %>%
  group_by(CHIP_group, LE8_cat) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(CHIP_group) %>%
  mutate(pct = n / sum(n) * 100) %>%
  print()

# A tibble: 6 x 4
# Groups:   CHIP_group [2]
#CHIP_group LE8_cat           n   pct
#<fct>      <fct>         <int> <dbl>
#1 Non_CHIP   high          33297  8.33
#2 Non_CHIP   low           43749 10.9 
#3 Non_CHIP   intermediate 322550 80.7 
#4 CHIP       high           1948  6.70
#5 CHIP       low            3492 12.0 
#6 CHIP       intermediate  23648 81.3 


le8_orig_items <- c(
  "diet_score",
  "PA_score",
  "final_smoking_score",
  "sleep_health_points",
  "bmi_points",
  "Blood_lipid_points",
  "HbA1c_score",
  "BP_Score"
)


merged_data2 %>%
  group_by(CHIP_group) %>%
  summarise(across(
    all_of(le8_orig_items),
    list(
      mean = ~ mean(.x, na.rm = TRUE),
      sd   = ~ sd(  .x, na.rm = TRUE)
    ),
    .names = "{col}_{fn}"
  )) %>%
  print(n = Inf)

le8_item_stats <- merged_data2 %>%
  group_by(CHIP_group) %>%
  summarise(across(
    all_of(le8_orig_items),
    list(
      mean = ~ mean(.x, na.rm = TRUE),
      sd   = ~ sd(  .x, na.rm = TRUE)
    ),
    .names = "{col}_{fn}"
  ))

as.data.frame(le8_item_stats)

#CHIP_group diet_score_mean diet_score_sd PA_score_mean PA_score_sd
#1   Non_CHIP        52.17391      32.26805      32.16559    28.08441
#2       CHIP        51.91093      32.25835      32.32433    28.63277
#final_smoking_score_mean final_smoking_score_sd sleep_health_points_mean
#1                 82.79312               31.01506                 89.65072
#2                 80.85809               32.25372                 89.58230
#sleep_health_points_sd bmi_points_mean bmi_points_sd Blood_lipid_points_mean
#1               19.06398        68.69222      28.57914                47.37780
#2               19.14836        68.12053      28.41161                47.01183
#Blood_lipid_points_sd HbA1c_score_mean HbA1c_score_sd BP_Score_mean
#1              30.23090         90.27250       21.45207      49.00866
#2              29.91851         88.74209       22.64544      45.68980
#BP_Score_sd
#1    27.44843
#2    26.69801




#### Add p-value
# ===================================================================
# Final Table 1 Construction with P-values (using tableone)
# ===================================================================
library(tableone)
library(dplyr)

# 1. Define the grouping variable
merged_data2 <- merged_data %>%
  mutate(
    CHIP_any = if_else(
      DNMT3A == 1 | nonDNMT3A == 1 | TET2 == 1 | ASXL1 == 1, 1, 0
    ),
    CHIP_group = factor(
      CHIP_any,
      levels = c(0, 1),
      labels = c("No CHIP", "CHIP")
    )
  )

# 2. Define the list of all variables for the table
le8_orig_items <- c(
  "diet_score", "PA_score", "final_smoking_score", "sleep_health_points",
  "bmi_points", "Blood_lipid_points", "HbA1c_score", "BP_Score"
)

table_variables <- c(
  "age", "Sex", "Race", "LE8_cat", le8_orig_items
)

# 3. Specify which of these variables are categorical (factors)
factor_variables <- c("Sex", "Race", "LE8_cat")

# 4. Create the TableOne object, stratified by CHIP status
table1_stratified <- CreateTableOne(
  vars = table_variables,
  strata = "CHIP_group", 
  data = merged_data2,
  factorVars = factor_variables,
  addOverall = FALSE 
)

# 5. Print the table to the console to see the results.
cat("--- Table 1: Baseline Characteristics by CHIP Status ---\n")
print(
  table1_stratified,
  nonnormal = c("age", le8_orig_items), 
  showAllLevels = TRUE,
  contDigits = 1,
  catDigits = 1,
  printToggle = TRUE, 
  explain = FALSE,
  varLabels = TRUE 
)

# 6. Convert the table to a data frame for saving or custom formatting
table1_df <- print(
  table1_stratified,
  nonnormal = c("age", le8_orig_items),
  showAllLevels = TRUE,
  contDigits = 1,
  catDigits = 1,
  printToggle = FALSE, 
  explain = FALSE,
  varLabels = TRUE
)

write.csv(table1_df, file = "/medpop/esp2/yxliu/CHIP_LE8/Cross_sectional_analysis/Table1_results.csv")


#============================ Table 1 CHIP & subtypes ================================#

# Define the CHIP subtype columns
chip_subtypes_to_analyze <- c("CHIP", "DNMT3A", "nonDNMT3A", "TET2", "ASXL1")

# Get the total number of participants in the merged_data
total_participants <- nrow(merged_data)

# Create an empty data frame to store results
chip_counts_df <- data.frame(
  Subtype = character(),
  Count = numeric(),
  Percentage = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each CHIP subtype and calculate count and percentage
for (subtype in chip_subtypes_to_analyze) {
  # Count participants with this CHIP subtype (where the column value is 1)
  count <- sum(merged_data[[subtype]] == 1, na.rm = TRUE)
  
  # Calculate the percentage
  percentage <- (count / total_participants) * 100
  
  # Add to the results data frame
  chip_counts_df <- rbind(chip_counts_df, data.frame(
    Subtype = subtype,
    Count = count,
    Percentage = percentage
  ))
}

cat("Number and Percentage of Patients by CHIP Subtype:\n")
print(chip_counts_df)


# ================= Step 3: Create Outcome Variables ================= #
data_CHIP      <- merged_data %>% filter(CHIP == 0 | CHIP == 1) %>%    mutate(CHIP_outcome      = ifelse(CHIP == 1, 1, 0))
data_DNMT3A    <- merged_data %>% filter(CHIP == 0 | DNMT3A == 1) %>%  mutate(DNMT3A_outcome    = ifelse(DNMT3A    == 1, 1, 0))
data_nonDNMT3A <- merged_data %>% filter(CHIP == 0 | nonDNMT3A == 1) %>% mutate(nonDNMT3A_outcome = ifelse(nonDNMT3A == 1, 1, 0))
data_TET2     <- merged_data %>% filter(CHIP == 0 | TET2      == 1) %>% mutate(TET2_outcome      = ifelse(TET2      == 1, 1, 0))
data_ASXL1    <- merged_data %>% filter(CHIP == 0 | ASXL1     == 1) %>% mutate(ASXL1_outcome     = ifelse(ASXL1     == 1, 1, 0))


# ================= Step 4: Logistic Regression ================= #
covariates    <- c("age","Sex_numeric", paste0("PC",1:10))
le8_variables <- c(
  "diet_score_scaled","diet_score_std",
  "BP_Score_scaled","BP_Score_std",
  "PA_score_scaled","PA_score_std",
  "final_smoking_score_scaled","final_smoking_score_std",
  "sleep_health_points_scaled","sleep_health_points_std",
  "bmi_points_scaled","bmi_points_std",
  "Blood_lipid_points_scaled","Blood_lipid_points_std",
  "HbA1c_score_scaled","HbA1c_score_std",
  "LE8_composite_scaled","LE8_composite_std","LE8_cat"
)
comparisons <- list(
  CHIP      = list(data=data_CHIP,      outcome="CHIP_outcome"),
  DNMT3A    = list(data=data_DNMT3A,    outcome="DNMT3A_outcome"),
  nonDNMT3A = list(data=data_nonDNMT3A, outcome="nonDNMT3A_outcome"),
  TET2      = list(data=data_TET2,      outcome="TET2_outcome"),
  ASXL1     = list(data=data_ASXL1,     outcome="ASXL1_outcome")
)
results_list <- list()

for(comp in names(comparisons)){
  comp_data   <- comparisons[[comp]]$data
  outcome_var <- comparisons[[comp]]$outcome
  for(var in le8_variables){
    mdl <- glm(as.formula(paste(outcome_var,"~",var,"+",paste(covariates,collapse="+"))),
               data=comp_data, family=binomial(link="logit"))
    res <- broom::tidy(mdl)
    res <- if(var=="LE8_cat") res %>% filter(grepl("^LE8_cat",term)) else res %>% filter(term==var)
    res$comparison <- comp; res$LE8_Component<-var
    results_list[[paste(comp,var,sep="_")]]<-res
  }
}
results_df <- bind_rows(results_list) %>%
  mutate(
    Odds_Ratio  = exp(estimate),
    Lower_CI    = exp(estimate-1.96*std.error),
    Upper_CI    = exp(estimate+1.96*std.error),
    Significance= ifelse(p.value<0.05,"Significant","Not Significant")
  )
write.csv(results_df,"/medpop/esp2/yxliu/CHIP_LE8/Cross_sectional_analysis/CHIP_LE8_association_results_VAF2.csv",row.names=FALSE)



############## Sequence: CHIP -> DNMT3A -> non-DNMT3A -> TET2 -> ASXL1 ####################

comp_levels <- c("CHIP","DNMT3A","nonDNMT3A","TET2","ASXL1")


# ================= Step 5: Visualization ================= #

library(ggplot2)
library(dplyr)
library(patchwork) 
library(scales) 

# ================= 1: Create Plot A (Continuous Score) ================= #

# --- a) Prepare data for Plot A ---
data_std <- results_df %>% 
  filter(LE8_Component == "LE8_composite_std") %>%
  mutate(comparison = factor(comparison, levels = comp_levels))

# --- b) Calculate optimal y-axis limits for Plot A ---
plot_a_min <- min(data_std$Lower_CI, na.rm = TRUE) * 0.95 
plot_a_max <- max(data_std$Upper_CI, na.rm = TRUE) * 1.05 

# --- c) Build Plot A ---
p_std <- ggplot(data_std, aes(x = comparison, y = Odds_Ratio, color = comparison)) +
  geom_pointrange(aes(ymin = Lower_CI, ymax = Upper_CI), size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  scale_y_continuous(limits = c(plot_a_min, plot_a_max), breaks = pretty_breaks(n = 5)) +
  scale_x_discrete(labels = italic_labels) +
  scale_color_manual(values = RColorBrewer::brewer.pal(5, "Set1")) +
  labs(
    x = "CHIP Subtype",
    y = "Odds Ratio (95% CI)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

# ================= 2: Create Plot B (Categorical Score) ================= #

# --- a) Prepare data for Plot B ---
data_cat <- results_df %>% 
  filter(LE8_Component == "LE8_cat") %>%
  mutate(
    comparison = factor(comparison, levels = comp_levels),
    category = factor(
      gsub("LE8_cat", "", term),
      levels = c("low", "intermediate"),
      labels = c("Low (<50)", "Intermediate (50-80)")
    )
  )

# --- b) Calculate optimal y-axis limits for Plot B ---
plot_b_min <- min(data_cat$Lower_CI, na.rm = TRUE) * 0.95
plot_b_max <- max(data_cat$Upper_CI, na.rm = TRUE) * 1.05

# --- c) Build Plot B ---
p_cat <- ggplot(data_cat, aes(x = comparison, y = Odds_Ratio, color = category)) +
  geom_pointrange(
    aes(ymin = Lower_CI, ymax = Upper_CI),
    size = 1.2,
    position = position_dodge(width = 0.7)
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  scale_y_continuous(limits = c(plot_b_min, plot_b_max), breaks = pretty_breaks(n = 5)) +
  scale_x_discrete(labels = italic_labels) +
  scale_color_brewer(
    palette = "Set2",
    name = "LE8 Category\n(vs. High)"
  ) +
  labs(
    x = "CHIP Subtype",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.title.y = element_blank(),
    legend.position = "right"
  )

# ================= 3: Combine Plots with Panel Labels ================= #

combined_plot <- (p_std | p_cat) +
  plot_annotation(
    title = "Association of Life's Essential 8 Score with CHIP Prevalence",
    subtitle = "Odds Ratios from logistic regression models",
    tag_levels = 'A',      
    tag_prefix = '(',      
    tag_suffix = ')'       
  ) &
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 16)
  )

# Display the combined plot
print(combined_plot)

# ================= 4: Save the Combined Plot ================= #
ggsave(
  filename = "/medpop/esp2/yxliu/CHIP_LE8/Cross_sectional_analysis/Plots/VAF2_LE8_composite_combined_final.png",
  plot     = combined_plot,
  width    = 14,
  height   = 7,
  dpi      = 300,
  bg       = "white"
)


# ================= Step 5B: Individual LE8 Component 1-SD Forest Plots ================= #

library(ggplot2)
library(dplyr)
library(patchwork) 
library(rlang)     

# ================= 1: Define Constants and Mappings ================= #

component_map_std <- c(
  Blood_lipid_points_std  = "Blood Lipid",
  bmi_points_std          = "BMI",
  sleep_health_points_std = "Sleep Health",
  HbA1c_score_std         = "Glucose",
  PA_score_std            = "Physical Activity",
  diet_score_std          = "Diet",
  BP_Score_std            = "Blood Pressure",
  final_smoking_score_std = "Smoking"
)

component_order <- results_df %>%
  filter(
    comparison   == "CHIP",
    grepl("_std$", LE8_Component),
    !is.na(Odds_Ratio)
  ) %>%
  mutate(Component = component_map_std[LE8_Component]) %>%
  arrange(Odds_Ratio) %>%
  pull(Component) %>%
  unique()

subtypes <- c("CHIP", "DNMT3A", "nonDNMT3A", "TET2", "ASXL1")

# ================= 2: Create All 5 Plots in a Loop ================= #

plot_list <- list()

for (sub in subtypes) {
  df_plot <- results_df %>%
    filter(
      comparison == sub,
      grepl("_std$", LE8_Component),
      !is.na(Lower_CI), !is.na(Upper_CI)
    ) %>%
    mutate(Component = component_map_std[LE8_Component]) %>%
    filter(!is.na(Component)) %>%                     
    mutate(
      # Apply the consistent ordering to the y-axis
      Component = factor(Component, levels = component_order)  
    )
  
  if (nrow(df_plot) == 0) {
    plot_list[[sub]] <- ggplot() + theme_void()
    next
  }
  
  # 2. Create the correctly formatted title for each subplot
  title_expr <- if (grepl("^non", sub)) {
    gene <- sub("non", "", sub)
    bquote("non-"*italic(.(gene)))
  } else {
    bquote(italic(.(sub)))
  }
  
  # 3. Create the plot object with all labels and improved aesthetics
  p <- ggplot(df_plot, aes(x = Odds_Ratio, y = Component)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI),
                   height = 0.2,
                   color = "#0072B2",
                   linewidth = 0.8) +
    geom_point(color = "#0072B2", size = 3) +
    coord_cartesian(xlim = c(0.7, 1.08)) + 
    labs(
      title = title_expr,
      x     = "Odds Ratio (95% CI)",
      y     = NULL 
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.y = element_text(size = 12),
      legend.position = "none"
    )
  
  # 4. Store the plot in the list
  plot_list[[sub]] <- p
}

# ================= 3: Arrange the Plots with Patchwork ================= #

# Arrange the 5 plots into a 3-column grid. 
combined_plot <- (plot_list$CHIP) + (plot_list$DNMT3A) + (plot_list$nonDNMT3A) +
  (plot_list$TET2) + (plot_list$ASXL1)  + plot_spacer() +
  plot_layout(ncol = 3)

# ================= 4: Add Final Overall Annotation ================= #

final_figure <- combined_plot +
  plot_annotation(
    title = "Association of Individual LE8 Components with CHIP Subtypes",
    subtitle = "Odds Ratios per 1-SD Increase in Each Component Score"
  ) &
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 16, margin = margin(b=10))
  )

print(final_figure)

# ================= 5: Save the Combined Figure ================= #

outdir <- "/medpop/esp2/yxliu/CHIP_LE8/Cross_sectional_analysis/Plots/"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

ggsave(
  filename = file.path(outdir, "forestplot_LE8_components_final.png"),
  plot     = final_figure,
  width    = 18, 
  height   = 9,  
  dpi      = 300,
  bg       = "white"
)





##########################################################################
#============================Supplementary Figures========================
##########################################################################


# Larger CHIP VAF10

# ================= Step 2: Process CHIP Data ================= #
# Identify *larger* CHIP clones (AF >= 0.1) — **no minAD filter**
CHIP        <- unique(chip_data[chip_data$AF >= 0.1,]$Broad_ID)
DNMT3A      <- unique(chip_data[chip_data$Gene.refGene == "DNMT3A" & chip_data$AF >= 0.1,]$Broad_ID)
nonDNMT3A   <- unique(chip_data[chip_data$AF >= 0.1 & chip_data$Gene.refGene != "DNMT3A",]$Broad_ID)
TET2        <- unique(chip_data[chip_data$Gene.refGene == "TET2"   & chip_data$AF >= 0.1,]$Broad_ID)
ASXL1       <- unique(chip_data[chip_data$Gene.refGene == "ASXL1"  & chip_data$AF >= 0.1,]$Broad_ID)

chip_genes <- c(
  "CHIP", "DNMT3A", "nonDNMT3A", "TET2","ASXL1")

# Merge LE8 dataset with phenotype data
merged_data <- merge(le8_data, pheno_filt, by = "id")

# Add binary indicators for each subtype
for (gene in chip_genes) {
  merged_data[[gene]] <- ifelse(merged_data$id %in% get(gene), 1, 0)
}

# Rename age.x to age if needed
if ("age.x" %in% colnames(merged_data)) {
  merged_data <- merged_data %>% rename(age = age.x)
}


# ================= Step 2.5: Create Composite & Component Scores ================= #
merged_data <- merged_data %>%
  mutate(
    LE8_composite        = (diet_score + BP_Score + PA_score + final_smoking_score +
                              sleep_health_points + bmi_points + Blood_lipid_points + HbA1c_score) / 8,
    LE8_composite_scaled = LE8_composite / 10,
    LE8_composite_std    = as.numeric(scale(LE8_composite)),
    LE8_cat = case_when(
      LE8_composite <  50 ~ "low",
      LE8_composite <  80 ~ "intermediate",
      TRUE                ~ "high"
    )
  ) %>%
  mutate(
    LE8_cat = factor(LE8_cat, levels = c("low", "intermediate", "high")),
    LE8_cat = relevel(LE8_cat, ref = "high")
  ) %>%
  mutate(
    # 10-point scaled
    diet_score_scaled          = diet_score / 10,
    BP_Score_scaled            = BP_Score / 10,
    PA_score_scaled            = PA_score / 10,
    final_smoking_score_scaled = final_smoking_score / 10,
    sleep_health_points_scaled = sleep_health_points / 10,
    bmi_points_scaled          = bmi_points / 10,
    Blood_lipid_points_scaled  = Blood_lipid_points / 10,
    HbA1c_score_scaled         = HbA1c_score / 10,
    # 1-SD standardized
    diet_score_std             = as.numeric(scale(diet_score)),
    BP_Score_std               = as.numeric(scale(BP_Score)),
    PA_score_std               = as.numeric(scale(PA_score)),
    final_smoking_score_std    = as.numeric(scale(final_smoking_score)),
    sleep_health_points_std    = as.numeric(scale(sleep_health_points)),
    bmi_points_std             = as.numeric(scale(bmi_points)),
    Blood_lipid_points_std     = as.numeric(scale(Blood_lipid_points)),
    HbA1c_score_std            = as.numeric(scale(HbA1c_score))
  )


# ================= Step 3: Create Outcome Variables ================= #
data_CHIP      <- merged_data %>% filter(CHIP == 0 | CHIP == 1) %>%    mutate(CHIP_outcome      = ifelse(CHIP == 1, 1, 0))
data_DNMT3A    <- merged_data %>% filter(CHIP == 0 | DNMT3A == 1) %>%  mutate(DNMT3A_outcome    = ifelse(DNMT3A    == 1, 1, 0))
data_nonDNMT3A <- merged_data %>% filter(CHIP == 0 | nonDNMT3A == 1) %>% mutate(nonDNMT3A_outcome = ifelse(nonDNMT3A == 1, 1, 0))
data_TET2     <- merged_data %>% filter(CHIP == 0 | TET2      == 1) %>% mutate(TET2_outcome      = ifelse(TET2      == 1, 1, 0))
data_ASXL1    <- merged_data %>% filter(CHIP == 0 | ASXL1     == 1) %>% mutate(ASXL1_outcome     = ifelse(ASXL1     == 1, 1, 0))


# ================= Step 4: Logistic Regression ================= #
covariates    <- c("age","Sex_numeric", paste0("PC",1:10))
le8_variables <- c(
  "diet_score_scaled","diet_score_std",
  "BP_Score_scaled","BP_Score_std",
  "PA_score_scaled","PA_score_std",
  "final_smoking_score_scaled","final_smoking_score_std",
  "sleep_health_points_scaled","sleep_health_points_std",
  "bmi_points_scaled","bmi_points_std",
  "Blood_lipid_points_scaled","Blood_lipid_points_std",
  "HbA1c_score_scaled","HbA1c_score_std",
  "LE8_composite_scaled","LE8_composite_std","LE8_cat"
)
comparisons <- list(
  CHIP      = list(data=data_CHIP,      outcome="CHIP_outcome"),
  DNMT3A    = list(data=data_DNMT3A,    outcome="DNMT3A_outcome"),
  nonDNMT3A = list(data=data_nonDNMT3A, outcome="nonDNMT3A_outcome"),
  TET2      = list(data=data_TET2,      outcome="TET2_outcome"),
  ASXL1     = list(data=data_ASXL1,     outcome="ASXL1_outcome")
)
results_list <- list()

for(comp in names(comparisons)){
  comp_data   <- comparisons[[comp]]$data
  outcome_var <- comparisons[[comp]]$outcome
  for(var in le8_variables){
    mdl <- glm(as.formula(paste(outcome_var,"~",var,"+",paste(covariates,collapse="+"))),
               data=comp_data, family=binomial(link="logit"))
    res <- broom::tidy(mdl)
    res <- if(var=="LE8_cat") res %>% filter(grepl("^LE8_cat",term)) else res %>% filter(term==var)
    res$comparison <- comp; res$LE8_Component<-var
    results_list[[paste(comp,var,sep="_")]]<-res
  }
}
results_df <- bind_rows(results_list) %>%
  mutate(
    Odds_Ratio  = exp(estimate),
    Lower_CI    = exp(estimate-1.96*std.error),
    Upper_CI    = exp(estimate+1.96*std.error),
    Significance= ifelse(p.value<0.05,"Significant","Not Significant")
  )
write.csv(results_df,"/medpop/esp2/yxliu/CHIP_LE8/Cross_sectional_analysis/CHIP_LE8_association_results_VAF10.csv",row.names=FALSE)



############## Sequence: CHIP -> DNMT3A -> non-DNMT3A -> TET2 -> ASXL1 ####################

comp_levels <- c("CHIP","DNMT3A","nonDNMT3A","TET2","ASXL1")


# ================= Step 5: Visualization ================= #

library(ggplot2)
library(dplyr)
library(patchwork) 
library(scales) 

# ================= 1: Create Plot A (Continuous Score) ================= #

# --- a) Prepare data for Plot A ---
data_std <- results_df %>% 
  filter(LE8_Component == "LE8_composite_std") %>%
  mutate(comparison = factor(comparison, levels = comp_levels))

# --- b) Calculate optimal y-axis limits for Plot A ---
plot_a_min <- min(data_std$Lower_CI, na.rm = TRUE) * 0.95 
plot_a_max <- max(data_std$Upper_CI, na.rm = TRUE) * 1.05 

# --- c) Build Plot A ---
p_std <- ggplot(data_std, aes(x = comparison, y = Odds_Ratio, color = comparison)) +
  geom_pointrange(aes(ymin = Lower_CI, ymax = Upper_CI), size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  scale_y_continuous(limits = c(plot_a_min, plot_a_max), breaks = pretty_breaks(n = 5)) +
  scale_x_discrete(labels = italic_labels) +
  scale_color_manual(values = RColorBrewer::brewer.pal(5, "Set1")) +
  labs(
    x = "CHIP Subtype",
    y = "Odds Ratio (95% CI)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

# ================= 2: Create Plot B (Categorical Score) ================= #

# --- a) Prepare data for Plot B ---
data_cat <- results_df %>% 
  filter(LE8_Component == "LE8_cat") %>%
  mutate(
    comparison = factor(comparison, levels = comp_levels),
    category = factor(
      gsub("LE8_cat", "", term),
      levels = c("low", "intermediate"),
      labels = c("Low (<50)", "Intermediate (50-80)")
    )
  )

# --- b) Calculate optimal y-axis limits for Plot B ---
plot_b_min <- min(data_cat$Lower_CI, na.rm = TRUE) * 0.95
plot_b_max <- max(data_cat$Upper_CI, na.rm = TRUE) * 1.05

# --- c) Build Plot B ---
p_cat <- ggplot(data_cat, aes(x = comparison, y = Odds_Ratio, color = category)) +
  geom_pointrange(
    aes(ymin = Lower_CI, ymax = Upper_CI),
    size = 1.2,
    position = position_dodge(width = 0.7)
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  scale_y_continuous(limits = c(plot_b_min, plot_b_max), breaks = pretty_breaks(n = 5)) +
  scale_x_discrete(labels = italic_labels) +
  scale_color_brewer(
    palette = "Set2",
    name = "LE8 Category\n(vs. High)"
  ) +
  labs(
    x = "CHIP Subtype",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.title.y = element_blank(),
    legend.position = "right"
  )

# ================= 3: Combine Plots with Panel Labels ================= #

combined_plot <- (p_std | p_cat) +
  plot_annotation(
    title = "Association of Life's Essential 8 Score with CHIP Prevalence",
    subtitle = "Odds Ratios from logistic regression models",
    tag_levels = 'A',      
    tag_prefix = '(',      
    tag_suffix = ')'       
  ) &
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 16)
  )

# Display the combined plot
print(combined_plot)

# ================= 4: Save the Combined Plot ================= #
ggsave(
  filename = "/medpop/esp2/yxliu/CHIP_LE8/Cross_sectional_analysis/Plots/VAF10_LE8_composite_combined_final.png",
  plot     = combined_plot,
  width    = 14,
  height   = 7,
  dpi      = 300,
  bg       = "white"
)

