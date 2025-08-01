library(dplyr)
library(tidyverse)
library(data.table)
library(ukbtools) 
library(survival)
library(survminer) 
library(ggrepel) 
library(gridExtra) 
library(plinkQC) 
library(tableone) 
library(ggplot2)
library(mice)

setwd("/medpop/esp2/")

### Read in Data 
pheno = fread("/medpop/esp2/mzekavat/UKBB/ukbb_PhenoFile.ALL_500k.UpdatedIncdPhenos_202020.REAL.txt")
pheno = data.frame(pheno)
pheno$IDs_toRemove_SampleQCv3 = ifelse((!(pheno$Submitted_Gender == pheno$Inferred_Gender) |  pheno$Non_Consented== 1),1,0) ## N:
pheno = pheno[which(pheno$IDs_toRemove_SampleQCv3==0),]

mypheno <- pheno[, c(1,12,13,14, # basic info
                     2273, #DASH diet score
                     2277,2278,2279,2280, #PA
                     1432, #smoke
                     2305, #sleep
                     18, #bmi
                     178,185,2174, #blood lipid score
                     2238,184, #glucose score
                     22,23,24)] #blood pressure score

# Initial missing summary
naniar::miss_var_summary(mypheno)

#### missing data
mypheno <- mypheno %>%
  mutate(
    # if either is answered, impute 0 for the missing one.
    days_10m_moderate_activity = ifelse(
      !is.na(days_10m_moderate_activity) | !is.na(days_10m_vigorous_activity),
      ifelse(is.na(days_10m_moderate_activity), 0, days_10m_moderate_activity),
      NA_real_
    ),
    days_10m_vigorous_activity = ifelse(
      !is.na(days_10m_moderate_activity) | !is.na(days_10m_vigorous_activity),
      ifelse(is.na(days_10m_vigorous_activity), 0, days_10m_vigorous_activity),
      NA_real_
    ),
    # if either days question is answered, impute 0 for missing duration.
    moderate_activity_duration = ifelse(
      !is.na(days_10m_moderate_activity) | !is.na(days_10m_vigorous_activity),
      ifelse(is.na(moderate_activity_duration), 0, moderate_activity_duration),
      NA_real_
    ),
    vigorous_activity_duration = ifelse(
      !is.na(days_10m_moderate_activity) | !is.na(days_10m_vigorous_activity),
      ifelse(is.na(vigorous_activity_duration), 0, vigorous_activity_duration),
      NA_real_
    )
  )

# missing summary
naniar::miss_var_summary(mypheno)

# Step 1: Impute categorical variable
mypheno$cholesterol_lowering_medication[is.na(mypheno$cholesterol_lowering_medication)] <- "no"

# Step 2: Define continuous variables to impute
vars_to_impute <- c("HDL.cholesterol", 
                    "SBP", "DBP", 
                    "Glycated.haemoglobin..HbA1c.", 
                    "Cholesterol", 
                    "days_10m_moderate_activity", 
                    "moderate_activity_duration", 
                    "days_10m_vigorous_activity", 
                    "vigorous_activity_duration")

# Step 3: Subset data for imputation (including other covariates if needed)
impute_data <- mypheno[, vars_to_impute]

# Step 4: Run MICE for imputation (default method = pmm for numeric)
mice_result <- mice(impute_data, m = 5, method = "pmm", seed = 123)

# Step 5: Complete the dataset using first imputed set
imputed_complete <- complete(mice_result, 1)

# Step 6: Replace the original variables in mypheno with the imputed ones
mypheno[, vars_to_impute] <- imputed_complete

# Optional: Check for any remaining missingness
sapply(mypheno[, vars_to_impute], function(x) sum(is.na(x)))

# missing summary
naniar::miss_var_summary(mypheno)


colSums(is.na(mypheno))
mypheno <- na.omit(mypheno)


############# LE8 Scores Calculation ############

### 1. DASH Diet Score ###
percentiles <- quantile(mypheno$dietscore_CONTINUOUS, probs = c(0.95, 0.75, 0.50, 0.25), na.rm = TRUE)
mypheno <- mypheno %>%
  mutate(
    diet_score = case_when(
      dietscore_CONTINUOUS >= percentiles[1] ~ 100,  # > 95th
      dietscore_CONTINUOUS >= percentiles[2] ~ 80,   # 75th–94th 
      dietscore_CONTINUOUS >= percentiles[3] ~ 50,   # 50th–74th
      dietscore_CONTINUOUS >= percentiles[4] ~ 25,   # 25th–49th
      TRUE ~ 0                                       # 1st–24th
    )
  )


### 2. Physical Activity Score (denote PA) ###
# Function to calculate PA
calculate_pa_points <- function(PA) {
  case_when(
    PA >= 150 ~ 100,
    PA >= 120 & PA <= 149 ~ 90,
    PA >= 90 & PA <= 119 ~ 80,
    PA >= 60 & PA <= 89 ~ 60,
    PA >= 30 & PA <= 59 ~ 40,
    PA >= 1 & PA <= 29 ~ 20,
    TRUE ~ 0
  )
}

# Apply
mypheno <- mypheno %>%
  mutate(
    total_moderate_minutes = (days_10m_moderate_activity * moderate_activity_duration)/10,
    total_vigorous_minutes = (days_10m_vigorous_activity * vigorous_activity_duration)/10,
    total_PA_minutes = total_moderate_minutes + (2 * total_vigorous_minutes),
    PA_score = calculate_pa_points(total_PA_minutes)
  )


### 3. Tobacco/nicotine Exposure Score ###
## Data Preprocess
# For Smoking Status
smoke_sentence_counts <- table(pheno$smok_detailed_)
# Convert to a data frame
smoke_sentence_counts_df <- as.data.frame(smoke_sentence_counts)

# Define the complete mapping
smok_detailed_mapping <- c(
  "current cigar pipe smoker, former cigarette smoker" = 0,
  "current cigar pipe smoker, not former cigarette smoker" = 0,
  "current cigarette smoker, <10/day" = 0,
  "current cigarette smoker, ≥40/day" = 0,
  "current cigarette smoker, 10 to <20/day" = 0,
  "current cigarette smoker, 20 to <40/day" = 0,
  "current occasional smoker, smoked <100 cigarettes in lifetime" = 0,
  "current occasional smoker, smoked ≥100 cigarettes in lifetime" = 0,
  "current occasional smoker, smoked cigarettes daily in past, <20/day" = 0,
  "current occasional smoker, smoked cigarettes daily in past, ≥20/day" = 0,
  "current occasional smoker, smoked cigars or pipes daily in past" = 0,
  "former cigarette smoker, <20/day, quit <1 year ago" = 25,
  "former cigarette smoker, <20/day, quit ≥20 year ago" = 75,
  "former cigarette smoker, <20/day, quit 1–5 year ago" = 50,
  "former cigarette smoker, <20/day, quit 10–20 year ago" = 75,
  "former cigarette smoker, <20/day, quit 5–10 year ago" = 75,
  "former cigarette smoker, ≥20/day, quit <1 year ago" = 25,
  "former cigarette smoker, ≥20/day, quit ≥20 year ago" = 75,
  "former cigarette smoker, ≥20/day, quit 1–5 year ago" = 50,
  "former cigarette smoker, ≥20/day, quit 10–20 year ago" = 75,
  "former cigarette smoker, ≥20/day, quit 5–10 year ago" = 75,
  "former daily cigar pipe smoker" = 75,
  "former occasional cigarette smoker, lifetime cigarette smoking unknown" = 50,
  "former occasional cigarette smoker, smoked <100 cigarettes in lifetime" = 75,
  "former occasional cigarette smoker, smoked ≥100 cigarettes in lifetime" = 50,
  "missing" = 100,
  "never smoker" = 100
)

# Apply
mypheno <- mypheno %>%
  mutate(
    smoking_score_detailed = case_when(
      smok_detailed_ %in% names(smok_detailed_mapping) ~ smok_detailed_mapping[smok_detailed_],
      TRUE ~ 100
    ),
    final_smoking_score = smoking_score_detailed)


### 4. Sleep Health Score ###
# Function to calculate *Sleep_health_points* 
sleep_health_points <- function(sleep_health) {
  if (sleep_health >= 7 && sleep_health <= 9) {
    sleep_health_points <- 100
  } else if (sleep_health >= 9 && sleep_health < 10) {
    sleep_health_points <- 90
  } else if (sleep_health >= 6 && sleep_health < 7) {
    sleep_health_points <- 70
  } else if ((sleep_health >= 5 && sleep_health < 6) || (sleep_health >= 10)) {
    sleep_health_points <- 40
  } else if (sleep_health >= 4 && sleep_health < 5) {
    sleep_health_points <- 20
  } else {
    sleep_health_points <- 0
  }
  
  return(sleep_health_points)
}

# Apply (store at col 24)
mypheno <- mypheno %>%
  mutate(sleep_health_points = sapply(sleepduration, sleep_health_points))


### 5. Body Mass Index Score ###
# Function to calculate *bmi* 
bmi_points <- function(bmi) {
  if (bmi < 25) {
    bmi_points <- 100
  } else if (bmi >= 25 && bmi <= 29.9) {
    bmi_points <- 70
  } else if (bmi >= 30 && bmi <= 34.9) {
    bmi_points <- 30
  } else if (bmi >= 35 && bmi <= 39.9) {
    bmi_points <- 15
  } else {
    bmi_points <- 0
  }
  
  return(bmi_points)
}

# Apply (store at col 25)
mypheno <- mypheno %>%
  mutate(bmi_points = sapply(BMI, bmi_points))


### 6. Blood Lipid Score (non-HDL cholesterol) ###
# Function to calculate *blood lipid score* 
calculate_blood_lipid_points <- function(Blood_lipid, Drug_treated) {
  if (Blood_lipid < 130) {
    Blood_lipid_points <- 100
  } else if (Blood_lipid >= 130 && Blood_lipid <= 159) {
    Blood_lipid_points <- 60
  } else if (Blood_lipid >= 160 && Blood_lipid <= 189) {
    Blood_lipid_points <- 40
  } else if (Blood_lipid >= 190 && Blood_lipid <= 219) {
    Blood_lipid_points <- 20
  } else {
    Blood_lipid_points <- 0
  }
  
  # Subtract 20 points if on drug-treated level (unless score is 0)
  if (Drug_treated == "yes" && Blood_lipid_points > 0) {
    Blood_lipid_points <- Blood_lipid_points - 20
  }
  
  return(Blood_lipid_points)
}

# Apply (store at col 27)
mypheno <- mypheno %>%
  mutate(
    non_HDL_cholesterol = (Cholesterol - HDL.cholesterol)*38.67,
    Blood_lipid_points = mapply(
      calculate_blood_lipid_points,
      Blood_lipid = non_HDL_cholesterol,
      Drug_treated = cholesterol_lowering_medication
    )
  )


### 7. Glucose Score ###
# Function to calculate *Glucose score* 
# should also consider the diabetes status
mypheno <- mypheno %>%
  mutate(HbA1c_percentage = (Glycated.haemoglobin..HbA1c. / 10.93) + 2.15)

# Function to calculate Glucose Score based on HbA1c percentage & diabetes status
calculate_hba1c_score <- function(hba1c_values, diabetes_status) {
  
  hba1c_scores <- numeric(length(hba1c_values))
  
  for (i in seq_along(hba1c_values)) {
    if (diabetes_status[i] == 0) {  # Non-diabetic
      if (hba1c_values[i] < 5.7) {
        hba1c_scores[i] <- 100
      } else if (hba1c_values[i] >= 5.7 && hba1c_values[i] <= 6.4) {
        hba1c_scores[i] <- 60
      } else {
        hba1c_scores[i] <- 0  
      }
    } else if (diabetes_status[i] == 1) {  # Diabetic
      if (hba1c_values[i] < 7.0) {
        hba1c_scores[i] <- 40
      } else if (hba1c_values[i] >= 7.0 && hba1c_values[i] <= 7.9) {
        hba1c_scores[i] <- 30
      } else if (hba1c_values[i] >= 8.0 && hba1c_values[i] <= 8.9) {
        hba1c_scores[i] <- 20
      } else if (hba1c_values[i] >= 9.0 && hba1c_values[i] <= 9.9) {
        hba1c_scores[i] <- 10
      } else {
        hba1c_scores[i] <- 0
      }
    } else {
      hba1c_scores[i] <- 0
    }
  }
  
  return(hba1c_scores)
}

# Apply the scoring function AFTER converting HbA1c to percentage
mypheno <- mypheno %>%
  mutate(HbA1c_score = calculate_hba1c_score(HbA1c_percentage, diabetes))


### 8. Blood Pressure Score ###

# Function to calculate *BP* 
# should consider systolic and diastolic BPs, as well as the treated level
calculate_bp_score <- function(systolic_bp, diastolic_bp, treated_level) {
  
  bp_points <- 0
  
  if (systolic_bp < 120 && diastolic_bp < 80) {
    bp_points <- 100
  } else if (systolic_bp >= 120 && systolic_bp <= 129 && diastolic_bp < 80) {
    bp_points <- 75
  } else if ((systolic_bp >= 130 && systolic_bp <= 139) || (diastolic_bp >= 80 && diastolic_bp <= 89)) {
    bp_points <- 50
  } else if ((systolic_bp >= 140 && systolic_bp <= 159) || (diastolic_bp >= 90 && diastolic_bp <= 99)) {
    bp_points <- 25
  } else if (systolic_bp >= 160 || diastolic_bp >= 100) {
    bp_points <- 0
  } else {
    bp_points <- 0  
  }
  
  # Subtract 20 points if treated level is true (unless score is 0)
  if (treated_level == 1 && bp_points > 0) {
    bp_points <- bp_points - 20
  }
  
  return(bp_points)
}

# Apply (store at col 29)
mypheno <- mypheno %>%
  mutate(
    BP_Score = mapply(
      calculate_bp_score,
      systolic_bp = SBP,
      diastolic_bp = DBP,
      treated_level = BP_Meds.x
    )
  )

#### Final LE8 score data ####
LE8score <- mypheno[, c(1,2,4,21,25,27,28,29,31,33,34)]
LE8score_df <- as.data.frame(LE8score)

colnames(LE8score_df)

write.csv(LE8score_df, "/medpop/esp2/yxliu/LE8score_calculation/updated_LE8score.txt", row.names = FALSE)
