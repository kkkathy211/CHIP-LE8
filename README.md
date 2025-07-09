# CHIP-LE8

**LE8 Calculation:**

/medpop/esp2/yxliu/LE8score_calculation/LE8scores.R

*Output file (txt):*

/medpop/esp2/yxliu/updated_LE8score.txt

**Cross-sectional Analysis:**

/medpop/esp2/yxliu/updated_cross_sectional.R

*Output file (csv):*

/medpop/esp2/yxliu/updated_5CHIP_LE8_association_results_VAF2.csv

*Output plots:*

LE8 composite

/medpop/esp2/yxliu/Plots/LE8_composite_cat_VAF2.png

/medpop/esp2/yxliu/Plots/LE8_composite_cat_VAF10.png

/medpop/esp2/yxliu/Plots/LE8_composite_scaled_VAF2.png

/medpop/esp2/yxliu/Plots/LE8_composite_scaled_VAF10.png

/medpop/esp2/yxliu/Plots/LE8_composite_std_VAF2.png

/medpop/esp2/yxliu/Plots/LE8_composite_std_VAF10.png

LE8 individual

/medpop/esp2/yxliu/Plots/forestplot_ASXL1_std_VAF2.png

/medpop/esp2/yxliu/Plots/forestplot_CHIP_std_VAF2.png

/medpop/esp2/yxliu/Plots/forestplot_DNMT3A_std_VAF2.png

/medpop/esp2/yxliu/Plots/forestplot_nonDNMT3A_std_VAF2.png

/medpop/esp2/yxliu/Plots/forestplot_TET2_std_VAF2.png

**Stratified Analysis:**

**1. Incidence**

  1.1. LE8 composite (“high”, “intermediate”, “low”)
  
/medpop/esp2/yxliu/updated_stratified.R

*Output file (txt):*

/medpop/esp2/yxliu/updated_stratified_cox_results_VAF2.txt

*Output plots:*

/medpop/esp2/yxliu/Plots/incidence_CAD_composite.png

/medpop/esp2/yxliu/Plots/incidence_CAD_death.png

/medpop/esp2/yxliu/Plots/incidence_CAD_intermediate.png

/medpop/esp2/yxliu/Plots/incidence_composite_CVD.png

  1.2.  LE8 each individual
  
/medpop/esp2/yxliu/updated_stratified.R

*Output file (txt):*

/medpop/esp2/yxliu/stratified_cox_results_by_component.txt

*Output plots:*

/medpop/esp2/yxliu/Plots_by_component

**2.  HR**

/medpop/esp2/yxliu/Stratified_HR.R

  2.1 Single-reference “No-CHIP & High” forest plots
  
*Output file (csv):*

/medpop/esp2/yxliu/Results_HR_stratified.txt

*Output plots:*

/medpop/esp2/yxliu/Plots/HR_CAD_Composite.png

/medpop/esp2/yxliu/Plots/HR_CAD_Death_HARD.png

/medpop/esp2/yxliu/Plots/HR_Composite_CVD.png

/medpop/esp2/yxliu/Plots/HR_intermediate_CAD.png

  2.2 Unique ref = No-CHIP_High + within-subgroup comparisons
  
*Output file (csv):*

/medpop/esp2/yxliu/Results_HR_stratified_v2.txt

*Output plots:*

/medpop/esp2/yxliu/Plots/V2_HR_CAD_Composite.png
/medpop/esp2/yxliu/Plots/V2_HR_CAD_Death_HARD.png
/medpop/esp2/yxliu/Plots/V2_HR_Composite_CVD.png
/medpop/esp2/yxliu/Plots/V2_HR_intermediate_CAD.png
