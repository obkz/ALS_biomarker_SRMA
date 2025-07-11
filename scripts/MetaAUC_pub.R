source("scripts/MetaFunc.R")
source("scripts/data_prep/DiagData.R")

# 📊 Study Count for AUC analysis  --------------------------------------------

AUC_count <- combined_df %>% 
  # Keep only the first row per unique DOI
  group_by(DOI) %>% arrange(desc(n1)) %>% ungroup() %>% 
  distinct(DOI, .keep_all = TRUE) %>% mutate(
    # Extract the numeric part after "#" from D_ID (e.g., D#12 → 12)
    D_num = as.numeric(sub("^[A-Z]#", "", D_ID)),
    # Extract the prefix (e.g., "D", "P", "N") from D_ID
    Prefix = sub("#.*", "", D_ID)
  ) %>%
  # Arrange rows: first by Prefix in the order D → P → N, then by numeric value
  arrange(
    factor(Prefix, levels = c("D", "P", "N")),  # D -> P -> N
    D_num  # ascending
  ) %>%
  # Drop the temporary helper columns
  select(-D_num, -Prefix) 

AUC_count%>% nrow()  

## Cases_count -------------------------------------------------------------

# Studies reporting NfL: 72.3 %

AUC_nfl_count <- combined_df %>% 
  filter(TargetBiomarker=="NfL") %>% 
  distinct(DOI, .keep_all = TRUE) %>% nrow()

AUC_nfl_count / nrow(AUC_count)

# All ALS, N = 5566 

AUC_count %>% pull(n1) %>% sum()

# mimics, N = 817 

N_auc_mimic <- combined_df %>% group_by(DOI) %>% 
  filter(Control %in% mdset) %>% 
  arrange(desc(n2)) %>% ungroup() %>% 
  distinct(DOI, .keep_all = TRUE) %>% pull(n2) %>% sum()
N_auc_mimic

# HC, N = 1776 (max)

N_auc_HC <- combined_df %>% group_by(DOI) %>% 
  filter(Control %in% hcset) %>% 
  arrange(desc(n2)) %>% ungroup() %>% 
  distinct(DOI, .keep_all = TRUE) %>% pull(n2) %>% sum()
N_auc_HC

# DC, N = 871 (max)

N_auc_DC <- combined_df %>% 
  group_by(DOI) %>% 
  filter(Control %in% dcset) %>% 
  arrange(desc(n2)) %>% ungroup() %>% 
  distinct(DOI, .keep_all = TRUE) %>% pull(n2) %>% sum()
N_auc_DC

# Control, Total

# Cons: 3464 
N_auc_mimic + N_auc_HC + N_auc_DC

# Total: 9030
sum(AUC_count$n1) + N_auc_mimic + N_auc_HC + N_auc_DC


# CSF ---------------------------------------------------------------------

# Run the combined data through metaAUC()

## 📊 CSF biomarkers vs ALS mimics -----------------------------------

csfAUCmd <- combined_df %>% filter(FluidType=="CSF", Control %in% mdset) %>% 
  group_by(TargetBiomarker) %>%filter(n() >= 2) %>% ungroup() %>%
  mutate(
    TargetBiomarker = factor(TargetBiomarker,
                             levels = c(bmsort, setdiff(unique(TargetBiomarker), bmsort))
    )) %>% arrange(TargetBiomarker, Study) 

#CairoPDF(file= "output/AUC/AUC_CSF_mimics.pdf", width=14, height=16)
csfAUCmd %>% metaAUC(xlab = "AUC (95%CI), CSF biomarkers, vs ALS mimics", subname = "CSF")
#dev.off()

## RegionCheck
RegionCheckAUC1 <- csfAUCmd %>% select(D_ID,Study,TargetBiomarker,Country,Region)
RegionCheckAUC1


## 📊 CSF biomarkers vs neurologically HC -----------------------------------

csfAUChc <- combined_df %>% filter(FluidType=="CSF", Control %in% hcset) %>% 
  group_by(TargetBiomarker) %>%filter(n() >= 2) %>% ungroup() %>%
  mutate(TargetBiomarker = factor(TargetBiomarker,
                                  levels = c(bmsort, setdiff(unique(TargetBiomarker), bmsort))
  )) %>% arrange(TargetBiomarker, Study)

#CairoPDF(file= "output/AUC/AUC_CSF_NHC.pdf", width=14, height=16)
csfAUChc %>% metaAUC(xlab = "AUC (95%CI), CSF biomarkers, vs NHC",subname="CSF")
#dev.off()

## Region Check
RegionCheckAUC2 <- csfAUChc %>% select(D_ID,Study,TargetBiomarker,Country,Region)
RegionCheckAUC2 


## 📊 CSF biomarkers vs neurological DC -----------------------------------


#CairoPDF(file= "output/AUC/AUC_CSF_NDC.pdf", width=14, height=7)

csfAUCdc <- combined_df %>% filter(FluidType=="CSF", Control %in% dcset) %>% 
  group_by(TargetBiomarker) %>%filter(n() >= 2) %>% ungroup() %>%
  mutate(TargetBiomarker = factor(TargetBiomarker,
                                  levels = c(bmsort, setdiff(unique(TargetBiomarker), bmsort))
  )) %>% arrange(TargetBiomarker, Study) 

csfAUCdc %>% metaAUC(xlab = "AUC (95%CI), CSF biomarkers, vs DC", subname="CSF")

#dev.off()

## Region Check
## Note) Shi 2022 is not included in univariate RE model, 
##  - but in bivariate Reitsma model (due to reporting only Sen/Spe w/o AUC)
RegionCheckAUC3 <- csfAUCdc %>% select(D_ID,Study,TargetBiomarker,Country,Region)
RegionCheckAUC3


## 📊 CSF biomarkers vs DC including ALS mimics -----------------------------------

csfAUCmddc <- combined_df %>% filter(FluidType=="CSF", Control %in% mddc) %>% 
  group_by(TargetBiomarker) %>%filter(n() >= 2) %>% ungroup() %>%
  mutate(TargetBiomarker = factor(TargetBiomarker,
                                  levels = c(bmsort, setdiff(unique(TargetBiomarker), bmsort))
  )) %>%
  arrange(TargetBiomarker, Study) %>%  
  filter(!(Study=="Gagliardi 2021" & TargetBiomarker=="NfL")) %>% #duplicates 2021 and 2024
  filter(!(Study=="Verde 2020" & TargetBiomarker=="pNfH")) %>% #duplicates, Gagliardi 2021 (>Verde2020), pNfH
  filter(!(Study=="Dreger 2021" & TargetBiomarker=="NfL" & Control=="ALSmimics")) %>% #duplicates DC and MD
  filter(!(Study=="Agnello 2021" & TargetBiomarker=="p-tau/t-tau" & Control=="DC"))   #duplicates DC and MD
  
#CairoPDF(file= "output/AUC/AUC_CSF_DC_and_mimics.pdf", width=14, height=17)
csfAUCmddc %>% metaAUC(xlab = "AUC (95%CI), CSF biomarkers, vs DC including ALS mimics", subname="CSF")
#dev.off()

## Region Check
RegionCheckAUC4 <- csfAUCmddc %>% select(D_ID,Study,TargetBiomarker,Country,Region)
RegionCheckAUC4


# blood -------------------------------------------------------------------


## 📊 blood, vs ALS mimics (ONLY NfL)-----------------------------------------

bloodAUCmd <- combined_df %>% filter(FluidType=="blood", Control %in% mdset) %>% 
  group_by(TargetBiomarker) %>%filter(n() >= 2) %>% ungroup() %>%
  mutate(TargetBiomarker = factor(TargetBiomarker,
                                  levels = c(bmsort, setdiff(unique(TargetBiomarker), bmsort))
  )) %>% arrange(TargetBiomarker, Study)

#CairoPDF(file= "output/AUC/AUC_blood_mimics.pdf", width=14, height=7)
bloodAUCmd %>% metaAUC(xlab = "AUC (95%CI), blood biomarkers, vs ALS mimics")
#dev.off()

# Region check
# Note) Davies 2023 is not included in univariate RE model, 
#  - but in Reitsma model (due to reporting only Sen/Spe without AUC)
RegionCheckAUC5 <- bloodAUCmd %>% select(D_ID,Study,TargetBiomarker,Country,Region)
RegionCheckAUC5 

## 📊 blood, vs NHC ----------------------------------------------------------

bloodAUChc <- combined_df %>% filter(FluidType=="blood", Control %in% hcset) %>% 
  group_by(TargetBiomarker) %>%filter(n() >= 2) %>% ungroup() %>%
  mutate(TargetBiomarker = factor(TargetBiomarker,
                                  levels = c(bmsort, setdiff(unique(TargetBiomarker), bmsort))
  )) %>% arrange(TargetBiomarker, Study)  

#CairoPDF(file= "output/AUC/AUC_blood_NHC.pdf", width=14, height=9)
bloodAUChc %>% metaAUC(xlab = "AUC (95%CI), blood biomarkers, vs NHC")
#dev.off()

## Region duplicates Check
## Note1) Falzone 2022 and Verde 2023–2 are both in Milan, but OTHER hospitals/groups
## Note2) Sugimoto 2020 is excluded in univariate RE model (since they reported AUC 1.00 without proper CI)

RegionCheckAUC6 <- bloodAUChc %>% select(D_ID,Study,TargetBiomarker,Country,Region)
RegionCheckAUC6 


## 📊 blood, vs DC -----------------------------------------------------------

bloodAUCdc <- combined_df %>% filter(FluidType=="blood", Control %in% dcset) %>% 
  group_by(TargetBiomarker) %>%filter(n() >= 2) %>% ungroup() %>%
  mutate(TargetBiomarker = factor(TargetBiomarker,
                                  levels = c(bmsort, setdiff(unique(TargetBiomarker), bmsort))
  )) %>% arrange(TargetBiomarker, Study) 


#CairoPDF(file= "output/AUC/AUC_blood_DC.pdf", width=14, height=7)
bloodAUCdc %>% metaAUC(xlab = "AUC (95%CI), blood biomarkers, vs DC")
#dev.off()

## Region Check
## Note) Shi 2022 is excluded in univariate RE model (see above)
RegionCheckAUC7 <- bloodAUCdc %>% select(D_ID,Study,TargetBiomarker,Country,Region)
RegionCheckAUC7 

## 📊 blood, vs DC including MD ---------------------------------------------

bloodAUCmddc <- combined_df %>% filter(FluidType=="blood", Control %in% mddc) %>% 
  mutate(TargetBiomarker = factor(TargetBiomarker,
                                  levels = c(bmsort, setdiff(unique(TargetBiomarker), bmsort))
  )) %>%
  arrange(TargetBiomarker, Study) %>%  
  filter(!(Study=="Falzone 2022" & TargetBiomarker=="NfL" & Control=="ALSmimics")) %>% #duplicates 2021 DC and MD
  filter(!(Study=="Verde 2023–1" & TargetBiomarker=="NfL")) %>% #duplicates 2021 DC and MD
  filter(!(Study=="Mondesert 2025" & TargetBiomarker=="NfL")) %>%  #duplicates 2023 Brousse, Montpelier Cohort
  filter(!(str_starts(Study, "Halbgebauer 2022") & TargetBiomarker=="NfL")) %>%  #duplicates 2019 Verde, Ulm Cohort
  filter(!Study=="Shi 2022") %>% # Shi2022 not reporting AUC
  group_by(TargetBiomarker) %>%filter(n() >= 2) %>% ungroup() 

#CairoPDF(file= "output/AUC/AUC_blood_DC_mimics.pdf", width=14, height=10)
bloodAUCmddc %>% metaAUC(xlab = "AUC (95%CI), blood biomarkers, vs DC including ALS mimics")
#dev.off()

## Region Check
## Note) Davies 2023 is excluded in univariate RE model (see above)
## Note) Halbgebauer 2022–2 may be including Verde 2019, Ulm Cohort (N, Verde 2019 >> Halbgebauer 2022-2)
## Note) Brosse 2023 may be including Mondesert 2025, Montpelier Cohort (N, Brousse 2023 >> Mondesert 2025)

RegionCheckAUC8 <- bloodAUCmddc %>% select(D_ID,Study,TargetBiomarker,Country,Region)
RegionCheckAUC8

