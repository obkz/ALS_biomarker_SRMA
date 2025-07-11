source("scripts/MetaFunc.R")
source("scripts/data_prep/DiagData.R")

## Meta Corr - DPR -----------------------------------------------------------

### 📊 CSF - DPR, Spearman ------------------------------------------------------

csfCORdpr <- pd %>% filterCOR(fluid_type = "CSF", filterN = 2) %>% 
  filter(!Study=="Costa 2019") %>%  #duplicates Costa 2021
  filter(!(Study=="Gagliardi 2021" & TargetBiomarker=="NfL")) %>% #duplicates, Gagliardi 2024 
  filter(!(Study=="Verde 2020" & TargetBiomarker=="pNfH"))        #duplicates, Gagliardi 2021 (>Verde 2020), Milan, IRCCS 

#CairoPDF(file= "output/Corr/DPR_CSF.pdf", width=14, height=18)
csfCORdpr %>% mcforest(sm="COR", subname="CSF", 
           xlab="Correlation of CSF biomarkers with DPR (Spearman's ρ, 95% CI)") 
#dev.off()

# Region Check
RegionCheckCorr <- csfCORdpr %>% select(TargetBiomarker, Study, Country, Region) %>% arrange(TargetBiomarker, Study)
RegionCheckCorr

# CSF DPR summary
csfDPRsum <- csfCORdpr %>% metaCorSUM(fluid_type = "CSF", cor="Corr_DPR", filterN = 2, sm="COR") 
csfDPRsum


### 📊 Blood - DPR, Spearman ---------------------------------------------------

bloodCORdpr <- pd %>% filterCOR(fluid_type = "blood", filterN = 2) %>% 
  filter(!(Study=="Verde 2019" & TargetBiomarker=="NfL")) #duplicates, Verde 2023-1

#CairoPDF(file= "output/Corr/DPR_blood.pdf", width=14, height=18)
bloodCORdpr %>% mcforest(sm="COR", #subname="blood",
  xlab="Correlation of blood biomarkers with DPR (Spearman's ρ, 95% CI)")
#dev.off()

## Region Check
RegionCheckCorr2 <- bloodCORdpr %>% 
  select(TargetBiomarker, Study, Country, Region) %>% arrange(TargetBiomarker, Study)
RegionCheckCorr2

# blood DPR summary
bloodDPRsum <- bloodCORdpr %>% metaCorSUM(fluid_type = "blood", cor="Corr_DPR", filterN = 2, sm="COR")
bloodDPRsum 



## Meta Corr - FRS ---------------------------------------------------------

### 📊 CSF - ALSFRS, Spearman -------------------------------------------------

csfCorFRS <- pd %>% filterCOR(fluid_type = "CSF", cor="Corr_ALSFRSR", filterN = 2)

#CairoPDF(file= "output/Corr/ALSFRS_CSF.pdf", width=14, height=14)
csfCorFRS %>% mcforest(sm="COR", cor="Corr_ALSFRSR",subname = "CSF",
           xlab="Correlation of CSF biomarkers with ALSFRS (Spearman's ρ, 95% CI)")
#dev.off()

## Region Check
RegionCheckCorr3 <- csfCorFRS %>%
  select(TargetBiomarker, Study, Country, Region) %>% arrange(TargetBiomarker, Study)
RegionCheckCorr3 

# CSF ALSFRSR summary
csfFRSRsum <- csfCorFRS %>% metaCorSUM(fluid_type = "CSF", cor="Corr_ALSFRSR", filterN = 2, sm="COR") 
csfFRSRsum


### 📊 Blood - ALSFRS, Spearman ------------------------------------------------

bloodCorFRS <- pd %>% 
  filterCOR(fluid_type = "blood", cor = "Corr_ALSFRSR",filterN = 2) %>% 
  filter(!str_detect(Study, "Sant Pau Cohort ALS-FTD")) %>% #duplicates
  filter(!(Study=="Sun 2021" & Corr_DPR=="NR"))             #duplicates, Sun 2021 IL-2

#CairoPDF(file= "output/Corr/ALSFRS_blood.pdf", width=14, height=21, family = "MyHelv")
bloodCorFRS %>% mcforest(sm="COR", cor="Corr_ALSFRSR",
           xlab="Correlation of Blood biomarkers with ALSFRS (Spearman's ρ, 95% CI)")
#dev.off()

# Region Check
RegionCheckCorr4 <- bloodCorFRS %>%
  select(TargetBiomarker, Study, Country, Region) %>% arrange(TargetBiomarker, Study)
RegionCheckCorr4

# blood ALSFRSR summary
bloodFRSRsum <- bloodCorFRS %>% metaCorSUM(fluid_type = "blood", cor="Corr_ALSFRSR", filterN = 2, sm="COR") 
bloodFRSRsum



# Merge Summary Tables ---------------------------------------------------------

bindDPR <- bind_rows(csfDPRsum,bloodDPRsum)
bindFRSR <-  bind_rows(csfFRSRsum,bloodFRSRsum)

# Merge the two tibbles
merged_Corr <- full_join(
  bindDPR, bindFRSR, by = c("fluid", "biomarker", "combined"), 
  suffix = c(", DPR", ", ALSFRS-R")) %>% 
  select(-combined)  %>% 
  arrange(factor(fluid, levels = c("blood", "CSF")))

# Create a formatted table with subgroups
merged_Corr %>% gt() %>%
  tab_row_group(label = "CSF Biomarkers", rows = fluid == "CSF") %>%
  tab_row_group(label = "Blood Biomarkers", rows = fluid == "blood") %>%
  cols_label(fluid = "Fluid Type",biomarker = "Biomarker") %>%
  tab_options(table.font.size = px(14),row_group.font.weight = "bold") %>% 
  tab_header(md("**pooled Spearman's rho to DPR and baseline ALSFRS-R**"))

## remove Egger's columns
merged_Corr %>%  select(-`Egger's test, DPR`, -`Egger's test, ALSFRS-R`) %>% 
  gt() %>%
  tab_row_group(label = "CSF Biomarkers", rows = fluid == "CSF") %>%
  tab_row_group(label = "Blood Biomarkers", rows = fluid == "blood") %>%
  cols_label(fluid = "Fluid Type",biomarker = "Biomarker") %>%
  tab_options(table.font.size = px(14),row_group.font.weight = "bold") %>% 
  tab_header(md("**pooled Spearman's rho to DPR and baseline ALSFRS-R**")) 


### 📊 Table 4 : Summary Corr Results ---------------------------------------

## FlexTable
Table4 <- merged_Corr %>%  select(-`Egger's test, DPR`, -`Egger's test, ALSFRS-R`) %>% 
  create_grouped_flextable(group_col = "fluid") %>% 
  flextable::fontsize(size = 10, part = "all") %>%   # make table more compact
  flextable::font(fontname = "Times New Roman", part = "all") %>% 
  flextable::width(j = 2:3, width = 0.7) %>%
  flextable::width(j = 4, width = 1.4) %>%
  flextable::width(j = 5:6, width = 0.7) %>%
  flextable::width(j = 7, width = 1.4) # %>% ExFT()   # output/ as .docx

Table4
# showtext_auto(FALSE)


## SampleSize calc -----------------------------------------------------------

### SampleSize DPR calc ------------------------------------------------------

# n=6341 (for DPR calc)

bind_rows(csfCORdpr,bloodCORdpr) %>% 
  group_by(DOI) %>% arrange(desc(SampleSize_ALS)) %>% ungroup() %>% 
  distinct(DOI, .keep_all = TRUE) %>% pull(SampleSize_ALS) %>% sum()

# Studies n=38
bind_rows(csfCORdpr,bloodCORdpr) %>% 
  group_by(DOI) %>% arrange(desc(SampleSize_ALS)) %>% ungroup() %>% 
  distinct(DOI, .keep_all = TRUE) %>% nrow()

### SampleSize FRS calc ------------------------------------------------------

# n=7095 (for FRS calc)

bind_rows(csfCorFRS,bloodCorFRS) %>% 
  group_by(DOI) %>% arrange(desc(SampleSize_ALS)) %>% ungroup() %>% 
  distinct(DOI, .keep_all = TRUE) %>% pull(SampleSize_ALS) %>% sum()

# Studies n=43
bind_rows(csfCorFRS,bloodCorFRS) %>% 
  group_by(DOI) %>% arrange(desc(SampleSize_ALS)) %>% ungroup() %>% 
  distinct(DOI, .keep_all = TRUE) %>% nrow()

### 📊 in Total (DPR + FRS) ----------------------------------------------------

# n=9945 (FRS and DPR)

bind_rows(csfCorFRS,bloodCorFRS) %>% bind_rows(csfCORdpr) %>% bind_rows(bloodCORdpr) %>% 
  group_by(DOI) %>% arrange(desc(SampleSize_ALS)) %>% ungroup() %>% 
  distinct(DOI, .keep_all = TRUE) %>% pull(SampleSize_ALS) %>% sum()

# Studies n=58
bind_rows(csfCorFRS,bloodCorFRS) %>% bind_rows(csfCORdpr) %>% bind_rows(bloodCORdpr) %>% 
  group_by(DOI) %>% arrange(desc(SampleSize_ALS)) %>% ungroup() %>% 
  distinct(DOI, .keep_all = TRUE) %>% nrow()



# Egger's Test ------------------------------------------------------------

## blood NfL, DPR

pd %>% filter(FluidType=="blood", TargetBiomarker=="NfL", CorType=="Spearman's rho") %>% 
  meta::metacor(
    cor = Corr_DPR,       # Extract the specified correlation column
    n = SampleSize_ALS,   # Sample size
    studlab = Study,      # Study labels
    method.tau = "REML",  # Random-effects model
    sm = "COR"
  ) %>% funnel()


## blood NfL, ALSFRS-R

pd %>% filter(FluidType=="blood", TargetBiomarker=="NfL",CorType=="Spearman's rho") %>% 
  meta::metacor(
    cor = Corr_ALSFRSR,   # Extract the specified correlation column
    n = SampleSize_ALS,   # Sample size
    studlab = Study,      # Study labels
    method.tau = "REML",  # Random-effects model
    sm = "COR"
  ) %>% funnel()

