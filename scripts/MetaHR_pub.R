source("scripts/MetaFunc.R")
source("scripts/data_prep/DiagData.R")

# Meta HR  ------------------------------------------------------------------

CategUniHR <- pdcox %>% filHR(analysis = "CategUni", filterN=2)
CategMultiHR <- pdcox %>% filHR(analysis = "CategMulti", filterN=2)

## Blood univariable, Categorical -------------------------------------------

BloodUniHR <- CategUniHR %>% 
  filter(FluidType=="blood") %>% 
  filter(!(Study=="Verde 2019" & TargetBiomarker=="NfL")) %>% #duplicates
  filter(!(str_starts(Study, "Zhu 2024") & TargetBiomarker=="CysC")) %>% #duplicates
  filter(!(str_starts(Study, "Zhu 2024") & TargetBiomarker=="Cre/CysC")) %>% #duplicates
  filter(!(Study=="Mondesert 2025"  & TargetBiomarker=="NfL" & AssayMethod=="Ella, serum")) %>% #duplicates
  filter(!(Study=="Mondesert 2025"  & TargetBiomarker=="NfL" & AssayMethod=="Lumipulse, serum")) %>% #duplicates
  filter(!(Study=="Mondesert 2025"  & TargetBiomarker=="NfL" & AssayMethod=="Elecsys, serum")) %>% #duplicates
  filter(!(Study=="Mondesert 2025"  & TargetBiomarker=="NfL" & AssayMethod=="Simoa, serum")) #duplicates Thouvenot 2020 (>Mondesert 2025), Montpelier Cohort

#CairoPDF(file= "output/HR/HR_blood_univariable.pdf", width=16, height=8)
BloodUniHR %>% metaHR(pattern = 1, #subname = "blood"
                      xlab="Univariable HR (95%CI), Blood Biomarkers, dichotomized comparison") 
#dev.off()

## Summary Table
bloodCategUniSum <- BloodUniHR %>% metaHRsum(fluid_type = "blood")
bloodCategUniSum

## NOTE) Montpelier Cohort: Mondesert 2025, Thouvenot 2020
## Region check
BloodUniHR %>% arrange(TargetBiomarker,Study) %>% select(FluidType, TargetBiomarker, Study, Country, Region)


## Blood multivariable, Categorical -------------------------------------------


BloodMultiHR <- CategMultiHR %>% 
  filter(FluidType=="blood") %>% 
  filter(!(str_starts(Study, "Zhu 2024") & TargetBiomarker=="Cre/CysC")) #duplicates

#CairoPDF(file= "output/HR/HR_blood_multivariable.pdf", width=16, height=14)
BloodMultiHR %>% metaHR(pattern = 1, xlab="Multivariable HR (95%CI), Blood Biomarkers, dichotomized comparison")
#dev.off()


## Summary Table 
bloodCategMultiSum <- BloodMultiHR %>% metaHRsum(fluid_type = "blood")

## Region check
BloodMultiHR %>% arrange(TargetBiomarker,Study) %>% select(FluidType, TargetBiomarker, Study, Country, Region)


## CSF univariable, Categorical -------------------------------------------

CSFUniHR <- CategUniHR %>% filter(FluidType=="CSF") 

#CairoPDF(file= "output/HR/HR_CSF_univariable.pdf", width=14, height=8)
CSFUniHR %>% metaHR(pattern = 1, subname = "CSF",
                    xlab="Univariable HR (95%CI), CSF Biomarkers, dichotomized comparison")
#dev.off()

## Region Check
CSFUniHR %>% arrange(TargetBiomarker,Study) %>% select(FluidType, TargetBiomarker, Study, Country, Region)

## Summary Table
csfCategUniSum <- CategUniHR %>% metaHRsum(fluid_type = "CSF")


## CSF multivariable, Categorical -------------------------------------------

CSFMultiHR <- CategMultiHR %>% filter(FluidType=="CSF")

#CairoPDF(file= "output/HR/HR_CSF_multivariable0623.pdf", width=14, height=10)
CSFMultiHR %>% metaHR(pattern = 1, xlab="Multivariable HR (95%CI), CSF Biomarkers, dichotomized comparison", subname = "CSF")
#dev.off()

# Region Check
CSFMultiHR %>% arrange(TargetBiomarker,Study) %>% select(FluidType, TargetBiomarker, Study, Country, Region)

## Summary Table
csfCategMultiSum <- CSFMultiHR %>% metaHRsum(fluid_type = "CSF")


## 📊 Supple Table 2 -----------------------------------------------------------------

### Covariates, Blood ------------------------------------------------------------

sT2a <- BloodMultiHR %>%
  select(Study, TargetBiomarker,  covar_aHR_death) %>%
  group_by(TargetBiomarker) %>% arrange(TargetBiomarker, Study)%>% ungroup() %>%
  create_grouped_flextable(group_col = "TargetBiomarker") %>%
  flextable::fontsize(size = 10, part = "all") %>%
  flextable::font(fontname = "Times New Roman", part = "all") %>%
  flextable::width(j = 1, width = 0.65) %>%
  flextable::width(j = 2, width = 1.4) %>%
  flextable::width(j = 3, width = 5) # %>% ExFT(landscape = FALSE)  # SuppleTable .docx

sT2a

### Covariates, CSF, Cox HR ----------------------------------------------

sT2b <- CSFMultiHR %>% 
  select(Study, TargetBiomarker, covar_aHR_death) %>% 
  group_by(TargetBiomarker) %>% arrange(TargetBiomarker, Study) %>% ungroup() %>% 
  create_grouped_flextable(group_col = "TargetBiomarker") %>% 
  flextable::fontsize(size = 10, part = "all") %>%   
  flextable::font(fontname = "Times New Roman", part = "all") %>% 
  width(j = 1, width = 0.65) %>%
  width(j = 2, width = 1.4) %>%
  width(j = 3, width = 5) #%>% ExFT(landscape = FALSE)

sT2b

# Make Summary Table ------------------------------------------------------

## Categ HR sum (n >=2) ----------------------------------------------------

CategUniSum <- bind_rows(bloodCategUniSum, csfCategUniSum)
CategMultiSum <- bind_rows(bloodCategMultiSum,csfCategMultiSum)

CategHRsum <- full_join(
  CategUniSum, CategMultiSum, 
  by = c("fluid", "biomarker", "combined"), 
  suffix = c(", univariable", ", multivariable")) %>% 
  select(-combined)  %>% 
  arrange(factor(fluid, levels = c("blood", "CSF"))) 

## 📊 Table 3 -----------------------------------------------------------------

CategHRsum %>% 
  select(-"Egger's test, univariable", -"Egger's test, multivariable") %>% 
  gtHR() 

Table3 <- CategHRsum %>% 
  select(-"Egger's test, univariable", -"Egger's test, multivariable") %>% 
  create_grouped_flextable(group_col = "fluid") %>% 
  align(align = "left", part = "footer") %>% 
  flextable::fontsize(size = 10, part = "all") %>%   # make table more compact
  flextable::font(fontname = "Times New Roman", part = "all") %>% 
  width(j = 1:3, width = 0.6) %>% 
  width(j = 5:6, width = 0.6)

Table3

## Export
# Table3 %>% ExFT(file = "output/SummaryHR.docx")


## Count the number of included studies -------------------------------------

### 📊 Uni & Multi ----------------------------------------------------------

## total study count (n=27 studies, Univariable and Multivariable Cox)

bind_rows(CSFMultiHR,BloodMultiHR) %>% 
  bind_rows(CSFUniHR) %>% bind_rows(BloodUniHR) %>% 
  ungroup %>% distinct(DOI) %>% nrow()


# n=5924 (in total)

bind_rows(CSFMultiHR,BloodMultiHR) %>% 
  bind_rows(CSFUniHR) %>% bind_rows(BloodUniHR) %>% 
  ungroup %>% 
  group_by(Study) %>% arrange(desc(SampleSize_ALS)) %>% ungroup() %>%
  distinct(Study, .keep_all = TRUE) %>% pull(SampleSize_ALS) %>% sum()


### Univariable Cox (n=13) ---------------------------------------------------

## n=13, total ALS cases included (univariable Cox), excluding Cohort duplication

bind_rows(CSFUniHR,BloodUniHR) %>% ungroup() %>% distinct(DOI, .keep_all = TRUE) %>%nrow()

### note) n=15, univariable Cox, study counts (before region check)
CategUniHR %>% ungroup() %>% distinct(DOI) %>% nrow()

## 2055 ALS cases (univariable Cox)

bind_rows(CSFUniHR,BloodUniHR) %>% ungroup() %>%
  group_by(DOI) %>% arrange(desc(SampleSize_ALS)) %>% ungroup() %>% 
  distinct(Study, .keep_all = TRUE) %>%
  summarise(total_sample_size = sum(SampleSize_ALS, na.rm = TRUE))


### Multivariable Cox -------------------------------------------------------

## n=24, total ALS cases included (Multivariable) 

CategMultiHR %>% ungroup() %>% distinct(DOI) %>% nrow()
bind_rows(CSFMultiHR,BloodMultiHR) %>% ungroup() %>% distinct(DOI) %>%nrow()

## 5508 ALS cases (multivariable Cox)

bind_rows(CSFMultiHR,BloodMultiHR) %>% ungroup() %>%
  group_by(DOI) %>% arrange(desc(SampleSize_ALS)) %>% ungroup() %>% 
  distinct(Study, .keep_all = TRUE) %>%   # Guo 2021 (male and female)
  summarise(total_sample_size = sum(SampleSize_ALS, na.rm = TRUE))

# note) Guo 2021 has 2 cohorts (same DOI)
bind_rows(CSFMultiHR,BloodMultiHR) %>% ungroup() %>% 
  select(Study,DOI,SampleSize_ALS) %>% 
  distinct() %>% group_by(DOI) %>% filter(n_distinct(Study) > 1)
