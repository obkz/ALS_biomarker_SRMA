source("scripts/MetaFunc.R")
source("scripts/data_prep/DiagData.R")

# Meta SMD -------------------------------------------------------------------

## CSF------------------------------------------------------------------------

### 📊 CSF, vs HC  --------------------------------------------------------------

# font_add( family = "MyHelv", regular = "/System/Library/Fonts/Helvetica.ttc")
# showtext_auto()

smdCSFhc <- sm %>% 
  filter(FluidType=="CSF" & con_sort=="HC") %>% 
  filter(!Study=="Kläppe 2022") %>%  #duplicates: Klappe2024
  filter(!TargetBiomarker=="CXCL12") %>% #from only one group
  filter(!(Study=="Lucchi 2024" & TargetBiomarker=="NfL")) %>%  # duplicates, Martinelli 2024 (> Lucchi2024), Modena University Hospital
  filter(!(Study=="Oeckl 2019 (Ulm gALS)" & TargetBiomarker=="NfL")) %>% #duplicates, Ulm, Halbgebauer 2022–2 (>Oeckl)
  filter(!(Study=="Lucchi 2024" & TargetBiomarker=="NfH")) %>%  # duplicates, Martinelli 2024 (> Lucchi2024), Modena University Hospital
  filter(!(Study=="Lucchi 2024" & TargetBiomarker=="YKL40")) %>%  # duplicates, Martinelli 2024 (> Lucchi2024), Modena University Hospital
  filter(!(Study=="Lucchi 2024" & TargetBiomarker=="SerpinA1")) %>%  # duplicates, Martinelli 2024 (> Lucchi2024), Modena University Hospital
  filter(!(Study=="Cousins 2022" & TargetBiomarker=="Aβ42")) %>%  # pens USA, same cohort (Olsson2019 > Cousins2022)
  group_by(TargetBiomarker) %>% filter(n() >= 2) %>%ungroup()

#CairoPDF(file = "output/SMD/SMD_CSF_HC.pdf", width=13, height=24, family = "MyHelv")
smdCSFhc %>% SMDplot(pattern = 2, Fluid="CSF", Con = "HC", filterN=2,
                     xlab="CSF biomarkers SMD (ALS vs HC), 95%CI", subname = "CSF")
#dev.off()


# Region check
# Note) Hafsteinsdóttir 2024, Rosen 2024 = Same region, but other time frame (Gothenburg, Sweden)
# Note) Abu-Rumeileh 2020 and Abu-Rumeileh2025 are different cohorts 

RegionCheckSMD1 <- smdCSFhc %>%
  group_by(TargetBiomarker) %>%arrange(TargetBiomarker, Study) %>% 
  select(Study, TargetBiomarker, Country, Region, SampleSize_ALS, SampleSize_Con)
RegionCheckSMD1

## Summary Table
csfSMDsumHC <- smdCSFhc %>%  
  metaSMDsum(fluid_type = "CSF", control = "HC", filterN = 2)

csfSMDsumHC


### 📊 CSF, vs ALS mimics  -----------------------------------------

smdCSFmd <- sm %>% 
  filter(FluidType=="CSF" & con_sort=="ALS mimics") %>% 
  filter(!Study=="Kläppe 2022") %>%  #duplicates: Klappe2024
  filter(!(Study == "Thompson 2019" & Control == "PLS")) %>% #duplicates
  filter(!(Study == "Sabbatini 2021" & Control == "PN")) %>%  #duplicates
  group_by(TargetBiomarker) %>% filter(n() >= 2) %>%ungroup()

#CairoPDF(file = "output/SMD/SMD_CSF_mimics.pdf", width=13, height=13, family = "MyHelv")
smdCSFmd %>% SMDplot(pattern = 2, Fluid="CSF", Con = "ALS mimics", filterN=2,
                     xlab="CSF biomarkers SMD (ALS vs mimics), 95%CI", subname = "CSF")
#dev.off()

#Region duplicates Check 
RegionCheckSMD2 <- smdCSFmd %>% group_by(TargetBiomarker) %>% arrange(TargetBiomarker, Study) %>%
    select(Study, TargetBiomarker, Country, Region, SampleSize_ALS, SampleSize_Con)
RegionCheckSMD2

## Summary Table
csfSMDsumMD <- smdCSFmd %>% 
  metaSMDsum(fluid_type = "CSF", control="ALS mimics",filterN = 2)

csfSMDsumMD


### 📊 CSF, vs DC & mimics  -----------------------------------------

smdCSFmddc <- sm %>% 
  filter(Control %in% mddc) %>% 
  filter(!Control=="AD") %>% 
  filter(!TargetBiomarker=="KYNA") %>% #only 1 study
  filter(!TargetBiomarker=="TDP43 misfolding") %>% #only 1 study
  filter(!TargetBiomarker=="UCHL1") %>% #only 1 study
  filter(!TargetBiomarker=="p-tau") %>% #only 1 retracted study (Gong 2022)
  filter(!(Study=="Kläppe 2022" & TargetBiomarker=="NfL")) %>%  #duplicates: Klappe2024
  filter(!(Study=="Vacchiano 2023" & Control=="AD" & TargetBiomarker=="NfL")) %>%  #duplicates: Vacciano 2021 mimics
  filter(!(Study=="Dreger 2021" & Control=="bvFTD" & TargetBiomarker=="NfL")) %>%  #duplicates: Dreger 2021 mimics
  filter(!(Study == "Thompson 2019" & Control == "PLS")) %>% #duplicates Thompson 2019 mimics
  filter(!(Study == "Sabbatini 2021" & Control == "Dementia" & TargetBiomarker=="NfL")) %>% #duplicates Sabbatini 2021 SMA
  filter(!(Study == "Sun 2020" & Control == "IMPN" & TargetBiomarker=="NfL")) %>% #duplicates Sun 2020 DC
  filter(!(Study == "Sun 2020" & Control == "NIMPN" & TargetBiomarker=="NfL")) %>% #duplicates Sun 2020 DC
  filter(!(Study == "Agnello 2021" & Control == "DC")) %>% #duplicates Agnello 2021 Tau proteins, mimics
  filter(!(Study == "Costa 2019" & TargetBiomarker == "pNfH")) %>% #duplicates Costa 2021 (>Costa 2019)
  filter(!(Study=="Halbgebauer 2022–1" & Control=="bvFTD")) %>% #t-tau
  filter(!(Study == "Shi 2022" & Control == "DC")) %>%   # med IQR => SMD, extreme value, unreliable
  group_by(TargetBiomarker) %>% filter(n() >= 2) %>%ungroup()


#Region duplication Check 
RegionCheckSMD3 <- smdCSFmddc %>%
  filter(FluidType=="CSF") %>%
  group_by(TargetBiomarker) %>% filter(n() >= 2) %>%
  arrange(TargetBiomarker, Study) %>%
  select(TargetBiomarker, Study,Country, Region, SampleSize_ALS, SampleSize_Con)

RegionCheckSMD3

# Export

#CairoPDF(file= "output/SMD/SMD_CSF_DC_and_mimics.pdf", width=14, height=18)
smdCSFmddc %>% mutate(TargetBiomarker = str_replace(TargetBiomarker, "^ptau_ttau_ratio$", "p-tau/t-tau")) %>%
  SMDplot(pattern = 2, Fluid="CSF", filterN=2, subname = "CSF",
          xlab="CSF biomarkers SMD (ALS vs DC, including ALS mimics), 95%CI",)
#dev.off()


csfSMDsumDC <-smdCSFmddc %>%  
  mutate(TargetBiomarker = str_replace(TargetBiomarker, "^ptau_ttau_ratio$", "p-tau/t-tau")) %>%
  metaSMDsum(fluid_type = "CSF", filterN = 2)

csfSMDsumDC


### CSF, vs Asymptomatic career -----------------------------------------

sm %>% SMDplot(pattern = 2, Fluid="CSF", Con = "asymptomatic carrier", filterN=2,
          xlab="CSF biomarkers SMD (ALS vs carrier), 95%CI")


## Blood -----------------------------------------------------------

### 📊 blood, vs HC ---------------------------------------------------

## Filter

smdBloodHC <- sm %>% filter(FluidType=="blood" & con_sort=="HC") %>% 
  # omit irrelevant biomarkers
  filter(
    !TargetBiomarker %in% 
      c("angiogenin","angiogenin in plasma exosome fraction",
        "CXCR3+EOMES+","TDP43","IL-18"), #only 1 study for each
    !Study=="Vacchiano 2021" #duplicates: Vacchiano 2023
  ) %>% 
  filter(!(Study=="Lucchi 2024" & TargetBiomarker=="NfL")) %>%  # duplicates, Martinelli 2024 (> Lucchi2024), Modena University Hospital
  filter(!(Study=="Lucchi 2024" & TargetBiomarker=="NfH")) %>%  # duplicates, Martinelli 2024 (> Lucchi2024), Modena University Hospital
  filter(!(Study=="Lucchi 2024" & TargetBiomarker=="YKL40")) %>%  # duplicates, Martinelli 2024 (> Lucchi2024), Modena University Hospital
  filter(!(Study=="Lucchi 2024" & TargetBiomarker=="SerpinA1")) %>%  # duplicates, Martinelli 2024 (> Lucchi2024), Modena University Hospital
  filter(!(Study=="Halbgebauer 2022–2" & TargetBiomarker=="NfL")) %>%  #duplicates, Witzel 2024 (>Halbgebauer 2022–2), Ulm Univ.
  filter(!(Study=="Mastrangelo 2023" & TargetBiomarker=="p-tau181")) %>%    #duplicates, Vacchiano 2023 (>Mastrangelo 2023), Institute of Neurological Sciences of Bologna
  filter(!(Study=="Mastrangelo 2023" & TargetBiomarker=="NfL")) %>%  #duplicates, Vacchiano 2023 (>Mastrangelo 2023), Institute of Neurological Sciences of Bologna  
  group_by(TargetBiomarker) %>% filter(n() >= 2) %>%ungroup()


#CairoPDF(file= "output/SMD/SMD_blood_HC.pdf", width=18, height=34)
smdBloodHC %>% SMDplot(pattern = 2, Fluid="blood", Con = "HC", filterN=2,
                       xlab="Blood biomarkers SMD (ALS vs NHC), 95%CI")
#dev.off()


## Region Check
RegionCheckSMD3 <- smdBloodHC %>%
  group_by(TargetBiomarker) %>% arrange(TargetBiomarker, Study) %>%
  select(Study, TargetBiomarker, Country, Region, SampleSize_ALS, SampleSize_Con)
RegionCheckSMD3 

## Summary Table
bloodSMDsumHC <- smdBloodHC %>% 
  metaSMDsum(fluid_type = "blood", control = "HC", filterN = 2)

bloodSMDsumHC


### 📊 blood, vs ALS mimics (-> NfL,GFAP only) ---------------------------

smdBloodmd <- sm %>%
  filter(FluidType=="blood" & con_sort=="ALS mimics") %>% 
  filter(!TargetBiomarker %in% c("peripherin")) %>%  #just 1 study reporting Peripherin
  filter(!(TargetBiomarker=="NfL" & Study=="Mondesert 2025" & !(AssayMethod=="Simoa, serum"))) %>% #duplicates, Mondesert 2025
  group_by(TargetBiomarker) %>% filter(n() >= 2) %>%ungroup()

#CairoPDF(file= "output/SMD/SMD_blood_mimics.pdf", width=14, height=7)
smdBloodmd %>% SMDplot(pattern = 2, Fluid="blood", Con = "ALS mimics", filterN=2,
          xlab="Blood biomarkers SMD (ALS vs mimics), 95%CI")
#dev.off()

# Region Check DONE 2025-05-30

RegionCheckSMD4 <- smdBloodmd %>%
    group_by(TargetBiomarker) %>% arrange(TargetBiomarker, Study) %>%
    select(Study, TargetBiomarker, Country, Region, SampleSize_ALS, SampleSize_Con)
RegionCheckSMD4 


## Summary Table

bloodSMDsumMD <- smdBloodmd %>% 
  metaSMDsum(fluid_type = "blood", control ="ALS mimics", filterN = 2)

bloodSMDsumMD


### 📊 blood, vs DC & mimics  -----------------------------------------

smdBloodmddc <- sm %>% 
  filter(FluidType=="blood" & Control %in% mddc) %>% 
  filter(!Control=="AD") %>%  #removing deLuna
  filter(!TargetBiomarker %in% c("C3", "C4", "IgG", "IgA", "IgM","peripherin")) %>%  #each reported by only 1 study
  filter(!(TargetBiomarker=="NfL" & Study=="Mondesert 2025" & !(AssayMethod=="Simoa, serum"))) %>% #duplicates, Mondesert 2025
  filter(!(Study == "Mondesert 2025" & TargetBiomarker == "NfL")) %>% #duplicates, Brousse 2023 (>Mondesert 2025)
  filter(!(TargetBiomarker=="NfL" & Study=="Ashrafzadeh-Kian 2024" & !(AssayMethod=="Simoa, plasma"))) %>% #duplicates, Ashrafzadeh-Kian 2024, other Assays
  filter(!(Study=="Vacchiano 2023" & Control=="AD" & TargetBiomarker=="NfL")) %>%  #duplicates: Vacciano 2021 mimics
  filter(!(str_starts(Study, "Chatterjee 2024") & Control=="bvFTD")) %>% #duplicates, vs PSP
  filter(!(Study == "Brousse 2023" & AssayMethod == "Ella, serum")) %>% #duplicates, Simoa serum
  filter(!(Study == "Falzone 2022" & Control == "DC")) %>% #duplicates, vs mimics
  filter(!(Study == "Simonini 2021" & Control == "DC" & TargetBiomarker=="pNfH")) %>% #duplicates, vs mimics
  group_by(TargetBiomarker) %>% filter(n() >= 2) %>%ungroup()

#CairoPDF(file= "output/SMD/SMD_blood_DC_and_mimics.pdf", width=14, height=15)
smdBloodmddc %>% SMDplot(pattern = 2, Fluid="blood", filterN=2,
          xlab="Blood biomarkers SMD (ALS vs DC, including ALS mimics), 95%CI")
#dev.off()

## RegionCheck
RegionCheckSMD5 <- smdBloodmddc %>%
  group_by(TargetBiomarker) %>% arrange(TargetBiomarker, Study) %>%
  select(Study, TargetBiomarker, Country, Region, SampleSize_ALS, SampleSize_Con)
RegionCheckSMD5

bloodSMDsumDC <- smdBloodmddc %>% metaSMDsum(fluid_type = "blood", filterN = 2)
bloodSMDsumDC


# Summary Table w/o Egger -------------------------------------------------

# HC
SMDsumHC <- bind_rows(bloodSMDsumHC,csfSMDsumHC) %>% select(-`Egger's test`)
# mimics
SMDsumMD <- bind_rows(bloodSMDsumMD,csfSMDsumMD) %>% select(-`Egger's test`)
# DC including mimics
SMDsumDC <- bind_rows(bloodSMDsumDC,csfSMDsumDC) %>% select(-`Egger's test`)

## overall (HC, ALS mimics, DC including mimics)
sumhcmd <- full_join(SMDsumHC,SMDsumMD,
                     by=c("fluid","biomarker","combined"),
                     suffix = c("†", "‡")) 

sumhcmd %>% full_join(SMDsumDC,by=c("fluid","biomarker","combined"),
                      suffix = c("", "⁑")) %>% select(-combined) %>% gtSMD()


# 📊 Supple T1 ----------------------------------------------------------------

sT1 <- sumhcmd %>% full_join(SMDsumDC,by=c("fluid","biomarker","combined"),
                             suffix = c("", "⁑")) %>% select(-combined) 


suppleT1 <- sT1 %>% create_grouped_flextable(group_col = "fluid") %>% 
  flextable::align(align = "left", part = "footer") %>% 
  flextable::fontsize(size = 10, part = "all") %>%   # make table more compact
  flextable::font(fontname = "Times New Roman", part = "all") %>% 
  flextable::width(j = 1:4, width = 0.55) %>% 
  flextable::width(j = 6:8, width = 0.55) %>% 
  flextable::width(j = 10:12, width = 0.55)

suppleT1 

#suppleT1 %>% ExFT("test.docx")


## Summary Table w/ Egger------------------------------------------------------

#HC
SMDsumHC <- bind_rows(bloodSMDsumHC,csfSMDsumHC)
#mimics
SMDsumMD <- bind_rows(bloodSMDsumMD,csfSMDsumMD)
#DC including mimics
SMDsumDC <- bind_rows(bloodSMDsumDC,csfSMDsumDC)

## overall (HC, ALS mimics, DC including mimics)
sumhcmd <- full_join(SMDsumHC,SMDsumMD, by=c("fluid","biomarker","combined"),
                     suffix = c("†", "‡")) 

sumhcmd %>% full_join(SMDsumDC,by=c("fluid","biomarker","combined"),
                      suffix = c("", "⁑")) %>% select(-combined) %>% gtSMD()

sumhcmd %>% full_join(SMDsumDC,by=c("fluid","biomarker","combined"),
                      suffix = c("", "⁑")) %>% select(-combined)



# 📊 SMD study-count ---------------------------------------------------------

doicount <- smdCSFhc %>% 
  bind_rows(smdCSFmd) %>% 
  bind_rows(smdCSFmddc) %>% 
  bind_rows(smdBloodHC) %>% 
  bind_rows(smdBloodmd) %>% 
  bind_rows(smdBloodmddc) %>% 
  group_by(DOI) %>% arrange(desc(SampleSize_ALS)) %>% ungroup() %>% 
  distinct(DOI, .keep_all = TRUE) 

# Study count, n=71 
doicount %>% nrow()
doicount %>% filter(FluidType=="CSF") %>% nrow()
doicount %>% filter(FluidType=="blood") %>% nrow()

# ALS, n=9025
sum(doicount$SampleSize_ALS)

# Mimics, n=712
smdmd_con <- smdCSFmd %>% 
  bind_rows(smdBloodmd) %>% 
  group_by(Study) %>% arrange(desc(SampleSize_Con)) %>% ungroup() %>% 
  distinct(Study, .keep_all = TRUE) %>% pull(SampleSize_Con) %>% sum()
smdmd_con

# mimics & DCs, n=2085
smdmddc_con <- smdCSFmddc %>% 
  bind_rows(smdBloodmddc) %>% 
  group_by(DOI) %>% arrange(desc(SampleSize_Con)) %>% ungroup() %>% 
  distinct(DOI, .keep_all = TRUE) %>% pull(SampleSize_Con) %>% sum()

smdmddc_con

# DC, n=1373  
smdmddc_con-smdmd_con

# HC, n=4694
smdCSFhc %>% 
  bind_rows(smdBloodHC) %>% 
  group_by(DOI) %>% arrange(desc(SampleSize_Con)) %>% ungroup() %>% 
  distinct(DOI, .keep_all = TRUE) %>% pull(SampleSize_Con) %>% sum()


