# Load shared functions and data preparation scripts

source("scripts/MetaFunc.R")
source("scripts/data_prep/DiagData.R")

# AUC/HSROC ---------------------------------------

# Combine all relevant diagnostic datasets into a single data frame
combined_df <- bind_rows(
  ## ========= CSF =========
  # CSF NfL 
  csf_nfl_mimics,  # CSF NfL vs ALS mimics
  csf_nfl_dc,      # CSF NfL vs disease controls
  csf_nfl_hc,      # CSF NfL vs healthy controls
  # CSF CHIT1
  CCM,             # CSF CHIT1 vs ALS mimics
  CCH,             # CSF CHIT1 vs healthy controls
  # CSF YKL40
  CYM,             # CSF YKL40 vs ALS mimics
  CYC,             # CSF YKL40 vs all controls
  # CSF pNfH
  CPM,             # CSF pNfH vs ALS mimics
  CPD,             # CSF pNfH vs disease controls
  CPH,             # CSF pNfH vs healthy controls
  # CSF NfH
  csf_nfh_hc,      # CSF NfH vs healthy controls
  # CSF t-tau
  csf_ttau_mimics, # CSF total‐tau vs ALS mimics
  csf_ttau_hc,
  # CSF p-tau181
  csf_ptau_con,    # CSF phosphorylated‐tau vs controls
  # CSF p-tau/t-tau
  csf_ptr_mimics,  # CSF p/t‐tau ratio vs ALS mimics
  csf_ptr_con,     # CSF p/t‐tau ratio vs controls
  # CSF MCP1
  csf_mcp1_mimics, # CSF MCP-1 vs ALS mimics
  # CSF TDP43
  csf_tdp43_hc,    # CSF TDP-43 vs healthy controls
  
  ## ========= Blood =========
  # Blood NfL
  BNM,             # Blood NFL vs ALS mimics
  BND,             # Blood NFL vs disease controls
  BNH,             # Blood NFL vs healthy controls
  # Blood pNfH
  blood_pnfh_hc,   # Blood pNfH vs healthy controls
  blood_pnfh_dc,   # Blood pNfH vs disease controls
  # Blood NfH
  blood_NfH_con,   # Blood NfH vs all controls
  # Blood GFAP
  blood_GFAP_hc,   # Blood GFAP vs healthy controls
  # Blood UCHL1
  blood_UCHL1_mddc # Blood UCH-L1 vs mixed disease & healthy controls
) %>%
  # Remove duplicate rows based on core study-level identifiers
  distinct(Study, FluidType, TargetBiomarker, Control, n1, n2, .keep_all = TRUE) %>% 
  # Harmonize biomarker naming
  mutate(TargetBiomarker = ifelse(TargetBiomarker == "ptau_ttau_ratio", "p-tau/t-tau", TargetBiomarker))


# Notes:
# - Studies not reporting AUC are excluded from univariate random-effects meta-analysis of AUC
# - These studies are retained for bivariate HSROC analysis (Reitsma et al. model)

combined_df %>% filter(is.na(AUC))


# SMD ------------------------------------------------------------

# Variables required for SMD-based meta-analysis
smdset <- c("Study","Control","Country","Region",
            "FluidType","TargetBiomarker","AssayMethod","MeasuredLevel",
            "Mean_ALS", "SD_ALS", "Mean_Con", "SD_Con",
            "SampleSize_ALS","SampleSize_Con","cases_with_FH","%FH",
            "genetic_cases",
            "ALS_age_mean","ALS_age_sd","ALS_age_median","ALS_age_1Q","ALS_age_3Q",
            "Con_age_mean","Con_age_SD","Con_age_median","Con_age_1Q","Con_age_3Q",
            "Med_ALS","Q1_ALS","Q3_ALS","Med_Con","Q1_Con","Q3_Con",
            "Country","Region",
            "ALS_setting_primary", "ALS_setting_secondary","DOI")

# Minimal dataset used for downstream SMD calculations
mini <- c("Study","Control","con_sort","Country","Region","FluidType","TargetBiomarker",
          "SampleSize_ALS","SampleSize_Con",
          "AssayMethod","MeasuredLevel",
          "Mean_ALS", "SD_ALS", "Mean_Con", "SD_Con",
          "Med_ALS","Q1_ALS","Q3_ALS","Med_Con","Q1_Con","Q3_Con","DOI")

smdata <- read_excel("data/diagnosis_prognosis_pub.xlsx") %>%
  mutate(across(where(is.character), ~na_if(.x, "NR"))) %>% 
  mutate(across(c(Mean_ALS, Mean_Con, SD_ALS, SD_Con,SampleSize_ALS,SampleSize_Con,
                  Med_ALS, Q1_ALS, Q3_ALS, Med_Con, Q1_Con, Q3_Con), as.numeric)) %>% 
  mutate(
    ALS_setting_primary = sub(",.*", "", `ALS setting`),   # Extract information before the comma
    ALS_setting_secondary = sub(".*, ", "", `ALS setting`) # Extract information after the comma
  ) %>% 
  select(all_of(smdset)) %>% 
  mutate(
    origFluid = FluidType,  
    FluidType = if_else(FluidType %in% c("serum", "plasma", "blood"), "blood", FluidType),
    AssayMethod = if_else(origFluid %in% c("serum", "plasma"),
                          paste0(AssayMethod, ", ", origFluid),
                          AssayMethod))  %>% 
  mutate(
    con_sort = case_when(
      Control %in% mdset ~ "ALS mimics",    
      Control %in% hcset ~ "HC", 
      Control %in% dcset~ "DC", 
      TRUE ~ Control # Others. as original value
    )) %>% 
  filter(
    (
      (!is.na(Mean_ALS) & !is.na(SD_ALS)) |
        (!is.na(Q1_ALS) & !is.na(Med_ALS) & !is.na(Q3_ALS) & Q1_ALS < Med_ALS & Med_ALS < Q3_ALS)
    ) &
      (
        (!is.na(Mean_Con) & !is.na(SD_Con)) |
          (!is.na(Q1_Con) & !is.na(Med_Con) & !is.na(Q3_Con) & Q1_Con < Med_Con & Med_Con < Q3_Con)
      )
  ) 

sm <- smdata %>% 
  select(all_of(mini)) %>%
  filter(!is.na(SampleSize_ALS) & SampleSize_ALS > 10 & !is.na(SampleSize_Con)) %>%
  mutate(con_sort = factor(con_sort, 
                           levels = c("ALS mimics", "DC", "HC", 
                                      setdiff(unique(con_sort), c("ALS mimics", "DC", "HC"))),
                           ordered = TRUE),
         TargetBiomarker = str_replace(TargetBiomarker, "^ptau_ttau_ratio$", "p-tau/t-tau"),
         TargetBiomarker = factor(TargetBiomarker,
                                  levels = c(bmsort, 
                                             setdiff(unique(TargetBiomarker), bmsort)),
                                  ordered = TRUE))


# HR and correlation -----------------------------------------------------------

# Variables required for prognostic analyses
sel <- c("Study","Country","Region","FluidType","TargetBiomarker", "AssayMethod",
         "Country","Region","SampleSize_ALS","HR_scaled",
         "n_HR_death","HR_def","HR_scaled","HR_death_cont","HR_death_cont_lCI","HR_death_cont_uCI",
         "HR_death_categ_def","HR_death_categ","HR_death_categ_lCI","HR_death_categ_uCI",
         "covar_aHR_death","n_aHR_death","aHR_death_cont","aHR_death_cont_lCI","aHR_death_cont_uCI",
         "aHR_death_categ","aHR_death_categ_lCI","aHR_death_categ_uCI","aHR_death_categ_def",
         "CorType",
         "Corr_DPR","p_cor_DPR","Corr_ALSFRSR","p_cor_ALSFRSR","Corr_SurvTime","p_cor_SurvT",
         "Corr_slope","p_cor_slope",
         "ALS_setting_primary", "ALS_setting_secondary","DOI")


prog <- read_excel("data/diagnosis_prognosis_pub.xlsx", sheet="Prognostic") %>%
  select(-Volume,-Issue,-Pages) %>% 
  mutate(across(where(is.character), ~na_if(.x, "NR"))) %>% 
  mutate(across(c(SampleSize_ALS,
                  ALS_age_mean,ALS_age_sd,ALS_age_median,ALS_age_1Q,ALS_age_3Q,
                  n_HR_death,HR_death_cont,HR_death_cont_lCI,HR_death_cont_uCI,
                  HR_death_categ,HR_death_categ_lCI,HR_death_categ_uCI,
                  n_aHR_death,aHR_death_cont,aHR_death_cont_lCI,aHR_death_cont_uCI,
                  aHR_death_categ,aHR_death_categ_lCI,aHR_death_categ_uCI,
                  Corr_DPR,p_cor_DPR,Corr_ALSFRSR,p_cor_ALSFRSR,Corr_SurvTime,p_cor_SurvT,
                  Corr_slope,p_cor_slope), as.numeric)) %>% 
  mutate(FluidGroup = case_when(
    FluidType == "CSF" ~ "CSF",
    FluidType == "serum" ~ "Serum",
    FluidType == "plasma" ~ "Plasma",
    TRUE ~ "Others" # others
  )) %>% 
  mutate(
    ALS_setting_primary = sub(",.*", "", `ALS setting`), # Extract information before the comma
    ALS_setting_secondary = sub(".*, ", "", `ALS setting`) # Extract information after the comma
  ) %>% select(all_of(sel))


pd <- prog   %>% 
  # set TargetBiomarker's order
  mutate(
    TargetBiomarker = str_replace(TargetBiomarker, "^ptau_ttau_ratio$", "p-tau/t-tau"),
    TargetBiomarker = factor(TargetBiomarker,
                             levels = c(bmsort, 
                                        setdiff(unique(TargetBiomarker), bmsort)),
                             ordered = TRUE)) %>% 
  # merge Plamsa AND Serum => Blood
  mutate(
    origFluid = FluidType,  
    FluidType = if_else(FluidType %in% c("serum", "plasma", "blood"), "blood", FluidType),
    AssayMethod = if_else(origFluid %in% c("serum", "plasma"),
                          paste0(AssayMethod, ", ", origFluid),
                          AssayMethod),
  ) 

# Harmonize HR scaling notation and fix legacy labels
pdcox <- pd %>% mutate(
  HR_scaled = if_else(is.na(HR_scaled), "n.r.", HR_scaled), 
  HR_scaled = if_else(HR_scaled == "logn", "log", HR_scaled)      # logn -> log
) 


