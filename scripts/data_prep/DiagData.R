
source("scripts/MetaFunc.R")

# Data  -------------------------------------------------------------------

num_cols <- c(
  "SampleSize_ALS","SampleSize_Con",
  "cases_with_FH", "genetic_cases",
  "ALS_age_mean",  "ALS_age_sd",   "ALS_age_median", "ALS_age_1Q",   "ALS_age_3Q",
  "Con_age_mean", "Con_age_SD",   "Con_age_median", "Con_age_1Q",   "Con_age_3Q",
  "Sensitivity",  "Sens_lowerCI", "Sens_upperCI",
  "Specificity",  "Spe_lowerCI",  "Spe_upperCI",
  "AUC",  "AUC_lowerCI",  "AUC_upperCI",
  "Mean_ALS",     "SD_ALS",      "SE_ALS",      "Med_ALS",   "Q1_ALS",   "Q3_ALS", "IQR_ALS",
  "Mean_Con",     "SD_Con",      "SE_Con",      "Med_Con",   "Q1_Con",   "Q3_Con", "IQR_Con"
)

data <- read_excel("data/diagnosis_prognosis_pub.xlsx") %>%
  mutate(across(where(is.character),
                ~ replace(.x, .x %in% c("NR", "–", "(SR w/o MA)"), NA)
  )) %>%
  # convert all specified columns to numeric
  mutate(across(all_of(num_cols), as.numeric))  %>%
  # compute additional variables
  mutate(FluidGroup = case_when(
    FluidType == "CSF" ~ "CSF",
    FluidType == "serum" ~ "Serum",
    FluidType == "plasma" ~ "Plasma",
    TRUE ~ "Others" # others
  )) %>%
  mutate(
    Sensitivity = as.numeric(Sensitivity)/100,
    Specificity = as.numeric(Specificity)/100,
    TP = round(Sensitivity * SampleSize_ALS),
    FN = SampleSize_ALS - TP,
    TN = round(Specificity * SampleSize_Con),
    FP = SampleSize_Con - TN,
    #Standard error of AUC (from CI when available)
    AUC_SE = (AUC_upperCI - AUC_lowerCI) / (2 * 1.96),
    Sensitivity_SE = sqrt((Sensitivity * (1 - Sensitivity)) / (TP + FN)),
    Specificity_SE = sqrt((Specificity * (1 - Specificity)) / (TN + FP)),
    Sensitivity_calc = CalcSENS(TP,FN),
    Specificity_calc = CalcSPEC(TN,FP)
  ) %>%
  mutate(
    n1 = TP + FN,  # positive cases
    n2 = TN + FP,  # negative cases
    n3 = n1 + n2,  # all cases
    calcAUC = ifelse(
      is.na(AUC),  # Check if AUC is NA
      NA,
      AUC          # Use the original AUC if it exists
    ),
    calcAUC_SE = ifelse(
      is.na(AUC_SE) & !is.na(AUC),
      #Hanley & McNeil (1982)
      sqrt((AUC * (1 - AUC) + (n1 - 1) * ((AUC / (2 - AUC)) - AUC^2) + (n2 - 1) * ((2 * AUC^2 / (1 + AUC)) - AUC^2)) / (n1 * n2)),
      AUC_SE
    )
  )%>%
  mutate(
    ALS_setting_primary = sub(",.*", "", `ALS setting`), # Extract information before the comma
    ALS_setting_secondary = sub(".*, ", "", `ALS setting`) # Extract information after the comma
  ) %>%
  # Merge plasma and serum into a single "blood" category
  mutate(
    origFluid = FluidType,
    FluidType = if_else(FluidType %in% c("serum", "plasma", "blood"), "blood", FluidType),
    AssayMethod = if_else(origFluid %in% c("serum", "plasma"),
                          paste0(AssayMethod, ", ", origFluid),
                          AssayMethod)
  ) %>%
  mutate(ControlCateg = case_when(
    Control %in% mdset ~ "mimics",
    Control %in% hcset ~ "NHC",
    Control %in% dcset ~ "DC",
    TRUE ~ NA_character_
  ))

## Variables used for downstream calculations ------------------------------

SET <- c("D_ID","Study","Country","Region","FluidType","TargetBiomarker","Control",
         "ControlCateg", "n1", "n2", "n3",
         "Sensitivity", "Specificity", "Sensitivity_calc","Specificity_calc",
         "Sensitivity_SE", "Specificity_SE","TP", "FN", "FP", "TN", "AUC", 
         "Cut_off_Value","AssayMethod",
         "AUC_lowerCI", "AUC_upperCI", "AUC_SE", "calcAUC", "calcAUC_SE", 
         "Mean_ALS", "SD_ALS", "Mean_Con", "SD_Con",
         "ALS_setting_primary", "ALS_setting_secondary", "DOI")


# CSF analysis ------------------------------------------------------------

## CSF NfL ----------------------------------------------------------------

csf_nfl <- data %>%
  filter(!is.na(Sensitivity) & !is.na(Specificity))%>%
  filter(FluidType == "CSF", TargetBiomarker == "NfL") %>% 
  select(any_of(SET)) %>% 
  mutate(
    Cut_off_Value = ifelse(is.na(Cut_off_Value), "not specified", Cut_off_Value)
  )

### CSF NfL vs mimics ---------------------------------------------

csf_nfl_mimics <- csf_nfl %>% filter(Control %in% mdset)
csf_nfl_mimics <- csf_nfl_mimics[csf_nfl_mimics$Sensitivity_SE > 0, ]
csf_nfl_mimics %>% select(Study, Country, Region, Control) %>% arrange(Country)


### CSF_NfL vs NHC ------------------------------------------------

csf_nfl_hc <- csf_nfl %>% filter(Control %in% hcset)
csf_nfl_hc %>% select(Study, Country, Region, Control) %>% arrange(Country)


### CSF_NfL vs DC --------------------------------------------------

csf_nfl_dc <- csf_nfl %>% filter(Control %in% dcset) %>% 
  filter(!(Study=="Gagliardi 2021" & TargetBiomarker=="NfL")) #duplicate cohort; 2021 and 2024

csf_nfl_dc %>% select(Study, Country, Region, Control) %>% arrange(Country)


## CSF CHIT1 ----------------------------------------------------------

csf_chit1 <- data %>%
  filter(!is.na(Sensitivity) & !is.na(Specificity))%>%
  filter(FluidType == "CSF", TargetBiomarker == "CHIT1") %>% 
  select(all_of(SET)) %>% 
  mutate(
    Cut_off_Value = ifelse(is.na(Cut_off_Value), "not specified", Cut_off_Value)
  )


### CSF CHIT1 vs mimics ------------------------------------------------

CCM <- csf_chit1 %>% filter(Control %in% mdset)
CCM %>% select(Study, Country, Region, Control) %>% arrange(Country)

### CSF CHIT1 vs mddc 
# ## No additional eligible data were available for this comparison.
# CCMDC <- csf_chit1 %>% filter(Control %in% mddc)
# CCH %>% select(Study, Country, Region, Control) %>% arrange(Country)


### CSF CHIT1 vs HC (+1 DC mixed) -------------------------------------------

CCH <- csf_chit1 %>% filter(Control %in% hcset)
CCH %>% select(Study, Country, Region, Control) %>% arrange(Country)



## CSF YKL40 ---------------------------------------------------------------

csf_ykl40 <- data %>%
  filter(!is.na(Sensitivity) & !is.na(Specificity))%>%
  filter(FluidType == "CSF", TargetBiomarker == "YKL40") %>% 
  select(all_of(SET)) %>% 
  mutate(
    Cut_off_Value = ifelse(is.na(Cut_off_Value), "not specified", Cut_off_Value)
  )

### YKL40 vs mimics --------------------------------------------------------

CYM <- csf_ykl40 %>% filter(Control %in% mdset)
CYM %>% select(Study, Country, Region, Control) %>% arrange(Country)

### YKL40 mddc 

## No additional eligible data were available for this comparison.
# CYMDC <- csf_ykl40 %>% filter(Control %in% mddc)
# CYMDC %>% select(Study, Country, Region, Control) %>% arrange(Country)

### YKL40 HC (+1 DC mixed) --------------------------------------------------

CYC <- csf_ykl40 %>% filter(Control %in% hcset)
CYC %>% select(Study, Country, Region, Control) %>% arrange(Country)



## CSF pNfH --------------------------------------------------------

csf_pnfh <- data %>%
  filter(!is.na(Sensitivity) & !is.na(Specificity))%>%
  filter(FluidType == "CSF", TargetBiomarker == "pNfH") %>% 
  select(all_of(SET)) %>% 
  mutate(
    Cut_off_Value = ifelse(is.na(Cut_off_Value), "not specified", Cut_off_Value)
  )

### CSF pNfH vs ALS mimics -----------------------------------------

CPM <- csf_pnfh %>% filter(Control %in% mdset) 
CPM %>% select(Study, Country, Region, Control) %>% arrange(Country)

### CSF pNfH vs MDDC --------------------------------------------------------

CPD <- csf_pnfh %>% filter(Control %in% mddc) %>% 
  filter(!(Study=="Verde 2020" & Control=="PLS")) #duplicate cohort; Gagliardi 2021 (> Verde2020) Milan IRCCS
CPD %>% select(Study,Region,Control,n3) %>% arrange(Region)

### CSF pNfH vs HC --------------------------------------------------------

CPH <- csf_pnfh %>% filter(Control %in% hcset) %>% 
  filter(!(Study=="Verde 2020" & TargetBiomarker=="pNfH"))  #duplicate cohort; Gagliardi 2021 = Verde2020, pNfH
CPH %>% select(Study, Country, Region, Control) %>% arrange(Country)  

## CSF NfH ------------------------------------------------------------

csf_nfh <- data %>%
  filter(!is.na(Sensitivity) & !is.na(Specificity))%>%
  filter(FluidType == "CSF", TargetBiomarker == "NfH") %>% 
  select(all_of(SET)) %>% 
  mutate(
    Cut_off_Value = ifelse(is.na(Cut_off_Value), "not specified", Cut_off_Value)
  )

csf_nfh  %>% select(Study, Country, Region, Control) %>% arrange(Country)


### CSF NfH vs HC ---------------------------------------------------

#all
# csf_nfh %>% metaAUC()
# csf_nfh %>% fHSROC()

## no mimics, no DC!!

#HC
csf_nfh_hc <- csf_nfh %>% filter(Control %in% hcset) 
csf_nfh_hc %>% select(Study, Country, Region, Control) %>% arrange(Country)

## CSF t-tau --------------------------------------------------------

csf_ttau <- data %>%
  filter(!is.na(Sensitivity) & !is.na(Specificity))%>%
  filter(FluidType == "CSF", TargetBiomarker == "t-tau") %>% 
  select(all_of(SET)) %>% 
  mutate(
    Cut_off_Value = ifelse(is.na(Cut_off_Value), "not specified", Cut_off_Value)
  )


### CSF Ttau vs mimics ------------------------------------------

csf_ttau_mimics <- csf_ttau %>% filter(Control %in% mdset) 
csf_ttau_mimics  %>% select(Study, Country, Region, Control) %>% arrange(Country)

### CSF Ttau vs NHC ------------------------------------------

csf_ttau_hc <- csf_ttau %>% filter(Control %in% hcset)
csf_ttau_hc  %>% select(Study, Country, Region, Control) %>% arrange(Country)

### CSF Ttau vs mddc ------------------------------------------

csf_ttau_mddc <- csf_ttau %>% filter(Control %in% mddc)%>% 
  filter(!(Study=="Agnello 2021" & Control=="DC")) %>% #duplicate cohort
  filter(!(Study=="Ye 2020" & Control=="AD")) #duplicate cohort

csf_ttau_mddc  %>% select(Study, Country, Region, Control) %>% arrange(Country)


## CSF p-tau --------------------------------------------------------

csf_ptau <- data %>%
  filter(!is.na(Sensitivity) & !is.na(Specificity))%>%
  filter(FluidType == "CSF", TargetBiomarker == "p-tau181") %>% 
  select(all_of(SET)) %>% 
  mutate(
    Cut_off_Value = ifelse(is.na(Cut_off_Value), "not specified", Cut_off_Value)
  )

### CSF ptau vs MDDC: just 1 study

csf_ptau %>% filter(Control %in% mddc) 

### CSF ptau vs HCDC
csf_ptau_con <- csf_ptau %>% filter(Control %in% hcdc) %>% filter(!Control=="AD") # Ye 2020 duplication (AD, DC)
csf_ptau_con %>% select(Study, Country, Region, Control) %>% arrange(Country)


## CSF p/t tau ratio --------------------------------------------------------

csf_ptr <- data %>%
  filter(!is.na(Sensitivity) & !is.na(Specificity))%>%
  filter(FluidType == "CSF", TargetBiomarker == "ptau_ttau_ratio") %>% 
  select(all_of(SET)) %>% 
  mutate(
    Cut_off_Value = ifelse(is.na(Cut_off_Value), "not specified", Cut_off_Value)
  )

### CSF p/t tau ratio vs ALS mimics -----------------------------------------

csf_ptr_mimics <- csf_ptr %>% filter(Control %in% mdset)
csf_ptr_mimics %>% select(Study, Country, Region, Control) %>% arrange(Country)


### CSF p/t tau ratio vs Cons (HC and DC) -----------------------------------

csf_ptr_con <- csf_ptr %>% filter(Control %in% hcdc)
csf_ptr_con %>% select(Study, Country, Region, Control) %>% arrange(Country)


## No additional eligible data were available for this comparison.


## CSF MCP1 --------------------------------------------------------

csf_mcp1 <- data %>%
  filter(!is.na(Sensitivity) & !is.na(Specificity))%>%
  filter(FluidType == "CSF", TargetBiomarker == "MCP1") %>% 
  select(all_of(SET)) %>% 
  mutate(
    Cut_off_Value = ifelse(is.na(Cut_off_Value), "not specified", Cut_off_Value)
  )

### CSF MCP1 vs ALS mimics ------------------------------------------
csf_mcp1_mimics <- csf_mcp1 %>% filter(Control %in% mdset) 

csf_mcp1_mimics %>% select(Study, Country, Region, Control) %>% arrange(Country)

## no HCs, no additional MDDC!




## CSF TDP43 --------------------------------------------------------

csf_tdp43 <- data %>%
  filter(!is.na(Sensitivity) & !is.na(Specificity))%>%
  filter(FluidType == "CSF", TargetBiomarker == "TDP43") %>% 
  select(all_of(SET)) %>% 
  mutate(
    Cut_off_Value = ifelse(is.na(Cut_off_Value), "not specified", Cut_off_Value)
  )

csf_tdp43 %>% select(Study, Country, Region, Control) %>% arrange(Country)

# no ALS mimics!, no DCs!

### CSF TDP43 vs NHC

csf_tdp43_hc <- csf_tdp43 %>% filter(Control %in% hcset) 
csf_tdp43_hc %>% select(Study, Country, Region, Control) %>% arrange(Country)


## CSF Ab42 --------------------------------------------------------

## Only Ye 2020

# csf_Ab42 <- data %>%
#   filter(!is.na(Sensitivity) & !is.na(Specificity))%>%
#   filter(FluidType == "CSF", TargetBiomarker == "Aβ42") %>%
#   select(all_of(SET)) %>%
#   mutate(
#     Cut_off_Value = ifelse(is.na(Cut_off_Value), "not specified", Cut_off_Value)
#   )
# 
# csf_Ab42



## CSF CXCL12 --------------------------------------------------------

## Same study group only:
## Andrés-Benito 2020 and Roca-Pereira 2024
## Spain (Barcelona), Bellvitge University Hospital (UFELA, Neurology Service)

# csf_CXCL12 <- data %>%
#   filter(!is.na(Sensitivity) & !is.na(Specificity))%>%
#   filter(FluidType == "CSF", TargetBiomarker == "CXCL12") %>%
#   select(all_of(SET)) %>%
#   mutate(
#     Cut_off_Value = ifelse(is.na(Cut_off_Value), "not specified", Cut_off_Value)
#   )
# 
# csf_CXCL12 %>% select(Study, Country, Region, Control)



# Blood analysis ----------------------------------------------------------

## Blood NfL ---------------------------------------------------------------

b_nfl <- data %>%
  filter(!is.na(Sensitivity) & !is.na(Specificity))%>%
  filter(FluidType == "blood" & TargetBiomarker == "NfL") %>% 
  select(all_of(SET)) %>% 
  mutate(
    Cut_off_Value = ifelse(is.na(Cut_off_Value), "not specified", Cut_off_Value)
  )

### blood NfL vs mimics ------------------------------

BNM <- b_nfl %>% 
  filter(Control %in% mdset) %>% 
  filter(!(AssayMethod == "MSD assay, serum" & Control == "HSP")) %>% #duplicate cohort; Gille 2019-2
  filter(!(AssayMethod == "MSD assay, serum" & Control == "PMA")) %>% #duplicate cohort; Gille 2019-2
  filter(!(AssayMethod == "MSD assay, serum" & Control == "PLS")) %>% #duplicate cohort; Gille 2019-2 
  distinct(Study, TargetBiomarker, Control, n1, n2, .keep_all = TRUE) # Mondesert 2025, Simoa

BNM %>% select(Study, Country, Region, Control) %>% arrange(Country)


### blood NfL vs DC ----------------------------------------------------

BND <- b_nfl %>% filter(Control %in% dcset) %>% 
  filter(!(AssayMethod=="Ella, serum" & Study=="Brousse 2023")) %>% #duplicate cohort (Ella & Simoa)
  filter(!(Control=="NonALS" & Study =="Verde 2019")) #duplicate cohort; Verde 2019, vs DC

BND %>% select(Study, Country, Region, Control) %>% arrange(Country)


### blood NfL vs HC ----------------------------------------------------

BNH <- b_nfl %>% filter(Control %in% hcset) %>% 
  filter(!(Study=="Verde 2019")) %>%  #duplicate cohort; Verde 2023-1
  filter(!(str_starts(Study, "Halbgebauer 2022"))) #duplicate cohort; Ulm cohort (Halbgebauer 2022–2  and Witzel 2024)
BNH %>% select(Study, Country, Region, Control) %>% arrange(Country)

## Note) Halbgebauer2022-2 may be included in Witzel 2024 (Ulm Cohort)


## Blood pNfH --------------------------------------------------------

blood_pnfh <- data %>%
  filter(!is.na(Sensitivity) & !is.na(Specificity))%>%
  filter(FluidType == "blood", TargetBiomarker == "pNfH") %>% 
  select(all_of(SET)) %>% 
  mutate(
    Cut_off_Value = ifelse(is.na(Cut_off_Value), "not specified", Cut_off_Value)
  )

## no vs mimics!
# blood_pnfh %>% filter(Control %in% mdset)

### blood pNfH vs HC  ----------------------------------------------

blood_pnfh_hc <- blood_pnfh %>% filter(Control %in% hcset) 
blood_pnfh_hc %>% select(Study, Country, Region, Control) %>% arrange(Country)

### blood pNfH vs DC ----------------------------------------------
blood_pnfh_dc <- blood_pnfh %>% filter(Control %in% mddc) 
blood_pnfh_dc %>% select(Study, Country, Region, Control) %>% arrange(Country)


## Blood NfH --------------------------------------------------------

blood_NfH <- data %>%
  filter(!is.na(Sensitivity) & !is.na(Specificity))%>%
  filter(FluidType == "blood", TargetBiomarker == "NfH") %>% 
  select(all_of(SET)) %>% 
  mutate(
    Cut_off_Value = ifelse(is.na(Cut_off_Value), "not specified", Cut_off_Value)
  )

### all controls (HCDC)  ----------------------------------------------

blood_NfH_con <- blood_NfH %>% filter(Control %in% hcset) 
blood_NfH_con %>% select(Study, Country, Region, Control) %>% arrange(Country)


## blood GFAP --------------------------------------------------------

blood_GFAP <- data %>%
  filter(!is.na(Sensitivity) & !is.na(Specificity))%>%
  filter(FluidType == "blood", TargetBiomarker == "GFAP") %>% 
  select(all_of(SET)) %>% 
  mutate(
    Cut_off_Value = ifelse(is.na(Cut_off_Value), "not specified", Cut_off_Value)
  )

### no dc, no md
# blood_GFAP %>% metaAUC()

### blood GFAP vs NHC ------------------------------------------------
blood_GFAP_hc <- blood_GFAP %>% filter(Control %in% hcset) 
blood_GFAP_hc %>% select(Study,Country,Region)


## blood UCHL1 --------------------------------------------------------

blood_UCHL1 <- data %>%
  filter(!is.na(Sensitivity) & !is.na(Specificity))%>%
  filter(FluidType == "blood", TargetBiomarker == "UCHL1") %>% 
  select(all_of(SET)) %>% 
  mutate(
    Cut_off_Value = ifelse(is.na(Cut_off_Value), "not specified", Cut_off_Value)
  )

### Sanity check
#blood_UCHL1 %>% select(Study, Country, Region, Control) %>% arrange(Country)

### serum UCHL1 vs mimics, DC --------------------------------------------

blood_UCHL1_mddc <- blood_UCHL1 %>% 
  filter(Control %in% mddc) %>% 
  filter(!(Study=="Falzone 2022" & Control=="DC")) #duplicate cohort

blood_UCHL1_mddc %>% select(Study, Country, Region, Control) %>% arrange(Country)

#no HC


# 📊 Summary Table -----------------------------------------------------------------

DiagTable <- list(
  
  ### blood analysis set
  BNM = list(type = "Blood", name = "NfL", comparison = "ALS mimics"),
  BNH = list(type = "Blood", name = "NfL", comparison = "NHC"),
  BND = list(type = "Blood", name = "NfL", comparison = "DC"),                 　# pure DC
  blood_pnfh_hc = list(type = "Blood", name = "pNfH", comparison = "NHC"),       # pure NHC
  blood_pnfh_dc = list(type = "Blood", name = "pNfH", comparison = "DC"),        # pure DC
  blood_NfH_con = list(type = "Blood", name = "NfH", comparison = "NHC§"),       # incl.DC (1 cohort)
  blood_UCHL1_mddc = list(type = "Blood", name = "UCHL1", comparison = "DC¶"),   # incl.MD (MD1, DC1)
  blood_GFAP_hc = list(type = "Blood", name = "GFAP", comparison = "NHC"),       # pure NHC
  
  ### CSF analysis set
  csf_nfl_mimics = list(type = "CSF", name = "NfL", comparison = "ALS mimics"),   # pure MD
  csf_nfl_dc = list(type = "CSF", name = "NfL", comparison = "DC"),               # pure DC
  csf_nfl_hc = list(type = "CSF", name = "NfL", comparison = "NHC"),              # pure NHC
  CPM = list(type = "CSF", name = "pNfH", comparison = "ALS mimics"),             # pure MD
  CPD = list(type = "CSF", name = "pNfH", comparison = "DC¶"),                    # incl.MD (mddc)
  CPH = list(type = "CSF", name = "pNfH", comparison = "NHC"),                    # pure NHC
  CCM = list(type = "CSF", name = "CHIT1", comparison = "ALS mimics"),            # pure MD
  CCH = list(type = "CSF", name = "CHIT1", comparison = "NHC§"),                  # incl.DC (1 cohort) **
  CYM = list(type = "CSF", name = "YKL40", comparison = "ALS mimics"),            # pure MD
  CYC = list(type = "CSF", name = "YKL40", comparison = "NHC§"),                  # incl.DC (1 cohort)
  csf_nfh_hc = list(type = "CSF", name = "NfH", comparison = "NHC"),               # pure NHC
  csf_mcp1_mimics = list(type = "CSF", name = "MCP1", comparison = "ALS mimics"),  # pure MD
  csf_ptr_mimics = list(type = "CSF", name = "p-tau/t-tau", comparison = "ALS mimics"), # pure MD
  csf_ptr_con = list(type = "CSF", name = "p-tau/t-tau", comparison = "NHC§"),     # incl.DC (1 cohort)
  csf_ttau_mimics = list(type = "CSF", name = "t-tau", comparison = "ALS mimics"), # pure mimics
  csf_ttau_hc = list(type = "CSF", name = "t-tau", comparison = "NHC"),            # pure NHC
  csf_ttau_mddc = list(type = "CSF", name = "t-tau", comparison = "DC¶"),          # incl.MD
  csf_ptau_con = list(type = "CSF", name = "p-tau181", comparison = "NHC§"),       # HC1 DC1
  csf_tdp43_hc = list(type = "CSF", name = "TDP43", comparison = "NHC")            # pure NHC
)


## Summarise cohorts included in the meta-analysis (HSROC sAUC and random-effects AUC)
## Results container
results <- list()

## Build the summary table
for (dataset_name in names(DiagTable)) {
  dataset <- get(dataset_name) # obtain the dataset
  type <- DiagTable[[dataset_name]]$type
  name <- DiagTable[[dataset_name]]$name
  comparison <- DiagTable[[dataset_name]]$comparison
  results[[dataset_name]] <- analyze_dataset(dataset, type, name, comparison)
}

# final results
final_results <- do.call(rbind, results)
tbl <- format_results(final_results) %>% select(-"AUC (phm model)")
#tbl %>% gt()


# Visualization Table 
tbl_vs_HC <- tbl %>% filter(Comparison %in% c("NHC", "NHC§")) %>% select(-Controls)
tbl_vs_mimics <- tbl %>% filter(Comparison == "ALS mimics")%>% select(-Controls)

# HC and mimics
tbl_merged <- full_join(tbl_vs_HC, tbl_vs_mimics, by = c("Type", "Biomarker"), 
                        suffix = c("", "*")) %>% 
  mutate(Biomarker = factor(Biomarker, levels = c(bmsort, setdiff(Biomarker, bmsort)))) %>% 
  arrange(Biomarker)

# Visual Check
# tbl_merged %>% gtROC() %>% tab_style(
#     style = cell_borders(sides = "left", weight = px(1), color="#D3D3D3", style="dashed"),
#     locations = cells_body(columns = c("Comparison","Comparison*"))) 


# Table 2 for Publication ---------------------------------------------------------------

Table2 <- tbl_merged %>% 
  select(-Cohorts_AUC,-`Cohorts_AUC*`,-`N, ALS_auc`,-`N, Con_auc`,-`N, ALS_auc*`,-`N, Con_auc*`) %>% 
  create_grouped_flextable(group_col = "Type") %>%
  set_header_labels(
    `Pooled Se.`            = "Summary Se.*",
    `Pooled Sp.`            = "Summary Sp.*",
    `Pooled Se.*`           = "Summary Se.*",
    `Pooled Sp.*`           = "Summary Sp.",
    `AUC (RE model)`        = "AUC (RE model)‡",
    `AUC (HSROC model)`     = "sAUC†",
    `Comparison*`           = "Comparison",
    `Cohorts*`              = "Cohorts",
    `N, ALS*`               = "N, ALS",
    `N, Con*`               = "N, Con",
    `AUC (HSROC model)*`    = "sAUC†",
    `AUC (RE model)*`       = "AUC (RE model)‡"
  ) %>%
  add_footer_lines(values = c(
    "*) Summary point on the HSROC curve (95% confidence interval).",
    "†) Summary AUC from the HSROC model (sAUC).",
    "‡) Random-effects model estimate (95% confidence interval).",
    "§) One study combines neurological disease controls with neurologically healthy controls."
  )) %>%
  align(align = "left", part = "footer") %>% 
  flextable::fontsize(size = 10, part = "all") %>%   # make table more compact
  flextable::font(fontname = "Times New Roman", part = "all") %>% 
  width(j = 1, width = 0.65) %>% 
  width(j = 2:5, width = 0.55) %>% 
  width(j = 6:9, width = 0.75) %>% 
  width(j = 10:13, width = 0.55) %>% 
  width(j = 14:17, width = 0.75)  #%>% ExFT(output/Table2.docx)
