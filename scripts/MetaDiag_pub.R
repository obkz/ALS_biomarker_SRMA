source("scripts/MetaFunc.R")
source("scripts/data_prep/DiagData.R")

# CSF ---------------------------------------------------------------------

## CSF NfL - mimics --------------------------------------------------------

csf_nfl_mimics %>% select(Study, Country, Region, Control) %>% arrange(Country)
csf_nfl_mimics %>% metaAUC()
fHSROC(csf_nfl_mimics)
DiagStbl(csf_nfl_mimics, title="CSF NfL, vs ALS mimics, Summary Table")


## CSF NfL - HC ------------------------------------------------------------

csf_nfl_hc %>% select(Study, Country, Region, Control) %>% arrange(Country)
csf_nfl_hc %>% metaAUC()
fHSROC(csf_nfl_hc)


## CSF NfL - DC ------------------------------------------------------------

csf_nfl_dc %>% select(Study, Country, Region, Control) %>% arrange(Country)
csf_nfl_dc %>% metaAUC()
fHSROC(csf_nfl_dc)


## CSF CHIT1 vs mimics ------------------------------------------------

CCM %>% select(Study, Country, Region, Control) %>% arrange(Country)
CCM %>% metaAUC()
CCM %>% fHSROC()

## CSF CHIT1 vs HC ----------------------------------------------------

CCH %>% select(Study, Country, Region, Control) %>% arrange(Country)
CCH %>% metaAUC()
CCH %>% fHSROC()


## CSF YKL40 vs mimics --------------------------------------------------------

CYM %>% select(Study, Country, Region, Control) %>% arrange(Country)
CYM %>% metaAUC()
CYM %>% fHSROC()

## CSF YKL40 HC --------------------------------------------------------

CYC %>% select(Study, Country, Region, Control) %>% arrange(Country)
CYC %>% metaAUC()
CYC %>% fHSROC()


## CSF pNfH vs ALS mimics -----------------------------------------

CPM %>% select(Study, Country, Region, Control) %>% arrange(Country)
CPM %>% metaAUC()
CPM %>% fHSROC()

## CSF pNfH vs MDDC --------------------------------------------------------

CPD %>% select(Study,Region,n3) %>% arrange(Region)
CPD %>% metaAUC()
CPD %>% fHSROC()

## CSF pNfH vs HC --------------------------------------------------------

CPH %>% select(Study, Country, Region, Control) %>% arrange(Country)  
CPH %>% metaAUC()
CPH %>% fHSROC()


## CSF NfH vs HC------------------------------------------------------------

csf_nfh_hc %>% metaAUC()
csf_nfh_hc %>% fHSROC()


## CSF Ttau vs Mimics -----------------------------------------------------

csf_ttau_mimics %>% select(Study, Country, Region, Control) %>% arrange(Country)
csf_ttau_mimics %>% metaAUC()
csf_ttau_mimics %>% fHSROC()

## CSF Ttau vs NHC -----------------------------------------------------
csf_ttau_hc  %>% select(Study, Country, Region, Control) %>% arrange(Country)

csf_ttau_hc %>% metaAUC()
csf_ttau_hc %>% fHSROC()

## CSF Ttau vs mddc -----------------------------------------------------

csf_ttau_mddc  %>% select(Study, Country, Region, Control) %>% arrange(Country)

csf_ttau_mddc %>% metaAUC()
csf_ttau_mddc %>% fHSROC()


## CSF p-tau --------------------------------------------------------

csf_ptau_con %>% select(Study, Country, Region, Control) %>% arrange(Country)

csf_ptau_con %>% metaAUC()
csf_ptau_con %>% fHSROC()


## CSF p/t tau ratio vs mimics -----------------------------

csf_ptr_mimics %>% select(Study, Country, Region, Control) %>% arrange(Country)
csf_ptr_mimics %>% metaAUC()
csf_ptr_mimics %>% fHSROC()

## CSF p/t tau ratio vs Cons -----------------------------

csf_ptr_con %>% select(Study, Country, Region, Control) %>% arrange(Country)
csf_ptr_con %>% metaAUC()
csf_ptr_con %>% fHSROC()

## NO MDDC (no additional data)


## CSF MCP1 vs ALS mimics -----------------------------------------

csf_mcp1_mimics %>% select(Study, Country, Region, Control) %>% arrange(Country)
csf_mcp1_mimics %>% metaAUC()
csf_mcp1_mimics %>% fHSROC()

## no HCs, no additional MDDC!


## CSF TDP43 vs HC --------------------------------------------------------

#no ALS mimics!, no DCs!

csf_tdp43_hc %>% select(Study, Country, Region, Control) %>% arrange(Country)
csf_tdp43_hc %>% metaAUC()
csf_tdp43_hc %>% fHSROC()


## CSF Ab42 --------------------------------------------------------

## Only Ye 2020


## CSF CXCL12 --------------------------------------------------------

## the same Study group, only
## Andrés-Benito 2020 & Roca-Pereira 2024
## Spain, Barcerona, Bellvitge University Hospital (UFELA of the Neurology Service)



# Blood ----------------------------------------------------------

## Blood NfL mimics------------------------------------------------------

BNM %>% select(Study, Country, Region, Control) %>% arrange(Country)
BNM %>% metaAUC()
BNM %>% fHSROC()

## blood NfL vs DC ----------------------------------------------------

BND %>% select(Study, Country, Region, Control) %>% arrange(Country)
BND %>% metaAUC()
BND %>% fHSROC()

## blood NfL vs HC ----------------------------------------------------

BNH %>% select(Study, Country, Region, Control) %>% arrange(Country)
BNH %>% metaAUC()
BNH %>% fHSROC()

DiagStbl(BNH)

## blood pNfH vs HC -----------------------------------------------

blood_pnfh_hc %>% select(Study, Country, Region, Control) %>% arrange(Country)
blood_pnfh_hc %>% metaAUC()
blood_pnfh_hc %>% fHSROC()

## blood pNfH vs DC -----------------------------------------------

blood_pnfh_dc %>% select(Study, Country, Region, Control) %>% arrange(Country)
blood_pnfh_dc %>% metaAUC() #Shi et al did not provide AUC, but provide Sens/Spec (==> bivariate Reitsma only)
blood_pnfh_dc %>% fHSROC()

## Blood NfH  (HC only) -----------------------------------------------

blood_NfH_con %>% select(Study, Country, Region, Control) %>% arrange(Country)
blood_NfH_con %>% metaAUC()
blood_NfH_con %>% fHSROC()

## blood GFAP vs NHC -----------------------------------------------

blood_GFAP_hc %>% select(Study,Country,Region,Control)
blood_GFAP_hc %>% metaAUC()
blood_GFAP_hc  %>% fHSROC()

## no DC, no mimics

## serum UCHL1 vs mimics, DC -----------------------------------------------

blood_UCHL1_mddc %>% select(Study,Country,Region,Control)
blood_UCHL1_mddc %>% metaAUC()
blood_UCHL1_mddc %>% fHSROC()


#no HC


# 📊 Summary Table ----------------------------------------------------------



# Visualization Table --------------------------------------------------------

tbl_vs_HC <- tbl %>% filter(Comparison %in% c("NHC", "NHC§")) %>% select(-Controls)
tbl_vs_mimics <- tbl %>% filter(Comparison == "ALS mimics")%>% select(-Controls)

tbl_merged <- full_join(tbl_vs_HC, tbl_vs_mimics, by = c("Type", "Biomarker"), 
                        suffix = c("", "*")) %>% 
  mutate(Biomarker = factor(Biomarker, levels = c(bmsort, setdiff(Biomarker, bmsort)))) %>% 
  arrange(Biomarker)

# 
# tbl1a <- tbl_merged %>% filter(Type=="Blood") %>% flextable()
# tble1b <- tbl_merged %>% filter(Type=="CSF") %>% flextable()
# 
# df_test <- bind_rows(
#   tbl_merged %>% filter(Type == "Blood"),
#   tbl_merged %>% filter(Type == "CSF")
# ) %>%arrange(Type)
# 
# df_test %>% flextable() %>% merge_v(j="Type") %>%  align(j = "Type", align = "left", part = "body") %>% autofit()
# 
# ?gt::as_word()

tbl_merged %>% gtROC() %>% 
  tab_style(
    style = cell_borders(sides = "left", weight = px(1), color="#D3D3D3", style="dashed"),
    locations = cells_body(columns = c("Comparison","Comparison*"))
  ) 


# Table 2 for Publication ---------------------------------------------------------------

create_grouped_flextable(tbl_merged, group_col = "Type") %>%
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
    "*) Summary Points on HSROC Curve, 95% Confidence Interval",
    "†) Summary AUC of HSROC Curve",
    "‡) Random Effect model, 95% Confidence Interval",
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


# ## remove RE model
# tbl_merged %>% select(-`AUC (RE model)`,-`AUC (RE model)*`) %>% 
#   create_grouped_flextable(group_col = "Type") %>%
#   set_header_labels(
#     `Pooled Se.*`           = "Summary Se.*",
#     `Pooled Sp.*`           = "Summary Sp.",
#     `Pooled Se.`            = "Summary Se.*",
#     `Pooled Sp.`            = "Summary Sp.*",
#     `AUC (HSROC model)`     = "sAUC†",
#     `Comparison*`           = "Comparison",
#     `Cohorts*`              = "Cohorts",
#     `N, ALS*`               = "N, ALS",
#     `N, Con*`               = "N, Con",
#     `AUC (HSROC model)*`    = "sAUC†"
#   ) %>%
#   add_footer_lines(values = c(
#     "*) Summary Points on HSROC Curve, 95% Confidence Interval",
#     "†) Summary AUC of HSROC Curve",
#     "‡) Random Effect model, 95% Confidence Interval"
#   )) %>%
#   align(align = "left", part = "footer") %>% 
#   flextable::fontsize(size = 10, part = "all") %>%   # make table more compact
#   flextable::font(fontname = "Times New Roman", part = "all") %>% 
#   width(j = 1:5, width = 0.55) %>% 
#   width(j = 6:8, width = 0.75) %>% 
#   width(j = 9:12, width = 0.55) %>% 
#   width(j = 13:15, width = 0.75) 


# Notes -----------------------------------------------------------

tbl %>% gt()

### Arrange
gtbl <- tbl %>% 
  group_by(Biomarker) %>%  # Group by Biomarker
  mutate(total_Cohorts = sum(Cohorts)) %>%  # Calculate total Cohorts per Biomarker
  ungroup() %>%  # Remove grouping temporarily to sort the entire dataset
  arrange(Type,desc(total_Cohorts), Biomarker) %>%  # Sort by total Cohorts and Biomarker name
  arrange(Type,Comparison,desc(Cohorts),desc(`AUC (HSROC model)` ))  %>% 
  select(-total_Cohorts)  #%>%select(-`AUC (RE model)`)

### All
#gtbl %>% gtt() %>% gtsave("output/test.html")

### ALS mimics
gtbl %>% filter(Comparison=="ALS mimics") %>% gtt(groupname_col = c("Type","Comparison"))

gtbl %>% filter(Comparison=="ALS mimics") %>% 
  filter(!(Biomarker == "NfL" & Type %in% c("Serum", "Plasma"))) %>% 
  mutate(Type = ifelse(Type %in% c("Serum", "Plasma"), "Blood", Type)) %>% 
  gtt(groupname_col = c("Type","Comparison"))


### HC
gtbl %>% filter(Comparison %in% c("HC")) %>% gtt(groupname_col = c("Type","Comparison"))

gtbl %>% filter(Comparison=="HC") %>% 
  filter(!(Biomarker == "NfL" & Type %in% c("Serum", "Plasma"))) %>% 
  mutate(Type = ifelse(Type %in% c("Serum", "Plasma"), "Blood", Type)) %>% 
  gtt(groupname_col = c("Type","Comparison"))



### DCs
gtbl %>% filter(Comparison %in% c("DC", "NHC§")) %>% gtt(groupname_col = c("Type","Comparison"))
gtbl %>% filter(Comparison %in% c("DC", "NHC§"), Type %in% c("Blood","CSF")) %>% gtt(groupname_col = c("Type","Comparison"))

#### all, grouped by biomarkers
gtbl %>% 
  gtt(groupname_col = "Biomarker")

#### all, grouped by biomarkers (excluding the biomarkers reported in only 1 comparison)
gtbl %>% 
  filter(!Biomarker %in% c("MCP1", "UCHL1","GFAP","TDP43","p-tau","p-tau/t-tau")) %>% 
  gtt(groupname_col = "Biomarker")

#### NfL
gtbl %>% 
  filter(Biomarker %in% c("NfL")) %>% 
  gtt(groupname_col = c("Type","Biomarker"))

#### PNfH, NfH
gtbl %>% 
  filter(Biomarker %in% c("pNfH","NfH")) %>% 
  gtt(groupname_col = "Biomarker")

#### Tau 
gtbl %>% 
  filter(Biomarker %in% c("T-tau","p-tau","p-tau/t-tau")) %>% 
  arrange(Biomarker) %>% 
  gtt(groupname_col = c("Type","Biomarker"))

#### CHIT
gtbl %>% 
  filter(Biomarker %in% c("CHIT1","YKL40")) %>% 
  gtt(groupname_col = c("Type","Biomarker"))







