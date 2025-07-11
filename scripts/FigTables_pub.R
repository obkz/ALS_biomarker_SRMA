source("scripts/MetaFunc.R")
source("scripts/data_prep/DiagData.R")

#meta
# source("scripts/MetaDiag_pub.R")
# source("scripts/MetaAUC_pub.R")
# source("scripts/MetaSMD_pub.R")
# source("scripts/MetaHR_pub.R")
# source("scripts/MetaCorr_pub.R")

# Main Figs & Tables ------------------------------------------------------

## Fig.1 Diagnostic Performance for Differentiating ALS from ALS Mimics --------

## Fig.1a

#CairoPDF(file= "output/AUC/AUC_blood_mimics.pdf", width=14, height=7)
#svg(file= "output/AUC/AUC_blood_mimics.svg", width=14, height=7)
bloodAUCmd %>% metaAUC(xlab = "AUC (95%CI), blood biomarkers, vs ALS mimics")
#dev.off()

## Fig.1b

#CairoPDF(file= "output/AUC/AUC_CSF_mimics.pdf", width=14, height=16)
#svg(file= "output/AUC/AUC_CSF_mimics.svg", width=14, height=16)
csfAUCmd %>% metaAUC(xlab = "AUC (95%CI), CSF biomarkers, vs ALS mimics", subname = "CSF")
#dev.off()

## Fig.1c
## vs mimics
#svg(file= "output/AUC/hsroc_compare.svg", width=7, height=7)
fHSROC_compare(csf_nfl_mimics, BNM, 
               plot_title = "SROC curves of blood and CSF NfL (vs ALS mimics)",
               label1="CSF NfL", label2="Blood NfL",
               col1="navy",   col2="darkgreen",
               pch1=19,       pch2=17)
#dev.off()

## vs NHC 
# fHSROC_compare(csf_nfl_hc, BNH, 
#                label1="CSF NfL", label2="Blood NfL",
#                col1="navy",   col2="darkgreen",
#                pch1=19,       pch2=17)


## Fig.2, Prognostic Performance Based on Multivariable Cox Models --------------

#CairoPDF(file= "output/HR/HR_blood_multivariable.pdf", width=16, height=14)
#svg(file= "output/HR/HR_blood_multivariable.svg", width=16, height=14)
BloodMultiHR %>% metaHR(pattern = 1, xlab="Multivariable HR (95%CI), Blood Biomarkers, dichotomized comparison")
#dev.off()

#CairoPDF(file= "output/HR/HR_CSF_multivariable0623.pdf", width=14, height=10)
#svg(file= "output/HR/HR_CSF_multivariable.svg", width=14, height=10)
CSFMultiHR %>% metaHR(pattern = 1, xlab="Multivariable HR (95%CI), CSF Biomarkers, dichotomized comparison", subname = "CSF")
#dev.off()


## Fig.3, Hypothetical Biomarker-informed Workflow for ALS ----------------
## BioRender


## Table 1 - overview -----------------------------------------------------
## MS word

## Table 2 - Summary of pooled diagnostic performance ---------------------
Table2

## Table 3 - Summary of pooled Cox Hazard ratios --------------------------
Table3

## Table 4 - Summary of pooled Corr coefficients --------------------------
Table4


# Supple --------------------------------------------------------------------

## SuppleFig.1 - PRISMA Flowchart---------------------------------------------
## MS Word

## SuppleFig.2 - Pooled AUCs (ALS vs NHC) ------------------------------------

#CairoPDF(file= "output/AUC/AUC_blood_NHC.pdf", width=14, height=9)
#svg(file= "output/supple/SuppleFig2a.svg", width=14, height=9)
bloodAUChc %>% metaAUC(xlab = "AUC (95%CI), blood biomarkers, vs NHC")
#dev.off()

#CairoPDF(file= "output/AUC/AUC_CSF_NHC.pdf", width=14, height=16)
#svg(file= "output/supple/SuppleFig2b.svg", width=14, height=16)
csfAUChc %>% metaAUC(xlab = "AUC (95%CI), CSF biomarkers, vs NHC",subname="CSF")
#dev.off()

## SuppleFig.3 - Pooled AUCs (ALS vs DC & mimics) ----------------------------

#CairoPDF(file= "output/AUC/AUC_blood_DC_mimics.pdf", width=14, height=10)
#svg(file= "output/supple/SuppleFig3a.svg", width=14, height=10)
bloodAUCmddc %>% metaAUC(xlab = "AUC (95%CI), blood biomarkers, vs DC including ALS mimics")
#dev.off()

#CairoPDF(file= "output/AUC/AUC_CSF_DC_and_mimics.pdf", width=14, height=17)
#svg(file= "output/supple/SuppleFig3b.svg", width=14, height=17)
csfAUCmddc %>% metaAUC(xlab = "AUC (95%CI), CSF biomarkers, vs DC including ALS mimics", subname="CSF")
#dev.off()


## SuppleFig.4 - QUADAS-2 ----------------------------------------------------

QUADASp

## SuppleFig.5 - Pooled SMDs (ALS vs mimics) ---------------------------------

#CairoPDF(file = "output/SMD/SMD_CSF_mimics.pdf", width=13, height=13, family = "MyHelv")
#svg(file = "output/supple/SuppleFig5a.svg", width=13, height=13)
smdCSFmd %>% mutate(Control = str_replace_all(Control, "ALSmimics", "mimics")) %>% 
  SMDplot(pattern = 2, Fluid="CSF", Con = "ALS mimics", filterN=2,
          xlab="CSF biomarkers SMD (ALS vs mimics), 95%CI", subname = "CSF")
#dev.off()
#dev.new()

#CairoPDF(file= "output/SMD/SMD_blood_mimics.pdf", width=14, height=7)
#svg(file = "output/supple/SuppleFig5b.svg", width=14, height=6)
smdBloodmd %>% mutate(Control = str_replace_all(Control, "ALSmimics", "mimics")) %>% 
  SMDplot(pattern = 2, Fluid="blood", Con = "ALS mimics", filterN=2,
          xlab="Blood biomarkers SMD (ALS vs mimics), 95%CI")
#dev.off()


## SuppleFig.6 - Pooled SMDs (ALS vs NHC) ------------------------------------

#CairoPDF(file= "output/SMD/SMD_blood_HC.pdf", width=18, height=34)
#svg(file = "output/supple/SuppleFig6a.svg", width=18, height=34)
smdBloodHC %>% SMDplot(pattern = 2, Fluid="blood", Con = "HC", filterN=2,
                       xlab="Blood biomarkers SMD (ALS vs NHC), 95%CI")
#dev.off()

#CairoPDF(file = "output/SMD/SMD_CSF_HC.pdf", width=13, height=24)
#svg(file = "output/supple/SuppleFig6b.svg", width=13, height=24)
smdCSFhc %>% SMDplot(pattern = 2, Fluid="CSF", Con = "HC", filterN=2,
                     xlab="CSF biomarkers SMD (ALS vs HC), 95%CI", subname = "CSF")
#dev.off()


## SuppleFig.7 - Pooled SMDs (ALS vs DC/mimics) ------------------------------

#CairoPDF(file= "output/SMD/SMD_blood_DC_and_mimics.pdf", width=14, height=15)
#svg(file = "output/supple/SuppleFig7a.svg", width=14, height=15)
smdBloodmddc %>% mutate(Control = str_replace_all(Control, "ALSmimics", "mimics")) %>% 
  SMDplot(pattern = 2, Fluid="blood", filterN=2,
                         xlab="Blood biomarkers SMD (ALS vs DC, including ALS mimics), 95%CI")
#dev.off()

#CairoPDF(file= "output/SMD/SMD_CSF_DC_and_mimics.pdf", width=14, height=18)
#svg(file = "output/supple/SuppleFig7b.svg", width=14, height=18)
smdCSFmddc %>% 
  mutate(Control = str_replace_all(Control, "ALSmimics", "mimics")) %>% 
  mutate(TargetBiomarker = str_replace(TargetBiomarker, "^ptau_ttau_ratio$", "p-tau/t-tau")) %>%
  SMDplot(pattern = 2, Fluid="CSF", filterN=2, subname = "CSF",
          xlab="CSF biomarkers SMD (ALS vs DC, including ALS mimics), 95%CI",)
#dev.off()


## SuppleFig.8 - Pooled HRs, univariable -------------------------------------

#CairoPDF(file= "output/HR/HR_blood_univariable.pdf", width=16, height=8)
#svg(file = "output/supple/SuppleFig8a.svg", width=16, height=8)
BloodUniHR %>% metaHR(pattern = 1, #subname = "blood"
                      xlab="Univariable HR (95%CI), Blood Biomarkers, dichotomized comparison") 
#dev.off()

#CairoPDF(file= "output/HR/HR_CSF_univariable.pdf", width=14, height=8)
#svg(file = "output/supple/SuppleFig8b.svg", width=14, height=8)
CSFUniHR %>% metaHR(pattern = 1, subname = "CSF",
                    xlab="Univariable HR (95%CI), CSF Biomarkers, dichotomized comparison")
#dev.off()


## SuppleFig.9 - QUIPS -------------------------------------------------------

QUIPSp

## SuppleFig.10 - Pooled rho (DPR) -------------------------------------------

#CairoPDF(file= "output/Corr/DPR_blood.pdf", width=14, height=18)
#svg(file = "output/supple/SuppleFig10a.svg", width=14, height=18)
bloodCORdpr %>% mcforest(sm="COR", #subname="blood",
                         xlab="Correlation of blood biomarkers with DPR (Spearman's ρ, 95% CI)")
#dev.off()

#CairoPDF(file= "output/Corr/DPR_CSF.pdf", width=14, height=18)
#svg(file = "output/supple/SuppleFig10b.svg", width=14, height=18)
csfCORdpr %>% mcforest(sm="COR", subname="CSF", 
                       xlab="Correlation of CSF biomarkers with DPR (Spearman's ρ, 95% CI)") 
#dev.off()

## SuppleFig.11 - Pooled rho (FRS) -------------------------------------------

#CairoPDF(file= "output/Corr/ALSFRS_blood.pdf", width=14, height=21)
#svg(file = "output/supple/SuppleFig11a.svg", width=14, height=21)
bloodCorFRS %>% mcforest(sm="COR", cor="Corr_ALSFRSR",
                         xlab="Correlation of Blood biomarkers with ALSFRS (Spearman's ρ, 95% CI)")
#dev.off()

#CairoPDF(file= "output/Corr/ALSFRS_CSF.pdf", width=14, height=14)
#svg(file = "output/supple/SuppleFig11b.svg", width=14, height=21)
csfCorFRS %>% mcforest(sm="COR", cor="Corr_ALSFRSR",subname = "CSF",
                       xlab="Correlation of CSF biomarkers with ALSFRS (Spearman's ρ, 95% CI)")
#dev.off()


# Tables ------------------------------------------------------------------

## Supple Table 1 - Summary of pooled Cox Hazard ratios --------------------

#sT1
suppleT1 


## Supple Table 2 - Covariates ---------------------------------------------

sT2a
sT2b

