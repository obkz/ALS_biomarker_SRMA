# Functions for meta-analysis
# Load libraries ----------------------------------------------------------
library(dplyr)
library(forestplot)
library(readxl)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(viridis)  
library(patchwork)
library(cowplot)
library(mada)
library(meta)
library(gt)
library(grid)
library(officer)
library(flextable)
library(clipr)
library(Cairo)

library(extrafont)
library(showtext)
library(ragg)


# Font --------------------------------------------------------------------

# font_add(family = "MyHelv",regular = "/System/Library/Fonts/Helvetica.ttc")
# showtext_auto()

# for data cleaning --------------------------------------------------------------

mdset <- c("ALSmimics", "SMA", "IPN", "IPMN", "IMPN","ID", 
           "MMN", "HSP", "SBMA", "PLS", "PMA", "CSM", "MD")

dcset <- c("DC", "AD", "NonALS", "iNPH", "MND vs DC", "PD", 
           "MSA", "PSP", "bvFTD", "FTD", "Dementia", "NIMPN", "MS")

hcset <- c("HC", "HC + DC", "HC+DC","DC+HC","NonND", "NonALS + NonND")

mddc <- unique(union(mdset,dcset))
hcdc <- unique(union(hcset,dcset))

bmsort <- c("NfL", "pNfH", "NfH","CHIT1", "YKL40", "CHI3L2","GFAP", "t-tau", "p-tau", "p-tau181",
            "ptau_ttau_ratio","p-tau/t-tau", "TDP43","sEV-TDP43","sEV TDP43","Aβ42",
            "UCHL1","sTREM2","Spp1","Creatinine","CysC","Cystatin C",
            "CK","Ferritin","IL-2","IL-6","IL-10","NLR",
            "Leptin", "Ghrelin", "GIP","GLP-1","C-peptide","insulin","Glucagon")


# basic calculation - Function ---------------------------------------------------

# Function to calculate Sensitivity (SE) for vector inputs
CalcSENS <- function(TP, FN) {
  ifelse((TP + FN) == 0, NA, TP / (TP + FN))
}

# Function to calculate Specificity (SP) for vector inputs
CalcSPEC <- function(TN, FP) {
  ifelse((TN + FP) == 0, NA, TN / (TN + FP))
}


# Make "Diagnostic Value Summary Table" - Function --------------------------------

DiagStbl <- function(data, title = "Biomarker Summary Table",
                     sourceNote = NULL) {
  
  # HSROC (Reitsma) model for Sensitivity & Specificity & sAUC
  hsroc_model <- mada::reitsma(data)
  hsroc_summary <- summary(hsroc_model)
  
  hsroc_sensitivity <- hsroc_summary$coefficients["sensitivity", "Estimate"]
  hsroc_sensitivity_lower <- hsroc_summary$coefficients["sensitivity", "95%ci.lb"]
  hsroc_sensitivity_upper <- hsroc_summary$coefficients["sensitivity", "95%ci.ub"]
  
  hsroc_specificity <- 1 - hsroc_summary$coefficients["false pos. rate", "Estimate"]
  hsroc_specificity_lower <- 1 - hsroc_summary$coefficients["false pos. rate", "95%ci.ub"]
  hsroc_specificity_upper <- 1 - hsroc_summary$coefficients["false pos. rate", "95%ci.lb"]
  
  # Random effects meta-analysis for AUC
  meta_model_auc <- meta::metagen(
    TE = data$calcAUC,
    seTE = data$calcAUC_SE,
    sm = "AUC",
    method.tau = "REML",
    method.random.ci = "HK",
    adhoc.hakn.ci = "se"
  )
  
  RE_AUC <- meta_model_auc$TE.random
  RE_AUC_LowerCI <- pmax(meta_model_auc$lower.random, 0)
  RE_AUC_UpperCI <- pmin(meta_model_auc$upper.random, 1)
  
  # Prepare individual study part
  table_data <- data %>%
    mutate(
      Sensitivity_CI = paste0(
        formatC(Sensitivity, format = "f", digits = 2), " (",
        formatC(pmax(0, Sensitivity - 1.96 * Sensitivity_SE), format = "f", digits = 2), "-",
        formatC(pmin(1, Sensitivity + 1.96 * Sensitivity_SE), format = "f", digits = 2), ")"
      ),
      Specificity_CI = paste0(
        formatC(Specificity, format = "f", digits = 2), " (",
        formatC(pmax(0, Specificity - 1.96 * Specificity_SE), format = "f", digits = 2), "-",
        formatC(pmin(1, Specificity + 1.96 * Specificity_SE), format = "f", digits = 2), ")"
      ),
      calcAUC_CI = paste0(
        formatC(calcAUC, format = "f", digits = 2), " (",
        formatC(pmax(0, calcAUC - 1.96 * calcAUC_SE), format = "f", digits = 2), "-",
        formatC(pmin(1, calcAUC + 1.96 * calcAUC_SE), format = "f", digits = 2), ")"
      )
    ) %>%
    select(Study, ALS_setting_primary, ALS_setting_secondary, Control, n1, n2, Sensitivity_CI, Specificity_CI, calcAUC_CI) %>%
    distinct()
  
  # Create Overall row
  overall_row <- tibble(
    Study = "Overall",
    ALS_setting_primary = NA,
    ALS_setting_secondary = NA,
    Control = NA,
    n1 = sum(data$n1, na.rm = TRUE),
    n2 = sum(data$n2, na.rm = TRUE),
    Sensitivity_CI = paste0(
      formatC(hsroc_sensitivity, format = "f", digits = 2), " (",
      formatC(pmax(0, hsroc_sensitivity_lower), format = "f", digits = 2), "-",
      formatC(pmin(1, hsroc_sensitivity_upper), format = "f", digits = 2), ")"
    ),
    Specificity_CI = paste0(
      formatC(hsroc_specificity, format = "f", digits = 2), " (",
      formatC(pmax(0, hsroc_specificity_lower), format = "f", digits = 2), "-",
      formatC(pmin(1, hsroc_specificity_upper), format = "f", digits = 2), ")"
    ),
    calcAUC_CI = paste0(
      formatC(RE_AUC, format = "f", digits = 2), " (",
      formatC(RE_AUC_LowerCI, format = "f", digits = 2), "-",
      formatC(RE_AUC_UpperCI, format = "f", digits = 2), ")"
    )
  )
  
  # Combine
  final_table <- bind_rows(table_data, overall_row) %>%
    arrange(factor(Study, levels = c(sort(unique(data$Study)), "Overall")))
  
  # Create gt table
  final_table %>%
    gt() %>%
    tab_header(
      title = title
    ) %>%
    cols_label(
      Study = "Study",
      ALS_setting_primary = "ALS diagnostic criteria",
      ALS_setting_secondary = "class",
      Control = "Control",
      n1 = "ALS (n)",
      n2 = "Controls (n)",
      Sensitivity_CI = "Sensitivity (95% CI)",
      Specificity_CI = "Specificity (95% CI)",
      calcAUC_CI = "AUC (95% CI)"
    ) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(
        rows = Study == "Overall"
      )
    ) %>%
    tab_footnote(
      footnote = "overall AUC is pooled using random-effects model (REML, Knapp-Hartung adjustment).",
      locations = cells_column_labels(columns = calcAUC_CI)
    ) %>%
    tab_options(
      table.font.size = "small"
    ) %>% 
    tab_source_note(
      source_note = sourceNote
    )
  
}


# plot HSROC curve - Function -------------------------------------------------

fHSROC <- function(data, plot_title = NULL) {
  # Save original graphical parameters
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))  
  # Set square plotting area and increase top margin for title + annotation
  par(pty = "s",mar = c(5.1, 4.1, 6, 2.1))
  
  # If no title is provided, generate one from FluidType, TargetBiomarker, and ControlCateg
  if (is.null(plot_title)) {
    fluid <- paste(unique(data$FluidType), collapse = ", ")
    biomarker <- paste(unique(data$TargetBiomarker), collapse = ", ")
    control <- paste(unique(data$ControlCateg), collapse = ", ")
    plot_title <- paste0("Summary ROC Curve (", fluid, " - ", biomarker, ", vs ", control, ")")
  }
  
  
  # Fit the bivariate model of Reitsma et al.(2005) - equivalent to the HSROC
  hsroc_result <- mada::reitsma(data)
  # Extract summary statistics (Se)
  SROCse <- summary(hsroc_result)$coefficients["sensitivity", "Estimate"]
  SROCse_ci_lower <- summary(hsroc_result)$coefficients["sensitivity", "95%ci.lb"]
  SROCse_ci_upper <- summary(hsroc_result)$coefficients["sensitivity", "95%ci.ub"]
  # Extract summary statistics (Sp)
  SROCsp <- 1 - summary(hsroc_result)$coefficients["false pos. rate", "Estimate"]
  SROCsp_ci_lower <- 1 - summary(hsroc_result)$coefficients["false pos. rate", "95%ci.ub"]
  SROCsp_ci_upper <- 1 - summary(hsroc_result)$coefficients["false pos. rate", "95%ci.lb"]
  # Extract the Summary AUC in SROC model
  SROCAUC <- AUC(hsroc_result)$AUC
  # Count sample sizes
  totalN <- sum(data$n3)
  totalALS <- sum(data$n1)
  totalCon <- sum(data$n2)
  # Plot the SROC curve
  plot(hsroc_result, lwd=1, asp=1, main = plot_title)
  lines(sroc(hsroc_result), lwd=2)
  # Point the original data in the same figure
  points(1 - SROCsp, SROCse, col = "blue", pch = 19, cex = 1.5)
  points(1 - data$Specificity, data$Sensitivity, col = "black", pch = 19, cex = 0.75)
  # Text the Summary AUC, sensitivity, and sepcificity
  mtext(paste0("Total N: ", totalN, " (ALS: ", totalALS, ", Con: ", totalCon, ")"), 
        side = 3, line = 0.5, adj = 0, cex = 0.9)  # Top margin (side=3), left(adj=0)
  text(0.55, 0.15, labels = paste0("HSROC AUC: ", round(SROCAUC, 3)), adj = 0)
  text(0.55, 0.10, 
       labels = paste0("Sensitivity: ", round(SROCse, 3), 
                       " (", round(SROCse_ci_lower, 3), ", ", round(SROCse_ci_upper, 3), ")"), 
       col = "blue", adj = 0)
  text(0.55, 0.05, 
       labels = paste0("Specificity: ", round(SROCsp, 3), 
                       " (", round(SROCsp_ci_lower, 3), ", ", round(SROCsp_ci_upper, 3), ")"), 
       col = "blue", adj = 0)
  #add legend
  legend("bottomleft",c("SROC", "conf.region"),lwd=c(2,1),bg="white")
  
  # Results
  return(list(
    Sensitivity = c(Estimate = SROCse, CI_Lower = SROCse_ci_lower, CI_Upper = SROCse_ci_upper),
    Specificity = c(Estimate = SROCsp, CI_Lower = SROCsp_ci_lower, CI_Upper = SROCsp_ci_upper),
    AUC = SROCAUC
  ))
}


## Function 

fHSROC_compare <- function(data1, data2, plot_title=NULL,
                           label1 = "CSF NfL", label2 = "Blood NfL",
                           col1   = "#F1C40F", col2   = "#FF69B4",
                           pch1   = 19,       pch2   = 17,
                           cex_sum = 1.8,     cex_pts = 0.75) {
  
  # Save original graphical parameters
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))  
  par(pty = "s")
  
  # Fit HSROC models
  hsroc1 <- reitsma(data1)
  hsroc2 <- reitsma(data2)
  
  # Extract summaries for dataset 1
  sum1 <- summary(hsroc1)
  se1   <- sum1$coefficients["sensitivity","Estimate"]
  se1_L <- sum1$coefficients["sensitivity","95%ci.lb"]
  se1_U <- sum1$coefficients["sensitivity","95%ci.ub"]
  sp1   <- 1 - sum1$coefficients["false pos. rate","Estimate"]
  sp1_L <- 1 - sum1$coefficients["false pos. rate","95%ci.ub"]
  sp1_U <- 1 - sum1$coefficients["false pos. rate","95%ci.lb"]
  auc1  <- AUC(hsroc1)$AUC
  
  # Extract summaries for dataset 2
  sum2 <- summary(hsroc2)
  se2   <- sum2$coefficients["sensitivity","Estimate"]
  se2_L <- sum2$coefficients["sensitivity","95%ci.lb"]
  se2_U <- sum2$coefficients["sensitivity","95%ci.ub"]
  sp2   <- 1 - sum2$coefficients["false pos. rate","Estimate"]
  sp2_L <- 1 - sum2$coefficients["false pos. rate","95%ci.ub"]
  sp2_U <- 1 - sum2$coefficients["false pos. rate","95%ci.lb"]
  auc2  <- AUC(hsroc2)$AUC
  
  # Totals
  tot1 <- sum(data1$n3); a1 <- sum(data1$n1); c1 <- sum(data1$n2)
  tot2 <- sum(data2$n3); a2 <- sum(data2$n1); c2 <- sum(data2$n2)
  
  # 1) Empty ROC plot region 
  plot(NA, NA,
       xlim = c(0,1), ylim = c(0,1),
       xlab = "1 - Specificity", ylab = "Sensitivity",
       asp  = 1, main = plot_title)
  
  # 2) Add HSROC curves 
  lines(sroc(hsroc1), col = col1, lwd = 1)
  lines(sroc(hsroc2), col = col2, lwd = 1, lty = 2)
  
  # 3) Summary points 
  points(1 - sp1, se1, col = col1, pch = 15, cex = cex_sum)
  points(1 - sp2, se2, col = col2, pch = 15, cex = cex_sum)
  
  # 4) Individual study points 
  points(1 - data1$Specificity, data1$Sensitivity, col = col1, pch = pch1, cex = cex_pts)
  points(1 - data2$Specificity, data2$Sensitivity, col = col2, pch = pch2, cex = cex_pts)
  
  # 5) Sample size annotations: two lines per dataset 
  
  # # CSF line
  # mtext(
  #   sprintf("%s N=%d (ALS=%d, Con=%d)", label1, tot1, a1, c1),
  #   side = 3, adj = 1, line = 2, col = col1, font = 1)
  
  # # Blood line
  # mtext(
  #   sprintf("%s N=%d (ALS=%d, Con=%d)", label2, tot2, a2, c2),
  #   side = 3, adj = 1, line = 1, col = col2, font = 1)
  
  # 6) AUC & CI text with blank line between
  
  ## Label 1 
  #text(0.55, 0.24, paste0(label1, " AUC: ", round(auc1,3)), col = col1, adj = 0)
  text(0.55, 0.24, sprintf("CSF, ALS %d vs Con %d", a1, c1), col = col1, adj = 0)
  text(0.55, 0.20, sprintf("- Se: %.3f (%.3f-%.3f)", se1, se1_L, se1_U), col = col1, adj = 0)
  text(0.55, 0.16, sprintf("- Sp: %.3f (%.3f-%.3f)", sp1, sp1_L, sp1_U), col = col1, adj = 0)
  
  # Label 2
  #text(0.55, 0.09, paste0(label2, " AUC: ", round(auc2,3)), col = col2, adj = 0)
  text(0.55, 0.09, sprintf("Blood, ALS %d vs Con %d", a2, c2), col = col2, adj = 0)
  text(0.55, 0.05, sprintf("- Se: %.3f (%.3f-%.3f)", se2, se2_L, se2_U), col = col2, adj = 0)
  text(0.55, 0.01, sprintf("- Sp: %.3f (%.3f-%.3f)", sp2, sp2_L, sp2_U), col = col2, adj = 0)
  
  # 7) Legend outside plot 
  
  legend("bottomleft", inset = c(0,0), xpd = TRUE,
         legend = c(label1, label2),
         col    = c(col1, col2),
         pch    = c(15, 15),
         lty    = c(1,2),
         bg     = "white")
  
  invisible(list(
    HSROC1 = list(sensitivity = c(Estimate=se1, CI_Lower=se1_L, CI_Upper=se1_U),
                  specificity = c(Estimate=sp1, CI_Lower=sp1_L, CI_Upper=sp1_U), AUC = auc1),
    HSROC2 = list(sensitivity = c(Estimate=se2, CI_Lower=se2_L, CI_Upper=se2_U),
                  specificity = c(Estimate=sp2, CI_Lower=sp2_L, CI_Upper=sp2_U), AUC = auc2)
  ))
}


# pool AUC - Function -----------------------------------------------------

metaAUC <- function(data,
                    xlab = NULL, 
                    pattern = 1,
                    subname = NULL,
                    filter_FluidType = NULL,
                    filter_TargetBiomarker = NULL) {
  
  # Calculate z-value for 95% CI
  z975 <- qnorm(0.975)
  
  # Filter and arrange the data
  data2 <- data %>%
    filter(
      !is.na(calcAUC),
      !is.na(calcAUC_SE),
      !is.na(Study),
      !is.na(TargetBiomarker)
    ) %>%
    { if (!is.null(filter_FluidType))        filter(., FluidType == filter_FluidType) else . } %>%
    { if (!is.null(filter_TargetBiomarker))  filter(., TargetBiomarker == filter_TargetBiomarker) else . } %>%
    arrange(Study) %>%
    mutate(
      TE_clipped = pmin(pmax(calcAUC, 0), 1),
      seTE_adj   = pmin(
        calcAUC_SE,
        (1 - TE_clipped) / z975,
        TE_clipped      / z975
      ),
      Cut_off_Value = ifelse(is.na(Cut_off_Value), "not specified", Cut_off_Value)
    )
  
  # Dynamically create xlab if not provided
  if (is.null(xlab)) {
    fluid <- paste(unique(data2$FluidType), collapse = ", ")
    biomarker <- paste(unique(data2$TargetBiomarker), collapse = ", ")
    control <- paste(unique(data2$ControlCateg), collapse = ", ")
    xlab <- paste0("AUC (95% CI), ", fluid, " - ", biomarker, " vs ", control)
  }
  
  # Run random-effects meta-analysis (classic CI only)
  meta_auc <- metagen(
    TE = TE_clipped,
    seTE = seTE_adj,
    studlab = Study,
    data = data2,
    sm = "AUC",
    method.tau = "REML",
    method.random.ci = "classic",  # <- unified to classic
    adhoc.hakn.ci = "se",
    subgroup = TargetBiomarker
  )
  
  # Clip pooled AUC CI to [0,1]
  meta_auc$lower.random    <- pmax(meta_auc$lower.random, 0)
  meta_auc$upper.random    <- pmin(meta_auc$upper.random, 1)
  meta_auc$lower.random.ci <- pmax(meta_auc$lower.random.ci, 0)
  meta_auc$upper.random.ci <- pmin(meta_auc$upper.random.ci, 1)
  
  #meta_auc$TE.random <- pmin(pmax(meta_auc$TE.random, 0), 1)
  
  
  # Set up left columns
  if (pattern == 1) {
    leftcols <- c("studlab", "AssayMethod", "n1", "n2")
    leftlabs <- c("Study", "Assay", "ALS (n)", "Controls (n)")
    just.addcols <- c("left", "right", "right")
  } else if (pattern == 2) {
    leftcols <- c("studlab", "AssayMethod", "n1")
    leftlabs <- c("Study", "Assay", "ALS (n)")
    just.addcols <- c("left", "right")
  } else if (pattern == 3) {
    leftcols <- c("studlab")
    leftlabs <- c("Study")
    just.addcols <- c("left")
  } else {
    stop("pattern must be 1, 2, or 3")
  }
  
  # Generate forest plot
  meta::forest(meta_auc,
               leftcols = leftcols,
               leftlabs = leftlabs,
               just.addcols = just.addcols,
               rightcols = c("effect", "ci", "w.random"),
               
               overall = FALSE,
               overall.hetstat = FALSE,
               common = FALSE,
               subgroup = TRUE,
               sep.subgroup = " - ",
               subgroup.name = subname,
               test.subgroup = FALSE,
               test.effect.subgroup = FALSE,
               
               print.pval.Q = FALSE,
               print.Q = TRUE,
               print.I2 = TRUE,
               print.Q.subgroup = TRUE,
               print.I2.subgroup = TRUE,
               
               type.study = "square",
               col.square = "#6B58A6",
               col.square.lines = "#6B58A6",
               col.diamond = "#6B58A6",
               col.diamond.lines = "#6B58A6",
               col.study = "#6B58A6",
               col.subgroup = "#333322",
               digits = 2,
               xlab = xlab,
               xlim = c(0, 1),
               ref = 0.5,
               xlab.pos = 0,
               fs.xlab = 14,
               ff.xlab = "bold.italic",
               colgap = "7mm",
               new = TRUE,
               header.line = "below",
  )
  
  invisible(meta_auc)
}


analyze_dataset <- function(data, type, name, comparison) {
  library(meta)
  
  # Data Prep
  z975 <- qnorm(0.975)
  data_auc <- data %>%
    filter(
      !is.na(calcAUC),
      !is.na(calcAUC_SE),
      !is.na(Study),
      !is.na(TargetBiomarker)
    ) %>%
    arrange(Study) %>%
    mutate(
      TE_clipped = pmin(pmax(calcAUC, 0), 1),
      seTE_adj   = pmin(
        calcAUC_SE,
        (1 - TE_clipped) / z975,
        TE_clipped      / z975
      ),
      Cut_off_Value = ifelse(is.na(Cut_off_Value), "not specified", Cut_off_Value)
    )
  
  meta_model <- meta::metagen(
    TE = data_auc$calcAUC,
    seTE = data_auc$calcAUC_SE,
    studlab = data_auc$Study,
    data = data_auc,
    sm = "AUC",
    method.tau = "REML",
    method.random.ci = "classic",
    adhoc.hakn.ci = "se"
  )
  
  phm_model <- mada::phm(data_auc)
  
  # HSROC model using mada::reitsma
  hsroc_model <- mada::reitsma(data)
  hsroc_summary <- summary(hsroc_model)
  
  # Extract AUC, Sensitivity, and Specificity with 95% CI
  hsroc_auc <- hsroc_summary$AUC$AUC
  hsroc_sensitivity <- hsroc_summary$coefficients["sensitivity", "Estimate"]
  hsroc_sensitivity_lower <- hsroc_summary$coefficients["sensitivity", "95%ci.lb"]
  hsroc_sensitivity_upper <- hsroc_summary$coefficients["sensitivity", "95%ci.ub"]
  hsroc_specificity <- 1 - hsroc_summary$coefficients["false pos. rate", "Estimate"]
  hsroc_specificity_lower <- 1 - hsroc_summary$coefficients["false pos. rate", "95%ci.ub"]
  hsroc_specificity_upper <- 1 - hsroc_summary$coefficients["false pos. rate", "95%ci.lb"]
  
  # Counting the kinds of Controls
  control_counts <- table(data$Control)
  
  # Dataframe results
  result <- data.frame(
    Type = type,
    Biomarker = name,
    Comparison = comparison,
    n1 = sum(data$n1),
    n2 = sum(data$n2),
    n1auc = sum(data_auc$n1),
    n2auc = sum(data_auc$n2),
    Cohorts = nrow(data),
    Cohorts_AUC = meta_model$k,
    Controls = paste(names(control_counts), control_counts, sep = "=", collapse = "; "),
    
    PHM_AUC = AUC(phm_model)$AUC,
    PHM_AUC_LowerCI = pmax(AUC(phm_model)$ci[2], 0),
    PHM_AUC_UpperCI = pmin(AUC(phm_model)$ci[1], 1),
    
    RE_AUC = meta_model$TE.random,
    RE_AUC_LowerCI = pmax(meta_model$lower.random, 0),
    RE_AUC_UpperCI = pmin(meta_model$upper.random, 1),
    
    HSROC_AUC = hsroc_auc,
    HSROC_Sensitivity = hsroc_sensitivity,
    HSROC_Sensitivity_LowerCI = hsroc_sensitivity_lower,
    HSROC_Sensitivity_UpperCI = hsroc_sensitivity_upper,
    HSROC_Specificity = hsroc_specificity,
    HSROC_Specificity_LowerCI = hsroc_specificity_lower,
    HSROC_Specificity_UpperCI = hsroc_specificity_upper
  )
  
  return(result)
}



### format results to make summary table
format_results <- function(results) {
  results %>%
    mutate(
      "N, ALS" = n1,"N, Con" = n2,
      "N, ALS_auc" = n1auc,"N, Con_auc" = n2auc,
      `AUC (RE model)` = sprintf("%.2f (%.2f, %.2f)", RE_AUC, RE_AUC_LowerCI, RE_AUC_UpperCI),
      `AUC (phm model)` = sprintf("%.2f (%.2f, %.2f)", PHM_AUC, PHM_AUC_LowerCI, PHM_AUC_UpperCI),
      `AUC (HSROC model)` = sprintf("%.2f", HSROC_AUC),
      `Pooled Se.` = sprintf("%.2f (%.2f, %.2f)", HSROC_Sensitivity, HSROC_Sensitivity_LowerCI, HSROC_Sensitivity_UpperCI),
      `Pooled Sp.` = sprintf("%.2f (%.2f, %.2f)", HSROC_Specificity, HSROC_Specificity_LowerCI, HSROC_Specificity_UpperCI)
    ) %>%
    select(
      Type,
      Biomarker,
      Comparison,
      Cohorts,
      Cohorts_AUC,
      "N, ALS","N, Con",
      "N, ALS_auc","N, Con_auc",
      Controls,
      `Pooled Se.`,
      `Pooled Sp.`,
      `AUC (HSROC model)`,
      `AUC (phm model)`,
      `AUC (RE model)`,
    ) %>% 
    arrange(desc(Cohorts),
            desc(`AUC (RE model)`),
            desc(`AUC (phm model)`)
    )
}


## gt () SROC func --------------------------------------------------------------

gtt <- function(data, groupname_col = "Type", row_group_as_column = FALSE) {
  data %>% 
    gt(groupname_col = groupname_col, row_group_as_column = row_group_as_column) %>%  # Add row_group_as_column
    # Remove unnecessary horizontal lines
    tab_style(
      style = cell_borders(
        sides = c("top", "bottom"),
        color = "transparent", # Make unwanted borders transparent
        weight = px(0)         # Set border width to 0
      ),
      locations = cells_body()
    ) %>%
    # Make the title row bold
    tab_style(
      style = list(cell_text(weight = "bold")),
      locations = cells_column_labels()
    ) %>% 
    cols_label(
      `Pooled Se.` = "Summary Se.*",
      `Pooled Sp.` = "Sumamry Sp.*",
      `AUC (RE model)` = "AUC (RE model)‡",
      `AUC (HSROC model)` = "sAUC†"
    ) %>%
    # Add a note below the table
    tab_source_note(
      source_note = c(
        "*) Summary Points on HSROC Curve, 95% Confidence Interval",
        "†) Summary AUC of HSROC Curve",
        "‡) Random Effect model, 95 % Confidence Interval"
      )
    )
}

gtt2 <- function(data, groupname_col = "Type", row_group_as_column = FALSE) {
  data %>% 
    gt(groupname_col = groupname_col, row_group_as_column = row_group_as_column) %>% 
    # Remove unnecessary horizontal lines
    tab_style(
      style = cell_borders(
        sides = c("top", "bottom"),
        color = "transparent", # Make unwanted borders transparent
        weight = px(0)         # Set border width to 0
      ),
      locations = cells_body()
    ) %>%
    # Make the title row bold
    tab_style(
      style = list(cell_text(weight = "bold")),
      locations = cells_column_labels()
    ) %>% 
    cols_label(
      `Pooled Se.` = "Summary Se.*",
      `Pooled Sp.` = "Summary Sp.*",
      `AUC (RE model)` = "AUC (RE model)‡",
      `AUC (HSROC model)` = "sAUC†"
    ) %>%
    # Add underlines for each Comparison group
    tab_style(
      style = cell_borders(
        sides = "bottom",
      ),
      locations = cells_body(
        rows = Comparison != lag(Comparison) | is.na(lag(Comparison)) # First row of each Comparison group
      )
    ) %>%
    # Add a note below the table
    tab_source_note(
      source_note = c(
        "*) Summary Points on HSROC Curve, 95% Confidence Interval",
        "†) Summary AUC of HSROC Curve"
        #"‡) Random Effect model, 95 % Confidence Interval"
      )
    )
}

gtROC <- function(data){
  data %>% gt() %>% 
    tab_row_group(label = "CSF Biomarkers", rows = Type == "CSF") %>%
    tab_row_group(label = "Blood Biomarkers", rows = Type == "Blood") %>%
    cols_label(Type = "Fluid Type") %>%
    tab_options(table.font.size = px(14),row_group.font.weight = "bold") %>% 
    tab_header(md("**Pooled Diagnostic values**")) %>% 
    cols_align(align = "left", columns = Biomarker) %>% 
    cols_hide("Type")
}




# pool SMD - Function ----------------------------------------------------------

SMDplot <- function(df, 
                    pattern = 2, 
                    Fluid = "blood", 
                    Con = NULL,  # If NULL, do not filter Control
                    Target = NULL, 
                    xlab = NULL,
                    filterN = 2,
                    subname = NULL) {
  # Validate pattern input
  if (!pattern %in% c(1, 2)) {
    stop("pattern must be either 1 or 2.")
  }
  
  # Set default xlab if not specified
  if (is.null(xlab)) {
    xlab <- paste0(Fluid, " biomarkers SMD (ALS vs Con), 95%CI")
  }
  
  # pattern 1 requires Target to be specified
  if (pattern == 1 && is.null(Target)) {
    stop("For pattern 1, Target must be specified.")
  }
  
  # --- Data filtering ---
  # Select necessary rows and filter by FluidType
  sm <- df %>%
    filter(FluidType == Fluid) %>% 
    mutate(TargetBiomarker = factor(TargetBiomarker,
                                    levels = c(bmsort, 
                                               setdiff(unique(TargetBiomarker), bmsort)),
                                    ordered = TRUE))
  
  if (pattern == 1) {
    # Pattern 1: Filter only the specified TargetBiomarker
    sm <- sm %>%
      filter(TargetBiomarker == Target) %>%
      group_by(Control) %>% filter(n() >= filterN) %>%ungroup()
    
  } else if (pattern == 2) {
    # Pattern 2: Filter by Control only if Con is specified
    if (!is.null(Con)) {
      sm <- sm %>% filter(con_sort == Con)
    }
    sm <- sm %>%
      group_by(FluidType, TargetBiomarker) %>% 
      filter(n() >= filterN) %>%  ungroup() %>% 
      arrange(Study) #Arrange alphabetically
  }
  
  # --- Determine subgroup variable ---
  # Pattern 1 uses con_sort for subgroup analysis, Pattern 2 uses TargetBiomarker
  subgroup_var <- if (pattern == 1) "con_sort" else "TargetBiomarker"
  
  # Extract the appropriate subgroup column
  subgroup_data <- sm[[subgroup_var]]
  
  # --- Perform meta-analysis ---
  meta_result <- metacont(
    data = sm,
    n.e = SampleSize_ALS,
    n.c = SampleSize_Con,
    mean.e = Mean_ALS, sd.e = SD_ALS,
    mean.c = Mean_Con, sd.c = SD_Con,
    median.e = Med_ALS, q1.e = Q1_ALS, q3.e = Q3_ALS,
    median.c = Med_Con, q1.c = Q1_Con, q3.c = Q3_Con,
    sm = "SMD",
    studlab = Study,
    subgroup = subgroup_data
  )
  
  # --- Generate forest plot ---
  
  meta::forest(meta_result,
               xlim = c(-3.5, 3.5),
               rightcols = c("effect", "ci", "w.random"),
               rightlabs = c("SMD", "95%CI", "Weight"),
               subgroup = TRUE, overall = FALSE, overall.hetstat = FALSE,
               subgroup.name = subname, sep.subgroup = " - ",
               pooled.totals = TRUE,
               test.effect.subgroup = FALSE,
               test.subgroup = FALSE,
               ff.test.subgroup = "bold", ff.random = "bold",
               type.study = "square", print.Q.subgroup = TRUE,
               col.square = "#6B58A6", col.square.lines = "#6B58A6",
               col.diamond = "#6B58A6", col.diamond.lines = "#6B58A6",
               col.study = "#6B58A6",   
               col.subgroup = "#333322",
               digits = 2, digits.mean = 2, digits.sd = 2,
               leftcols = c("studlab", "Control", "AssayMethod", "n.e", "n.c"),
               just.addcols.left = c("right", "left"),
               leftlabs = c("Study", "Controls def.", "Assay", "N, ALS", "N, Con"),
               xlab = xlab,               
               xlab.pos = 0, fs.xlab = 14, 
               ff.xlab = "bold.italic",  # X-axis label adjustments
               header.line = "below",
               colgap = "7mm",
               common = FALSE,
               new = TRUE
  )
  return(meta_result)
}


### To make summary table
metaSMDsum <- function(data, 
                       biomarker = NULL, 
                       fluid_type = "CSF", 
                       control = NULL,  # Allow NULL to disable control filtering
                       filterN = 2, 
                       sort = bmsort) {
  # Load required packages
  require(dplyr)
  require(tibble)
  require(meta)
  
  # Define base filtering (excluding con_sort filter if control is NULL)
  filtered_data_base <- data %>%
    filter(FluidType == fluid_type) %>%
    { if (!is.null(control)) filter(., con_sort == control) else . } %>% 
    group_by(FluidType, TargetBiomarker) %>%
    filter(n() >= filterN) %>%
    ungroup()
  
  # If no biomarker vector is provided, extract unique TargetBiomarker values
  if (is.null(biomarker)) {
    biomarker <- filtered_data_base %>%
      pull(TargetBiomarker) %>%
      unique()
  }
  
  # If a sort order is provided, reorder the biomarker vector accordingly
  if (!is.null(sort)) {
    biomarker <- c(intersect(sort, biomarker), setdiff(biomarker, sort))
  }
  
  # For each biomarker, perform meta-analysis using meta::metacont and create a summary row
  results_list <- lapply(biomarker, function(bio) {
    # Filter data for the specified fluid_type, control, and biomarker
    filtered_data <- filtered_data_base %>%
      filter(TargetBiomarker == bio)
    
    # If no valid data remains, return NULL
    if (nrow(filtered_data) == 0) return(NULL)
    
    # Perform meta-analysis using metacont
    meta_obj <- meta::metacont(
      n.e = filtered_data$SampleSize_ALS,
      n.c = filtered_data$SampleSize_Con,
      mean.e = filtered_data$Mean_ALS,
      sd.e = filtered_data$SD_ALS,
      mean.c = filtered_data$Mean_Con,
      sd.c = filtered_data$SD_Con,
      median.e = filtered_data$Med_ALS, 
      q1.e = filtered_data$Q1_ALS, 
      q3.e = filtered_data$Q3_ALS,
      median.c = filtered_data$Med_Con, 
      q1.c = filtered_data$Q1_Con, 
      q3.c = filtered_data$Q3_Con,
      sm = "SMD",
      studlab = filtered_data$Study
    )
    
    ## Egger et al (1997) Linear regression test of funnel plot asymmetry
    ## Reporting bias
    egger_str <- tryCatch({
      egger <- meta::metabias(meta_obj, k.min=4)  # default k.min=10
      if (is.null(egger$estimate) || 
          is.null(egger$estimate["bias"]) || 
          is.null(egger$estimate["se.bias"])) {
        NA
      } else {
        egger_intercept <- egger$estimate["bias"]
        egger_se <- egger$estimate["se.bias"]
        egger_lci <- egger_intercept - 1.96 * egger_se
        egger_uci <- egger_intercept + 1.96 * egger_se
        sprintf("%.2f [%.2f; %.2f], p=%.2f", 
                egger_intercept, egger_lci, egger_uci, egger$pval)
      }
    }, error = function(e) { 
      NA 
    })
    
    # Create the formatted SMD result string
    smd_str <- sprintf("%.2f [%.2f; %.2f]",
                       meta_obj$TE.random,
                       meta_obj$lower.random,
                       meta_obj$upper.random)
    
    # Create a summary tibble row
    tibble(
      fluid = fluid_type,
      biomarker = bio,
      combined = paste(fluid_type, "-", bio),
      `Number of cohorts` = meta_obj$k.all,
      `N(ALS)` = sum(meta_obj$n.e, na.rm = TRUE),
      `N(con)` = sum(meta_obj$n.c, na.rm = TRUE),
      `SMD [95%CI]` = smd_str,
      `Egger's test` = egger_str,
    )
  })
  
  # Combine the list of tibbles into one tibble and return
  bind_rows(results_list)
}


## gt() SMD

gtSMD <- function(data){
  data %>% gt() %>% 
    tab_row_group(label = "CSF Biomarkers", rows = fluid == "CSF") %>%
    tab_row_group(label = "Blood Biomarkers", rows = fluid == "blood") %>%
    cols_label(fluid = "Fluid Type",biomarker = "Biomarker") %>%
    tab_options(table.font.size = px(14),row_group.font.weight = "bold") %>% 
    tab_header(md("**Biomarker level's SMD, vs HC†, ALS mimics‡, Disease controls**"))
}


# pool Cox HR - Function ------------------------------------------------------------


metaHR <- function(data, 
                   xlab = "Hazard ratio (95%CI)", 
                   pattern = 1,
                   subname = NULL) {
  # Set left column specifications based on the chosen pattern
  if (pattern == 1) {
    leftcols    <- c("studlab", "AssayMethod", "SampleSize_ALS", "def")
    leftlabs    <- c("Study", "Assay", "N", "Comparison")
    just.addcols <- c("left", "right", "left")
  } else if (pattern == 2) {
    leftcols    <- c("studlab", "AssayMethod", "SampleSize_ALS", "HR_scaled")
    leftlabs    <- c("Study", "Assay", "N", "scale")
    just.addcols <- c("left", "right", "right")
  } else if (pattern == 3){
    leftcols    <- c("studlab", "SampleSize_ALS")
    leftlabs    <- c("Study",  "N")
    just.addcols <- c("right")      
  } else {
    stop("The 'pattern' argument must be either 1-3.")
  }  
  # Sort by Study alphabetically
  data <- data %>% arrange(Study)
  # Perform meta-analysis using the metagen function from the meta package.
  meta_HR <- meta::metagen(
    data         = data,
    TE           = logHR,             # column with log-transformed HR
    seTE         = selogHR,           # column with standard error of logHR
    studlab      = Study,           # study label column
    sm           = "HR",              # summary measure: hazard ratio
    method.tau   = "PM",              # between-study variance estimator (Paule and Mandel, 1982)
    # other options: REML, DL, ML, HS, SJ, HE, EB. see help(meta) and search "method.tau"
    method.random.ci = "classic",      # use classic CI calculation (no HK adjustment)
    subgroup     = TargetBiomarker    # subgroup variable
  )
  # Create the forest plot using the meta::forest function with the specified parameters.
  meta::forest(meta_HR,
               leftcols    = leftcols, 
               leftlabs    = leftlabs,
               just.addcols = just.addcols,
               rightcols = c("effect", "ci", "w.random"),
               subgroup    = TRUE, sep.subgroup = " - ", subgroup.name = subname, 
               overall     = FALSE, overall.hetstat = FALSE,
               pooled.totals = TRUE,
               common      = FALSE,
               test.effect.subgroup = FALSE,
               test.subgroup = FALSE,
               ff.test.subgroup = "bold", 
               ff.random   = "bold",
               type.study  = "square", 
               print.Q.subgroup = TRUE,
               col.square  = "#6B58A6", 
               col.square.lines = "#6B58A6",
               col.diamond = "#6B58A6", 
               col.diamond.lines = "#6B58A6",
               col.study   = "#6B58A6",   
               col.subgroup = "#333322",
               digits      = 2,
               xlab        = xlab,
               xlab.pos    = 0, 
               fs.xlab     = 14, 
               ff.xlab     = "bold.italic",
               colgap      = "7mm",
               new         = TRUE,
               #print.Q     = TRUE,
               print.pval.Q =TRUE,
               header.line = "below"
  )
  # Return the meta-analysis object invisibly in case further inspection is needed.
  invisible(meta_HR)
}

## DataCleaning Function
filHR <- function(data, analysis, filterN=2) {
  # analysis_type: one of "ContUni", "ContMulti", "CategUni", or "CategMulti"
  # Set column names based on the analysis type
  if (analysis == "ContUni") {
    outcome <- "HR_death_cont"
    lowerCI <- "HR_death_cont_lCI"
    upperCI <- "HR_death_cont_uCI"
    include_def <- FALSE
  } else if (analysis == "ContMulti") {
    outcome <- "aHR_death_cont"
    lowerCI <- "aHR_death_cont_lCI"
    upperCI <- "aHR_death_cont_uCI"
    include_def <- FALSE
  } else if (analysis == "CategUni") {
    outcome <- "HR_death_categ"
    lowerCI <- "HR_death_categ_lCI"
    upperCI <- "HR_death_categ_uCI"
    def_col <- "HR_death_categ_def"
    include_def <- TRUE
  } else if (analysis == "CategMulti") {
    outcome <- "aHR_death_categ"
    lowerCI <- "aHR_death_categ_lCI"
    upperCI <- "aHR_death_categ_uCI"
    def_col <- "aHR_death_categ_def"
    include_def <- TRUE
  } else {
    stop("Invalid analysis type. Use one of 'ContUni', 'ContMulti', 'CategUni', or 'CategMulti'.")
  }
  
  # Define the required columns that must not be NA
  required_cols <- c("FluidType", "TargetBiomarker", outcome, lowerCI, upperCI)
  
  # Filter data: remove rows with missing values in "required_columns" above,
  # group by FluidType and TargetBiomarker (keeping groups with at least 2 studies),
  # and calculate selogHR and logHR.
  cleaned_data <- data %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(required_cols), ~ !is.na(.))) %>%
    dplyr::group_by(FluidType, TargetBiomarker) %>%
    dplyr::filter(dplyr::n() >= filterN) %>%
    #dplyr::ungroup() %>%
    dplyr::mutate(
      selogHR = (log(.data[[upperCI]]) - log(.data[[lowerCI]])) / (2 * 1.96),
      logHR = log(.data[[outcome]])
    )
  
  # For categorical analyses, add the 'def' column to the data
  if (include_def) {
    cleaned_data <- cleaned_data %>%
      dplyr::mutate(def = .data[[def_col]])
  }
  
  return(cleaned_data)
}


## summary table function
metaHRsum <- function(data, 
                      biomarker = NULL, 
                      fluid_type = "blood", 
                      sort = bmsort) {
  # If no biomarker vector is provided, extract unique TargetBiomarker values for the specified fluid_type
  if (is.null(biomarker)) {
    biomarker <- data %>%
      filter(FluidType == fluid_type) %>%
      pull(TargetBiomarker) %>%
      unique()
  }
  
  # If a sort order is provided, reorder the biomarker vector accordingly
  if (!is.null(sort)) {
    biomarker <- c(intersect(sort, biomarker), setdiff(biomarker, sort))
  }
  
  # For each biomarker, perform meta-analysis using meta::metagen and create a summary row
  results_list <- lapply(biomarker, function(bio) {
    # Filter data for the specified fluid_type and biomarker
    filtered_data <- data %>%
      filter(FluidType == fluid_type, TargetBiomarker == bio)
    
    # If there is no data for this biomarker, return NULL
    if (nrow(filtered_data) == 0) return(NULL)
    
    # Run meta-analysis using meta::metagen (using logHR and selogHR columns)
    meta_obj <- meta::metagen(
      TE = filtered_data$logHR,
      seTE = filtered_data$selogHR,
      studlab = filtered_data$Study,
      sm = "HR",
      method.tau = "PM",
      method.random.ci = "classic"
    )
    
    # Back-transform the log-scale estimates to obtain HR values
    pooled_HR_str <- sprintf("%.2f [%.2f; %.2f]",
                             exp(meta_obj$TE.random),
                             exp(meta_obj$lower.random),
                             exp(meta_obj$upper.random))
    
    ## Egger et al (1997) Linear regression test of funnel plot asymmetry
    ## Reporting bias
    egger_str <- tryCatch({
      egger <- meta::metabias(meta_obj, k.min=4)  # default k.min=10
      if (is.null(egger$estimate) || 
          is.null(egger$estimate["bias"]) || 
          is.null(egger$estimate["se.bias"])) {
        NA
      } else {
        egger_intercept <- egger$estimate["bias"]
        egger_se <- egger$estimate["se.bias"]
        egger_lci <- egger_intercept - 1.96 * egger_se
        egger_uci <- egger_intercept + 1.96 * egger_se
        sprintf("%.2f [%.2f; %.2f], p=%.2f", 
                egger_intercept, egger_lci, egger_uci, egger$pval)
      }
    }, error = function(e) { 
      NA 
    })
    
    
    tibble(
      fluid = fluid_type,
      biomarker = bio,
      combined = paste(fluid_type, "-", bio),
      `Number of cohorts` = meta_obj$k.all,
      `Total sample size` = sum(filtered_data$SampleSize_ALS, na.rm = TRUE),
      `pooled HR [95%CI]` = pooled_HR_str,
      `Egger's test` = egger_str
    )
  })
  
  # Combine the list of tibbles into one tibble and return
  bind_rows(results_list)
}


## gt() HR

gtHR <- function(data){
  data %>% gt() %>% 
    tab_row_group(label = "CSF Biomarkers", rows = fluid == "CSF") %>%
    tab_row_group(label = "Blood Biomarkers", rows = fluid == "blood") %>%
    cols_label(fluid = "Fluid Type",biomarker = "Biomarker") %>%
    tab_options(table.font.size = px(14),row_group.font.weight = "bold") %>% 
    tab_header(md("**Cox hazard model, biomarker levels categorized (dichotomized)**"))
}



# pool Corr coefficients - Function -------------------------------------------

## Filter Function 

filterCOR <- function(df, 
                      cor_type = "Spearman's rho", 
                      fluid_type = "CSF", 
                      cor = "Corr_DPR",
                      filterN = 2) {
  cor_sym <- rlang::sym(cor) # convert a column name to a symbol
  df %>%
    filter(!is.na(FluidType) & !is.na(TargetBiomarker) & !is.na(.data[[cor]])) %>% 
    filter(CorType == cor_type) %>%  
    filter(FluidType == fluid_type) %>%  
    group_by(FluidType, TargetBiomarker) %>%
    filter(n() >= filterN) %>%
    ungroup()
}

## Function to perform meta-analysis and generate a forest plot
mcforest <- function(
    data, cor = "Corr_DPR", 
    xlab = "CSF biomarkers correlation to DPR, Spearman's rho (95% CI)",
    subname = NULL,
    sm = "ZCOR"  # "ZCOR" = Fisher's Z transformation
) {
  # Sort by Study alphabetically
  data <- data %>% 
    mutate(TargetBiomarker = factor(TargetBiomarker,
                                    levels = c(bmsort, 
                                               setdiff(unique(TargetBiomarker), bmsort)),
                                    ordered = TRUE)) %>% 
    arrange(TargetBiomarker, Study)
  
  # Ensure that the correlation column exists in the dataset
  if (!cor %in% colnames(data)) {
    stop(paste("Error: Column", cor, "not found in the dataset."))
  }
  # Extract study labels
  studlab <- data %>% pull(Study)
  
  # Perform meta-analysis
  meta_result <- meta::metacor(
    data = data,         # Data frame
    cor = data[[cor]],   # Extract column by name instead of symbol
    n = SampleSize_ALS,  # Sample size
    studlab = studlab,   # Study labels
    method.tau = "REML", # Random-effects model
    sm = sm,   
    subgroup = TargetBiomarker # Perform subgroup analysis by biomarker
  )
  # Generate the forest plot
  meta::forest(meta_result,
               xlim = c(-1, 1),
               rightcols = c("effect", "ci", "w.random"),
               subgroup = TRUE, overall = FALSE, overall.hetstat = FALSE,
               subgroup.name = subname, pooled.totals = TRUE, sep.subgroup = " - ",
               test.effect.subgroup = FALSE,
               test.subgroup = FALSE,
               ff.test.subgroup = "bold", ff.random = "bold",
               type.study = "square", print.Q.subgroup = TRUE,
               col.square = "#6B58A6", col.square.lines = "#6B58A6",
               col.diamond = "#6B58A6", col.diamond.lines = "#6B58A6",
               col.study = "#6B58A6",   
               col.subgroup = "#333322",
               digits=2,
               leftcols = c("studlab", "AssayMethod", "n"),  # Add AssayMethod column
               leftlabs = c("Study", "Assay", "N"),  # Column labels
               xlab = xlab, # X-axis label
               xlab.pos = 0, fs.xlab = 14, ff.xlab = "bold.italic",  # X-axis label adjustments
               just.addcols = c("left", "right"),
               header.line = "below",
               colgap = "7mm",
               common = FALSE,    # Report only the random-effects model
               #sort.subgroup = TRUE #Sort subgroup names, alphabetically
               new = TRUE
  )
  return(meta_result)
}

metaCorSUM <- function(data, 
                       biomarker = NULL, 
                       fluid_type = "CSF", 
                       cor = "Corr_DPR", 
                       sm = "COR", 
                       filterN = 2, 
                       cor_type = "Spearman's rho",
                       sort = bmsort) {
  # If no biomarker vector is provided, create one from the filtered data
  if (is.null(biomarker)) {
    biomarker <- data %>%
      filterCOR(cor = cor, filterN = filterN, fluid_type = fluid_type, cor_type = cor_type) %>% 
      pull(TargetBiomarker) %>% 
      unique()
  }
  
  # Reorder the biomarker vector:
  biomarker <- c(intersect(sort, biomarker), setdiff(biomarker, sort))
  
  # For each biomarker in the list, perform meta-analysis and create a summary row
  results_list <- lapply(biomarker, function(bio) {
    # Filter data based on the provided criteria and the current biomarker
    filtered_data <- data %>%
      filterCOR(cor = cor, filterN = filterN, 
                fluid_type = fluid_type, cor_type = cor_type) %>%
      filter(TargetBiomarker == bio)
    
    # If there is no data for this biomarker, return NULL
    if (nrow(filtered_data) == 0) return(NULL)
    
    # Run the meta-analysis using meta::metacor
    meta_obj <- meta::metacor(
      data = filtered_data,
      cor = filtered_data[[cor]],       # Extract the specified correlation column
      n = filtered_data$SampleSize_ALS,   # Sample size
      studlab = filtered_data$Study,      # Study labels
      method.tau = "REML",                # Random-effects model
      sm = sm
    )
    
    # Extract the correlation name dynamically (e.g., "DPR" from "Corr_DPR")
    cor_name <- gsub("Corr_", "", cor) # Remove "Corr_" prefix
    # Define the dynamic column name
    corr_col_name <- paste("pooled Corr to", cor_name, "[95%CI]")
    
    ## Egger et al (1997) Linear regression test of funnel plot asymmetry
    ## Reporting bias
    egger_str <- tryCatch({
      egger <- meta::metabias(meta_obj, k.min=4)  # default k.min=10
      if (is.null(egger$estimate) || 
          is.null(egger$estimate["bias"]) || 
          is.null(egger$estimate["se.bias"])) {
        NA
      } else {
        egger_intercept <- egger$estimate["bias"]
        egger_se <- egger$estimate["se.bias"]
        egger_lci <- egger_intercept - 1.96 * egger_se
        egger_uci <- egger_intercept + 1.96 * egger_se
        sprintf("%.2f [%.2f; %.2f], p=%.2f", 
                egger_intercept, egger_lci, egger_uci, egger$pval)
      }
    }, error = function(e) { 
      NA 
    })
    
    # Create a summary tibble row for this biomarker
    tibble(
      fluid = fluid_type,
      'biomarker' = bio,
      combined = paste(fluid_type, "-", bio),
      'Number of cohorts' = meta_obj$k.all,
      'Total sample size' = sum(meta_obj$n),
      !!corr_col_name := sprintf("%.2f [%.2f; %.2f]", 
                                 meta_obj$TE.random,
                                 meta_obj$lower.random,
                                 meta_obj$upper.random),
      `Egger's test` = egger_str
    )
  })
  
  # Remove any NULL results (biomarkers with no data) and combine rows
  results_list <- results_list[!sapply(results_list, is.null)]
  summary_tbl <- bind_rows(results_list)
  
  return(summary_tbl)
}

# For Publication -------------------------------------------------------------

## Make Tables for Publication ----------------------------------------

create_grouped_flextable <- function(data, group_col = "Type") {
  # Get unique group values (preserving order in the data)
  groups <- unique(data[[group_col]])
  
  # Initialize output data frame and vector for storing group sizes
  df_out <- tibble()
  group_lengths <- integer()
  
  for (gr in groups) {
    # Filter by group, drop group_col, and convert all columns to character
    df_group <- data %>%
      filter(.data[[group_col]] == gr) %>%
      select(-all_of(group_col)) %>%
      mutate(across(everything(), as.character))
    
    # Create a heading row: first column is the group name, rest are empty strings
    empty_row <- as_tibble(
      matrix("", nrow = 1, ncol = ncol(df_group)),
      .name_repair = "minimal"
    )
    colnames(empty_row) <- colnames(df_group)
    empty_row[[1]] <- gr
    
    # Combine heading row and group data
    df_out <- bind_rows(df_out, empty_row, df_group)
    group_lengths <- c(group_lengths, nrow(df_group))
  }
  
  # Calculate row indices for each group heading
  heading_indices <- integer()
  idx <- 1
  for (len in group_lengths) {
    heading_indices <- c(heading_indices, idx)
    idx <- idx + 1 + len
  }
  
  # Create flextable and add a top border
  border <- fp_border(color = "#000000", width = 1)
  ft <- flextable(df_out) %>%
    hline_top(part = "body", border = border)
  
  # Add borders above and below each group heading
  for (i in seq_along(heading_indices)) {
    hi <- heading_indices[i]
    
    # Border below the heading row
    ft <- ft %>% hline(i = hi, part = "body", border = border)
    
    # Border above the heading row (starting from the second group)
    if (i > 1) {
      ft <- ft %>% hline(i = hi - 1, part = "body", border = border)
    }
  }
  
  # Adjust column widths automatically
  ft %>% autofit()
}


## Export as .docx file -----------------------------------------------

ExFT <- function(ft,
                 title = "Table Title",
                 file = "output/test.docx",
                 fontname = "Times New Roman",
                 fontsize = 12,
                 landscape = TRUE) {
  # Define custom formatted text (Times New Roman, size 12, bold)
  my_par <- fpar(
    ftext(title, fp_text(font.size = fontsize, font.family = fontname, bold = TRUE))
  )
  
  # Start a new Word document and add the formatted title and table
  doc <- read_docx() %>%
    body_add_fpar(my_par) %>%
    body_add_flextable(ft)
  
  # If landscape mode is requested, add landscape section
  if (landscape) {
    doc <- body_end_section_landscape(doc)
  }
  
  # Save the document
  print(doc, target = file)
}

