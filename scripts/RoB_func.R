library(tidyverse)
library(ggplot2)

# Function ----------------------------------------------------------------

make_traffic_light <- function(
    df,
    id_col         = "Study",       # Column name to use as the study identifier
    domain_cols,                   # Character vector of column names representing domains
    domain_labels  = domain_cols,  # Labels to display on the x-axis for each domain
    annotations     = NULL,        # Character vector of notes to include as a caption
    risk_levels     = c("LOW", "UNCLEAR", "HIGH"),
    symbols         = c(
      "LOW"     = "+",
      "UNCLEAR" = "-",
      "HIGH"    = "×"
    ),
    fill_colors     = c(
      "LOW"     = "#4CAF50",  # green
      "UNCLEAR" = "#FFEB3B",  # yellow
      "HIGH"    = "#F44336"   # red
    ),
    fill_chr = "Assessment",
    circle_size     = 8,        # Size of the circle (shape = 21)
    text_size       = 4,        # Size of the symbol inside the circle
    caption_size    = 10,       # Caption text size (in points)
    caption_margin  = 20        # Bottom margin to make room for the caption (in points)
) {
  # 1) Convert to long format
  df_long <- df %>%
    select(all_of(c(id_col, domain_cols))) %>%
    pivot_longer(
      cols      = -all_of(id_col),
      names_to  = "Domain",
      values_to = "Risk"
    ) %>%
    mutate(
      Domain = factor(Domain, levels = domain_cols),
      Study  = factor(.data[[id_col]], levels = unique(.data[[id_col]])),
      Risk   = factor(Risk, levels = risk_levels)
    )
  
  # 2) Build the plot
  p <- ggplot(df_long, aes(x = Domain, y = Study)) +
    # White tiles with grey borders as the grid background
    geom_tile(colour = "grey80", fill = "white") +
    # Colored circles inside each tile
    geom_point(aes(fill = Risk),
               shape  = 21,
               size   = circle_size,
               colour = "grey70") +
    # Symbol inside each circle
    geom_text(aes(label = symbols[Risk]),
              size     = text_size,
              colour   = "grey10",
              fontface = "bold") +
    # Replace x-axis tick labels with custom domain labels
    scale_x_discrete(labels = domain_labels, position="top") +
    # Use custom fill colors
    scale_fill_manual(values = fill_colors) +
    # Remove axis labels and set legend title
    labs(x = NULL, y = NULL, fill = fill_chr) +
    # Apply a clean theme and adjust margins for the caption
    theme_minimal(base_size = 14, base_family="MyHelv") +
    theme(
      text                   = element_text(family="MyHelv"),
      axis.text.x            = element_text(angle = 0, hjust = 0.5, face = "bold"),
      axis.text.y            = element_text(size = 10),
      legend.position        = "bottom",
      panel.grid             = element_blank(),
      plot.margin            = margin(5, 5, caption_margin, 5),
      plot.caption.position  = "plot",
      plot.caption           = element_text(
        hjust = 0,
        vjust = 1,
        size  = caption_size,
        lineheight = 1.15
      )
    )
  
  # 3) Add caption if annotations are provided
  if (!is.null(annotations)) {
    p <- p + labs(caption = paste(annotations, collapse = "\n"))
  }
  
  return(p)
}




# QUADAS2 df ---------------------------------------------------------------

# Define the original column names for the 7 QUADAS-2 domains
domains_q2 <- c(
  "RoB_PatientSelection",
  "RoB_IndexTest",
  "RoB_ReferenceStandard",
  "RoB_FlowAndTiming",
  "Applicability_PatientSelection",
  "Applicability_IndexTest",
  "Applicability_ReferenceStandard"
)

# Define the labels to show along the x-axis
domain_labels <- c(
  "RoB 1", "RoB 2", "RoB 3", "RoB 4",
  "App 1", "App 2", "App 3"
)

# Define the notes to appear as a caption at the bottom
annotations <- c(
  "RoB 1: Patient Selection - Risk of Bias",
  "RoB 2: Index Test - Risk of Bias",
  "RoB 3: Reference Standard - Risk of Bias",
  "RoB 4: Flow & Timing - Risk of Bias",
  "App 1: Patient Selection - Concerns regarding applicability",
  "App 2: Index Test - Concerns regarding applicability",
  "App 3: Reference Standard - Concerns regarding applicability"
)


# QUIPS df ----------------------------------------------------------------

domains_quips <- c(
  "1_StudyParticipation",
  "2_StudyAttrition",
  "3_PF_measurement",
  "4_Outcome_Measurement",
  "5_Study_Confounding",
  "6_StatAnalysis_Reporting"
)

# domain_labels_quips <- c(
#   "Study\nParticipation",
#   "Study\nAttrition",
#   "PF\nMeasurement",
#   "Outcome\nMeasurement",
#   "Study\nConfounding",
#   "Stat Analysis\nReporting"
# )

domain_labels_quips <- c(
  "Domain 1",
  "Domain 2",
  "Domain 3",
  "Domain 4",
  "Domain 5",
  "Domain 6"
)

#risk levels
risk_levels_quips <- c("LOW", "MODERATE", "HIGH", "UNCLEAR")

#marks 
symbols_quips <- c(
  "LOW"      = "+",
  "MODERATE" = "±",
  "HIGH"     = "×",
  "UNCLEAR"  = "?"
)

# colors
fill_colors_quips <- c(
  "LOW"      = "#4CAF50",  
  "MODERATE" = "#FFC107",  
  "HIGH"     = "#F44336",  
  "UNCLEAR"  = "#9E9E9E"   
)

# caption
annotations_quips <- c(
  "Domain 1: Study participation",
  "Domain 2: Study attrition",
  "Domain 3: Prognostic factor measurement",
  "Domain 4: Outcome measurement",
  "Domain 5: Study confounding",
  "Domain 6: Statistical analysis & reporting"
)

