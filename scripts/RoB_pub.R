
# RoB traffic lights ------------------------------------------------------

source("scripts/RoB_func.R")

# font_import(prompt = FALSE)
# loadfonts(device = "pdf")
# fonts()


# QUADAS2 -----------------------------------------------------------------

quadas_final <- read_excel("data/Rob_pub.xlsx", sheet = "QUADAS2") %>% arrange(desc(Study))

# Create the plot
QUADASp <- make_traffic_light(
  df            = quadas_final,
  id_col        = "Study",         # Use the Study column for labels
  domain_cols   = domains_q2,
  domain_labels = domain_labels,
  annotations   = annotations,
  circle_size   = 6,               # Adjust circle size
  text_size     = 3,               # Adjust symbol size
  caption_size  = 9,               # Adjust caption text size
  caption_margin= 25               # Increase bottom margin
)

# Display the plot
print(QUADASp)

# QUADASp <- QUADASp + theme(text = element_text(family = "Helvetica"))
# QUADASp


# ggsave(
#   filename = "output/QUADAS2_fig.pdf",
#   plot     = QUADASp,
#   device   = cairo_pdf,
#   width    = 8,      # inches
#   height   = 20,     # make this as tall as you need
#   units    = "in",   # units for width/height
#   dpi      = 300
# )


# QUIPS -------------------------------------------------------------------

quips_final <- read_excel("data/Rob_pub.xlsx", sheet = "QUIPS") 

# Create the plot
QUIPSp <- make_traffic_light(
  df            = quips_final,
  id_col        = "Study",
  domain_cols   = domains_quips,
  domain_labels = domain_labels_quips,
  annotations   = annotations_quips,
  # QUIPS style
  fill_chr      = "Risk of Bias",   
  risk_levels   = risk_levels_quips,
  symbols       = symbols_quips,
  fill_colors   = fill_colors_quips,
  # Size modeification
  circle_size   = 6,
  text_size     = 3,
  caption_size  = 10.5,
  caption_margin= 25
)

print(QUIPSp)


# QUIPSp <- QUIPSp + theme(text = element_text(family = "Helvetica"))
# QUIPSp

# save as PDF file
# ggsave(
#   filename = "output/QUIPS_fig.pdf",
#   plot     = QUIPSp,
#   device   = cairo_pdf,
#   width    = 8,
#   height   = 20,
#   units    = "in",
#   dpi      = 300
# )

















