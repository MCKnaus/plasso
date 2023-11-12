library(showtext)
library(hexSticker)

font_add_google("Inconsolata", "incon")

hexSticker::sticker(
  
  # text
  package = "plasso",
  p_family = "incon",
  p_size = 16,
  p_y = 0.3,
  p_x = 1,
  
  # plot
  subplot = "docs/figures/raw.png",
  s_y = 1.15,
  s_x = 1,
  s_width = 0.7,
  
  # spotlight
  spotlight = TRUE,
  l_y = 1.3,
  l_x = 0.6,
  
  # background
  h_fill = "#EECE63",
  h_color = "#AD890F",
  
  # resolution
  dpi = 1200,
  
  filename = "docs/figures/plasso.png"
)