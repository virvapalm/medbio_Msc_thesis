### Libraries ----
library(dplyr)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(viridis)
library(parallel)
library(docstring)
library(ggrepel)
library(patchwork)
### Set up ----

## ggplot theme
plot_theme <-  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))



#### Parallel env ----
# Core numbers
mc.cores <- detectCores()*.5 %>% round(1)

#### Set directory ----

# directory frpm local file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))

# Main directory 
main_dir <- getwd()


# Output path
out_path <- paste(main_dir, "/Export/", sep = "")
# Figure directory
fig_path <-  paste(main_dir, "/Fig/", sep = "")

# Create directory if not exist
check_dir <- function(path) if(!dir.exists(file.path(path))) dir.create(file.path(path))
sapply(list(out_path, fig_path, pdf_path), check_dir)

### functions ----
# GGplot save output
save_fn <- function(name, h =2, w = 2, ratio = 10, dpi = 300){
  ggsave( paste(fig_path, name, ".png", sep = ""), 
          height = ratio*h, width = ratio*w, limitsize = FALSE,
          dpi = dpi, units = "cm", device = "png")}

