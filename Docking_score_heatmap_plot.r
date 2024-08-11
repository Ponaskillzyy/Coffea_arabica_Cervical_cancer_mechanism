########
#load libraries
########
#set new wkdir
setwd("/Users/prosperchukwuemeka/Movies/Ziglar/New_CC_arabica")
suppressWarnings(suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(viridis); library(ggExtra); library(ggpubr);
  library(tidyr); library(ggplotify); library(ComplexHeatmap); library(Seurat); library(limma); library(SCpubr); library(pheatmap);
  library(ggrastr); library(stringr); library(enrichplot); library(data.table); library("RColorBrewer");
  library(tibble); library(stats); library(ggfortify); library(devtools); library(ggsignif)}))

#@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@
#load docking score data
ccna2.doc.res <- read.csv(file = "/Users/prosperchukwuemeka/Movies/Ziglar/New_CC_arabica/Docking_scores.csv.csv", sep =",", row.names = 1)

#remove redundant docked compounds
colnames(ccna2.doc.res) <- "CCNA2"

#convert dataframe to matrix
ccna2_docking_res <- as.matrix(ccna2.doc.res)

#@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@
# Define heatmap hyperparameters
heatmap_hp <- ht_opt(
  legend_title_gp = gpar(fontsize = 14, fontface = "bold"), 
  legend_labels_gp = gpar(fontsize = 16, fontface = "bold"), 
  heatmap_column_names_gp = gpar(fontsize = 14, fontface = "bold"),
  heatmap_column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  heatmap_row_title_gp = gpar(fontsize = 14, fontface = "bold"), 
  heatmap_border = TRUE
)

# Define column title
column_title <- "Predicted binding affinity of C.arabica compounds against CCNA2"

# Custom function to add text to each cell
cell_text_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(round(ccna2_docking_res[i, j], 2), x, y, 
            gp = gpar(fontsize = 16, fontface = "bold", col = "white"))
}

# Create heatmap
tx_deg_heatmap <- ComplexHeatmap::Heatmap(
  ccna2_docking_res, 
  name = "Predicted binding affinity (kcal/mol)", 
  raster_device = "png", 
  row_names_gp = gpar(fontsize = 14, fontface = "bold"), 
  cluster_rows = TRUE, 
  row_gap = unit(1, "mm"), 
  border = TRUE, 
  heatmap_legend_param = list(
    legend_direction = "horizontal", 
    title_position = "topcenter", 
    legend_width = unit(8, "cm")
  ), 
  show_row_names = TRUE, 
  row_title_rot = 45, 
  show_column_names = TRUE, 
  column_names_rot = 40, 
  width = unit(120, "mm"), 
  column_title = column_title,
  cluster_columns = FALSE, 
  col = rev(brewer.pal(n = 11, name ="RdBu")), 
  show_row_dend = TRUE,
  cell_fun = cell_text_fun
) + heatmap_hp 

# Plot heatmap
draw(tx_deg_heatmap, heatmap_legend_side = "top")
 
