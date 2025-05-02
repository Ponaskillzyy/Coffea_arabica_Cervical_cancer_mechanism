########
#load libraries
########
suppressWarnings(suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(msigdbr); library(viridis); library(ggExtra);
  library(tidyr); library(ggplotify); library(ComplexHeatmap); library(Seurat); library(limma);
  library(ggrastr); library(stringr); library(enrichplot); library(data.table); library(fgsea); library("org.Hs.eg.db"); 
  library(tibble); library(stats); library(ggfortify); library(clusterProfiler); library(edgeR)}))

#setwd 
setwd("/Users/prosperchukwuemeka/Movies/Ziglar/New_CC_arabica/Revision")

#@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@
#isolate interesting genesets for KEGG and Reactome and gene ontology analysis
c2.genesets <- msigdbr(species="Homo sapiens", category="C2") 

#get KEGG LEGACY and Reactome genesets in c2 from msigdb
kegg.indices <- grepl("KEGG", c2.genesets$gs_name, ignore.case = TRUE)
Reactome.indices <- grepl("REACTOME", c2.genesets$gs_name, ignore.case = TRUE)

##########
#Get the rows that match the pattern
##########
#kegg
c2.kegg_genesets <- c2.genesets[kegg.indices, ]
length(unique(c2.kegg_genesets$gs_name))
#total kegg gse == "186"

#reactome
c2.reactome_genesets <- c2.genesets[Reactome.indices, ]
length(unique(c2.reactome_genesets$gs_name))
#total reactome gse == "1615"

#isolate hallmark genesets 
hallmark.genesets <- msigdbr(species="Homo sapiens", category="H") 
length(unique(hallmark.genesets$gs_name))
#total gse == "50"

#combine both genesets
genesets <- rbind(c2.reactome_genesets, c2.kegg_genesets, hallmark.genesets)

#number of genesets
length(unique(genesets$gs_name))
#total gse == "1851"

#print
print(genesets)

#make each genes in each geneset unique
genesets <- genesets %>% group_by(gs_name) %>% summarise(gene_symbol = unique(gene_symbol))

#Function to format row names
format_row_names <- function(names) {
  names <- str_to_upper(gsub("_", " ", names))   # Replace underscores with spaces and Convert to title case
  return(names)
}

#Apply the formatting function to names of geneset
genesets$gs_name <- format_row_names(genesets$gs_name)

# #Remove the first word from each row name
# c2.genesets$gs_name <- gsub("^[^ ]+ ", "", c2.genesets$gs_name)

#@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@
#load compound-disease shared targets
compound_disease_targets <- fread("Common Targets.csv") %>% 
  pull(`Common Targets`) 
  
#exclude specific genes
compound_disease_targets <- setdiff(compound_disease_targets, c("AKR1C1", "AKR1C2", "NOS2"))

##########
#Run enrichment analysis using compound_disease_targets and curated genesets from kegg, reactome, and hallmark database
##########
compound_disease_targets.enrichment <- enricher(compound_disease_targets, minGSSize = 10, maxGSSize = 500, TERM2GENE = genesets, pAdjustMethod = "BH")
compound_disease_targets.enrichment.res <- compound_disease_targets.enrichment@result

#subset significant pathway
compound_disease_targets.enrichment.res <- subset(compound_disease_targets.enrichment.res, subset = p.adjust < 0.05)

# Function to format genes
format_genes <- function(gene_string) {
  # Split the string by "/"
  gene_list <- unlist(strsplit(gene_string, "/"))
  # Combine into the desired format
  formatted_genes <- paste0('"', gene_list, '"', collapse = ", ")
  return(formatted_genes)
}

#Apply the function to the column
compound_disease_targets.enrichment.res$format_genes <- sapply(compound_disease_targets.enrichment.res$geneID, format_genes)

#calculate fold enrichment
compound_disease_targets.enrichment.res <- compound_disease_targets.enrichment.res %>%
  separate(GeneRatio, into = c("GeneNum", "GeneDenom"), sep = "/") %>%
  separate(BgRatio, into = c("BgNum", "BgDenom"), sep = "/") %>%
  mutate(
    #Convert character columns to numeric
    GeneNum = as.numeric(GeneNum),
    GeneDenom = as.numeric(GeneDenom),
    BgNum = as.numeric(BgNum),
    BgDenom = as.numeric(BgDenom),
    
    #Calculate Gene Ratio and Background Ratio
    GeneRatioValue = GeneNum / GeneDenom,
    BgRatioValue = BgNum / BgDenom,
    
    # Calculate Fold Enrichment
    `Fold Enrichment` = GeneRatioValue / BgRatioValue
  ) 

#View the result
View(compound_disease_targets.enrichment.res)

#save csv
write.csv(compound_disease_targets.enrichment.res, file = "cc_carabica_enriched_pathways.csv")

#######
#plot
#######
#filter top 20 pathways
compound_disease_targets.enrichment.res <- compound_disease_targets.enrichment.res %>% 
  slice_max(n = 20, order_by = `Fold Enrichment`) %>% as.data.frame()

# Use RdBu color scale from RColorBrewer
color_scale <- scale_fill_distiller(
  palette = "RdBu",  # Use the RdBu palette
  direction = -1,    # Reverse the palette if needed
  name = expression(bold("-Log"[10]*"(FDR)"))  # Custom legend title
)

# Your ggplot code
ggplot(compound_disease_targets.enrichment.res, aes(x = reorder(rownames(compound_disease_targets.enrichment.res), `Fold Enrichment`), y = `Fold Enrichment`, fill = -log10(p.adjust))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "Enriched Genesets", y= "Fold Enrichment",
       fill = expression(bold("-Log"[10]*"(FDR)"))) +
  geom_text(aes(label = format_genes, fontface = "bold"),
            position = position_stack(vjust = 0.5),
            size = 3.2, color = "black") +
  ggtitle("Significant Genesets Associated with Shared Compound-Cervical Cancer Targets") +
  color_scale + theme_classic() + themes

# #@@@@@@@@@@@@@@@@@@@@
# #@@@@@@@@@@@@@@@@@@@@
# ##################
# ## ggplot theme ##
# ##################
# theme_mrl <- function(x = 1) {
#   theme_minimal() +
#     theme(
#       axis.line = element_line(),
#       axis.ticks.x = element_line(),
#       axis.ticks.y = element_line(),
#       axis.text.x = element_text(size = 12*x,face = "bold", angle = 0, vjust = 0.6),
#       axis.text.y = element_text(size = 12*x,face = "bold"),
#       axis.title.x = element_text(size = 12*x,face = "bold"),
#       axis.title.y = element_text(size = 12*x,face = "bold"),
#       strip.background = element_rect(fill="gray20", colour="gray20", linetype="solid"),
#       strip.text = element_text(size=14*x, colour="white", face="bold"),
#       legend.title = element_text(size=14*x, face = "bold"),
#       legend.text = element_text(size=12*x, color="gray20", face="bold"),
#       legend.background = element_rect(fill = "transparent", colour = "transparent"),
#       plot.title =  element_text(hjust=0.5, vjust=2, face="bold"),
#       plot.subtitle = element_text(hjust=0.5, vjust=3, face="italic"),
#       plot.caption = element_text(hjust = 0, face = "italic")
#     )
# }
# 
# 
# #plot
# ggplot(compound_disease_targets.enrichment.res, aes(x = reorder(ID, `Fold Enrichment`), y = `Fold Enrichment`, fill = -log10(p.adjust), size = Count)) +
#   geom_point(shape = 21, stroke = 0.1) +
#   coord_flip() +
#   theme_mrl(1.3) +
#   theme(strip.text = element_text(face = "bold", size = 8), panel.border = element_rect(color = "black", fill = NA, size = 1),
#         legend.text = element_text(face = "bold", size = 15),legend.title = element_text(face = "bold", size = 13),
#         plot.title = element_text(face = "bold", size = 15), axis.title = element_text(face = "bold", size = 15),
#         axis.text.x = element_text(face = "bold", size = 10, colour = "black", angle = 15, hjust = 1),
#         axis.text.y = element_text(face = "bold", size = 12, colour = "black"),
#         aspect.ratio = 2, legend.key.size = unit(0.4, "cm"),
#         panel.background = element_rect(fill = "white"), axis.line = element_line(linewidth = 0.8, colour = "black"))+
#   scale_x_discrete(position = "bottom",labels = function(x) str_wrap(x, width = 30)) +
#   labs(x = "Enriched Genesets", y = "Fold Enrichment", fill = expression(bold("-Log"[10]*"(FDR)"))) +
#   color_scale + #scale_fill_viridis() +
#   guides(size= guide_legend(title = "Number of genes")) +
#   ggtitle("Significant Genesets Associated with Shared Compound-Cervical Cancer Targets")
# 
