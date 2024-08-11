########
#load libraries
########
suppressWarnings(suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(msigdbr); library(viridis); library(ggExtra);
  library(tidyr); library(ggplotify); library(ComplexHeatmap); library(Seurat); library(limma);
  library(ggrastr); library(stringr); library(enrichplot); library(data.table); library(fgsea); library("org.Hs.eg.db"); 
  library(tibble); library(stats); library(ggfortify); library(clusterProfiler); library(edgeR)}))

#setwd 
setwd("/Users/prosperchukwuemeka/Movies/Ziglar/New_CC_arabica")

#@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@
#isolate interesting genesets for KEGG and Reactome and gene ontology analysis
c2.genesets <- msigdbr(species="Homo sapiens", category="C2") 

# #get KEGG LEGACY genesets in c2 from msigdb
# kegg.indices <- grepl("KEGG", c2.genesets$gs_name, ignore.case = TRUE)
Reactome.indices <- grepl("REACTOME", c2.genesets$gs_name, ignore.case = TRUE)

# #Get the rows that match the pattern
# c2.kegg_genesets <- c2.genesets[kegg.indices, ]
c2.reactome_genesets <- c2.genesets[Reactome.indices, ]

#isolate hallmark genesets 
hallmark.genesets <- msigdbr(species="Homo sapiens", category="H") 

#combine both genesets
genesets <- rbind(c2.reactome_genesets, hallmark.genesets)

#number of genesets
length(unique(genesets$gs_name))


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
##########
#read protein targets
##########
cc_carabica.proteins <- c("PLAU", "TPX2", "TPI1", "FABP5", "KIF11", "GPI", "ESRRA", "AURKA", "MTHFD1", 
                          "UMPS", "GM2A", "RANBP1", "TYMS", "GART", "CCNA2", "STAT1", "DOT1L", "MMP12", 
                          "PARP1", "VDR", "HMGCR", "TK1", "MDM2", "SYK", "DTYMK", "HPRT1", "ADK", 
                          "CA2", "CDK2", "ACP3", "SEC14L2", "GALE", "YARS1", "SULT2B1", "TYMP", 
                          "S100A9", "UCK2", "MMP9")

##########
#Run enrichment analysis using cc_carabica.proteins and curated c2 genesets
##########
cc_carabica.proteins.enrichment <- enricher(cc_carabica.proteins, minGSSize = 38, maxGSSize = 300, TERM2GENE = genesets, pAdjustMethod = "BH")
cc_carabica.proteins.enrichment.res <- cc_carabica.proteins.enrichment@result

#subset significant pathway
cc_carabica.proteins.enrichment.res <- subset(cc_carabica.proteins.enrichment.res, subset = p.adjust < 0.05)

# Function to format genes
format_genes <- function(gene_string) {
  # Split the string by "/"
  gene_list <- unlist(strsplit(gene_string, "/"))
  # Combine into the desired format
  formatted_genes <- paste0('"', gene_list, '"', collapse = ", ")
  return(formatted_genes)
}

#Apply the function to the column
cc_carabica.proteins.enrichment.res$format_genes <- sapply(cc_carabica.proteins.enrichment.res$geneID, format_genes)

#calculate fold enrichment
cc_carabica.proteins.enrichment.res <- cc_carabica.proteins.enrichment.res %>%
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

#filter top 20 pathways
cc_carabica.proteins.enrichment.res <- cc_carabica.proteins.enrichment.res %>% 
  slice_max(n = 20, order_by = `Fold Enrichment`) %>% as.data.frame()

#View the result
View(cc_carabica.proteins.enrichment.res)

#save csv
write.csv(cc_carabica.proteins.enrichment.res, file = "cc_carabica_enriched_pathways.csv")


#######
#plot
#######
#plot parameter
themes <- theme(strip.text = element_text(face = "bold", size = 8), panel.border = element_rect(color = "black", fill = NA, size = 1),
                legend.text = element_text(face = "bold", size = 15),legend.title = element_text(face = "bold", size = 13),
                plot.title = element_text(face = "bold", size = 15), axis.title = element_text(face = "bold", size = 15),
                axis.text.x = element_text(face = "bold", size = 10, colour = "black", angle = 15, hjust = 1),
                axis.text.y = element_text(face = "bold", size = 12, colour = "black"),
                panel.background = element_rect(fill = "white"), axis.line = element_line(linewidth = 0.8, colour = "black"))

#plot
plot <- ggplot(cc_carabica.proteins.enrichment.res, aes(x = reorder(ID, `Fold Enrichment`), y = `Fold Enrichment`, fill = -log10(p.adjust))) +
  geom_bar(stat = "identity") +
  coord_flip() + labs(x = "Enriched Genesets", y= "Fold Enrichment",
                      fill = expression(bold("-Log"[10]*"(FDR)"))) +
  geom_text(aes(label = sprintf("%.1f", `Fold Enrichment`), fontface = "bold"),
                                                   position = position_stack(vjust = 0.5),
                                                   size = 5, color = "white") +
  themes + ggtitle("Significant genesets associated with target overlap between CC DEGs and C.arabica") 

#save fig
# Save the plot
ggsave(filename = "cc_carabica_pathway_plot.png", plot = plot, width = 16, height = 8, dpi = 600)


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
# ggplot(cc_carabica.proteins.enrichment.res, aes(x = reorder(ID, `Fold Enrichment`), y = `Fold Enrichment`, fill = -log10(p.adjust), size = Count)) +
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
#   scale_fill_viridis_c() +
#   guides(size= guide_legend(title = "Number of genes")) +
#   ggtitle("Significant genesets associated with target overlap between CC and C.arabica") 
# 
