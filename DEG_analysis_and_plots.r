#NB: TCGA data can be assessed using tcgabiolink or https://xenabrowser.net

###################
#Set work directory
###################
setwd("/Users/prosperchukwuemeka/Movies/Ziglar/New_CC_arabica")

###############
#Load libraries
###############
suppressWarnings(suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(tidyr); 
  library("DESeq2"); library(biomaRt); library(survival);
  library(stringr); library(edgeR); 
  library(survminer); library(data.table);
  library(tidyverse)}))

####################################################
#Download gene expression data of cervical from TCGA
####################################################
#download link
cervical.exp_mat.link <- "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-CESC.htseq_counts.tsv.gz"

#retrieve gene expression data
download.file(cervical.exp_mat.link, destfile= paste(getwd(), "cer_dataset.tsv.gz", sep ="/"))
cervical.exp_mat <- fread(paste0(getwd(), "/cer_dataset.tsv.gz"))

#View(cervical.exp_mat)

#extract Ensembl_ID
Ensembl_ID <- as.data.frame(cervical.exp_mat$Ensembl_ID)
colnames(Ensembl_ID) <- "Ensembl_ID"

##remove Ensembl_ID from columns
cervical.exp_mat <- cervical.exp_mat[,-1]

#transform gene expression data from pseudo-count to raw count matrix
pseudocnt_count <- function(count){
  raw_count <- as.integer(2^(count)-1)
  return(raw_count)
}

#apply function on dataframe to get raw count
cervical.exp_mat <- apply(cervical.exp_mat, 2, FUN = pseudocnt_count)

View(cervical.exp_mat)

#cbind Ensembl_ID with cervical.exp.mat
cervical.exp_mat  <- cbind(Ensembl_ID, cervical.exp_mat) %>% as.data.frame() %>% column_to_rownames(var = "Ensembl_ID")

#View(cervical.exp_mat)

#mapping gene symbol to Ensembl_ID and feature reduction
#remove version ID from rownames in count matrix
rownames(cervical.exp_mat) <- gsub("\\..*","",rownames(cervical.exp_mat))

#create object for annotating features (i.e Ensembl_ID)
mart <- biomaRt::useMart(biomart = "ensembl", dataset =  "hsapiens_gene_ensembl")
transcript.info <- biomaRt::getBM(attributes = c("hgnc_symbol", "ensembl_gene_id","entrezgene_id",
                                                 "transcript_biotype",  "transcript_length"), mart = mart)

# #save transcript information
# saveRDS(transcript.info, paste0(getwd(), "/transcript.info.rds"))

#select protein coding genes only from biomaRt annotation data
transcript.info <- filter(transcript.info, transcript_biotype %in% c("protein_coding"))

#remove duplicates in transcript.info based on "hgnc_symbol", 
transcript.info <- transcript.info %>% distinct(hgnc_symbol, .keep_all = TRUE) #in this case I kept only the first occurrence of duplicates

#annotate "ensembl_ID" with "hgnc_symbol"
annotation_table <- transcript.info[,1:2]

#create a column containing ensembl_ID in brca raw_count 
cervical.exp_mat <- rownames_to_column(cervical.exp_mat)

#rename new column as "ensembl_gene_id"
colnames(cervical.exp_mat)[1] <- "ensembl_gene_id"

#map cervical.exp_mat to created annotation_table by "ensembl_gene_id"
cervical.exp_mat <- left_join(cervical.exp_mat, annotation_table, 
                              by='ensembl_gene_id') #note: this create a column called "hgnc_symbol" at the end of the raw_count df

#rm rows containing NA
cervical.exp_mat <- na.omit(cervical.exp_mat)

#make "hgnc_symbol" the new rownames and delete "ensembl_gene_id" column
rownames(cervical.exp_mat) <- cervical.exp_mat$hgnc_symbol

#delete the column containing "ensembl_gene_id" & "hgnc_symbol"
cervical.exp_mat <- cervical.exp_mat[, !names(cervical.exp_mat) %in% c("ensembl_gene_id", "hgnc_symbol")]


#@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@
#download link for phenotype data
phenoData_link <- "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-CESC.GDC_phenotype.tsv.gz"

#retrieve phenodata
download.file(phenoData_link, destfile= paste(getwd(), "cer_phenoData.tsv.gz", sep ="/"))
cer_phenodata <- fread(paste0(getwd(), "/cer_phenoData.tsv.gz"))

#subset only samples found in cervical.exp_mat
cer_phenodata <- cer_phenodata %>% filter(submitter_id.samples %in% colnames(cervical.exp_mat))

#subset "Normal" $ "primary Tumor"
cer_phenodata <- cer_phenodata %>% 
  filter(sample_type.samples %in% c("Primary Tumor", "Solid Tissue Normal")) %>% as.data.frame() %>% 
  column_to_rownames(var = "submitter_id.samples")

#subset only samples that are "Metastatic" & "Primary Tumor" in cervical.exp_mat
cervical.exp_mat <- cervical.exp_mat[, colnames(cervical.exp_mat) %in% rownames(cer_phenodata)] %>% as.data.frame()

#check if desired sample type is appropriately retrieved
#table(cer_phenodata$sample_type.samples)

#View phenodata 
head(cer_phenodata)

#change level
cer_phenodata$sample_type.samples <- ifelse(cer_phenodata$sample_type.samples == "Solid Tissue Normal",
                                            "Normal tissue", "Primary Tumor")

#@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@
#subset samples in "cervical.exp_mat" that is contained in "cer_phenodata" (i.e samples with phenotype data)
cervical.exp_mat <- cervical.exp_mat[, colnames(cervical.exp_mat) %in% rownames(cer_phenodata)] %>% as.data.frame()

#check dimension of "cervical.exp_mat" & "cer_phenodata"; note that n_samples must be the same in all three dataframe
dim(cervical.exp_mat)
dim(cer_phenodata)

#retrieve transcript_length from biomart
transcript <- data.frame(transcript.info[, c("hgnc_symbol", "transcript_length")])

#make a dataframe of genes present in cervical.exp_mat
hgnc_symbol <- as.data.frame(rownames(cervical.exp_mat), row.names = rownames(cervical.exp_mat))
colnames(hgnc_symbol) <- "hgnc_symbol"

#join hgnc_symbol dataframe $ transcript dataframe to get 
transcript <- left_join(hgnc_symbol, transcript, 
                        by='hgnc_symbol') #%>% column_to_rownames(var = "hgnc_symbol")

rownames(transcript) <- transcript$hgnc_symbol

#View transcript
#head(transcript)

#@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@
#Assume cervical.exp_mat and cer_phenodata are your data objects
# Get the column names of cervical.exp_mat
exp_mat_names <- colnames(cervical.exp_mat)

#Reorder the rows of cer_phenodata to match the column names of cervical.exp_mat
cer_phenodata_ordered <- cer_phenodata[match(exp_mat_names, rownames(cer_phenodata)), ]

# Check if the reordering was successful
head(rownames(cer_phenodata_ordered))

##############################
#Normalization by TMM in edgeR
##############################
#Organize data in a summarized experiment (se) for downstream
cer_se <- SummarizedExperiment(assays = cervical.exp_mat,
                               rowData = rownames(transcript),
                               colData = cer_phenodata_ordered)

# #save summarized experiment as RDS for posterity
# saveRDS(cer_se, "/Users/prosperchukwuemeka/Movies/Ziglar/cer_se.rds")
# 
# #read summarized experiment
# cer_se <- readRDS("cer_se.rds")
# #assay(cer_se)
# #colData(cer_se)
# #rowData(cer_se)

#normalization by trimmed mean of M (TMM)
#create factor for all samples
group <- factor(cer_se$sample_type.samples)

#create dgeList
dge <- DGEList(counts=assay(cer_se), group = group,
               samples=colData(cer_se),
               genes=as.data.frame(rowData(cer_se)))

#re-filter data using edgeR
keep <- rowSums(edgeR::cpm(dge)>100) >= 2 #we're only keeping a gene if it has a cpm of 100 or greater for at least two samples.
dge <- dge[keep, , keep.lib.sizes=FALSE]

#free memory
rm(keep)

#reset library size
dge$samples$lib.size <- colSums(dge$counts)
head(dge$samples)

# Normalization (by TMM)
dge <- calcNormFactors(dge, method="TMM")

# #save dge as RDS
# saveRDS(object = dge,
#         file = "cer_dge.RDS",
#         compress = FALSE)

#read dge
dge <- readRDS("cer_dge.RDS")

#get tmm normalized count
cervicalData_tmm_normalize <- edgeR::cpm(dge, log = FALSE) %>% as.matrix.default()

############
#Primary tumor vs Normal tissue
############
#create contrasts for pair-wise comparisons
design <- model.matrix(~0+group, data = dge$samples)
colnames(design) <- levels(dge$samples$group)
colnames(design) <- c("Normal_Tissue", "Primary_Tumor")
model.contrasts <- makeContrasts("Primary_vs_Primary_Tumor" = Primary_Tumor-Normal_Tissue, levels = design)

#Estimating dispersions for quantile-adjusted conditional maximum likelihood
dge <- estimateDisp(dge, design)

#Fits a negative binomial GLM and creates a DGELM object
fit <- glmQLFit(dge, design)

#Perform Quasi-likelihood model test on the samples according to the mentioned pairs
qlf.Primary_vs_Primary_Tumor <- glmQLFTest(fit, contrast = model.contrasts[,"Primary_vs_Primary_Tumor"])

#extract deg table
res <- qlf.Primary_vs_Primary_Tumor$table

#correct pvalue @ 5% FDR
res$FDR <- p.adjust(res$PValue, method = "fdr")

#add gene column
res$Gene <- rownames(res)

#make expression argument
res$Expression = ifelse(res$FDR < 0.05 & abs(res$logFC) >= 1, 
                        ifelse(res$logFC>= 1 ,'Up-regulated in Tumor','Down-regulated in Tumor'),
                        'Stable')

# #save DEG result
# saveRDS(res, "cc_deg_result.rds")
# 
#readRDS
res <- readRDS("cc_deg_result.rds")

#@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@
# #subset upregulated genes
# up_in_cc <- res %>% 
#   filter(Expression == "Up-regulated in Tumor")
# 
# #arrange based on logfold change
# up_in_cc <- up_in_cc %>% arrange(desc(logFC))
# 
# #write csv
# write.csv(up_in_cc, "up_in_cc.csv")

#@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@
#############
#VOLCANO PLOT 
#############
#volcano plot hyper parameters
volcano.theme <- theme(strip.text = element_text(face = "bold", size = 20),
                       strip.background = element_blank(),
                       panel.grid.major = element_blank(), 
                       panel.border = element_rect(colour = "black", fill=NA, size=1),
                       panel.grid.minor = element_blank(), 
                       legend.text = element_text(face = "bold", size = 17),
                       legend.title = element_text(face = "bold", size = 20),
                       plot.title = element_text(face = "bold", size = 20),
                       axis.title = element_text(face = "bold", size = 20),
                       axis.text.x = element_text(face = "bold", size = 20, colour = "black"),
                       axis.text.y = element_text(face = "bold", size = 20, colour = "black"),
                       panel.background = element_rect(fill = "white"),
                       axis.line = element_line(linewidth = 0.8, colour = "black"))

#remove genes with NA as adjpvalue
res <- na.omit(res)

#Apply a small offset to p-value
offset <- 0 + 1e-50 #prevent having a pvalue of exactly zero since log10(0) is undefined

#index genes with padj with exactly 0
ndx.0.padj <- which(res$FDR==0)
res[ndx.0.padj,] <- offset

#Plot with ggplot
res.plot <- ggplot(data = res, 
                   aes(x = logFC, 
                       y = -log(FDR, 10), 
                       colour= Expression)) + 
  geom_point(size=3.5) +
  xlim(c(-6, 6)) +
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log(0.05, 10),lty=4,col="black",lwd=0.8) +
  labs(x=expression(bold("Log"[2]*"(Fold Change)")),
       y= expression(bold("-Log"[10]*"(FDR)")),
       title="DEGs in Tumor samples vs Normal Tissue in Cervical Cancer")  + 
  theme_classic() + volcano.theme + theme(legend.position = "top")

#view
print(res.plot)

#add labels and connectors for top enriched genes
top.genes <- 50
top_genes <- rbind(
  res %>%
    filter(Expression == 'Up-regulated in Tumor') %>%
    arrange(FDR, desc(abs(logFC))) %>%
    head(top.genes),
  res %>%
    filter(Expression == 'Down-regulated in Tumor') %>%
    arrange(FDR, desc(abs(logFC))) %>%
    head(top.genes)
)

#add gene name to column
top_genes$gene_name <- rownames(top_genes)

#show gene of interest
res <- res.plot + 
  geom_label_repel(data = top_genes,
                   mapping = aes(logFC, -log(FDR,10), label = Gene, fontface = "bold"), direction = "both",
                   size = 3, nudge_x = 0.1, box.padding = 0.5,  show.legend = FALSE, max.overlaps = Inf,
                   nudge_y = 0.2) + guides(color = guide_legend(override.aes = list(size = 5))) 

# Display the plot
print(res) 


########
#PLOT 2
########
#readRDS
res <- readRDS("cc_deg_result.rds")

#select core compound targeted genes in cervical cancer
compound_disease_signature <- c("PLAU", "TPX2", "TPI1", "FABP5", "KIF11", "GPI", "ESRRA", "AURKA", "MTHFD1", 
                                "UMPS", "GM2A", "RANBP1", "TYMS", "GART", "CCNA2", "STAT1", "DOT1L", "MMP12", 
                                "PARP1", "VDR", "HMGCR", "TK1", "MDM2", "SYK", "DTYMK", "HPRT1", "ADK", 
                                "CA2", "CDK2", "ACP3", "SEC14L2", "GALE", "YARS1", "SULT2B1", "TYMP", 
                                "S100A9", "UCK2", "MMP9")

#remove genes with NA as adjpvalue
res <- na.omit(res)

#Apply a small offset to p-value
offset <- 0 + 1e-50 #prevent having a pvalue of exactly zero since log10(0) is undefined

#index genes with padj with exactly 0
ndx.0.padj <- which(res$FDR==0)
res[ndx.0.padj,] <- offset

#Plot with ggplot
res.plot <- ggplot(data = res, 
                   aes(x = logFC, 
                       y = -log(FDR, 10), 
                       colour= Expression)) + 
  geom_point(size=3.5) +
  xlim(c(-6, 6)) +
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log(0.05, 10),lty=4,col="black",lwd=0.8) +
  labs(x=expression(bold("Log"[2]*"(Fold Change)")),
       y= expression(bold("-Log"[10]*"(FDR)")),
       title="Overlapping gene targets between cervical cancer DEGs and C. arabica compound-derived targets")  + 
  theme_classic() + volcano.theme + theme(legend.position = "top")

#view
print(res.plot)

#add labels and connectors for top enriched genes
compound_disease_signature_length <- 38

#extract the data of compound targeted genes in cervical cancer
compound_disease_signature.dat <- res %>% 
  filter(res$Gene %in% compound_disease_signature)

#show gene of interest
res <- res.plot + 
  geom_label_repel(data = compound_disease_signature.dat,
                   mapping = aes(logFC, -log(FDR,10), label = Gene, fontface = "bold"), direction = "both",
                   size = 3, nudge_x = 0.1, box.padding = 0.5,  show.legend = FALSE, max.overlaps = Inf,
                   nudge_y = 0.2) + guides(color = guide_legend(override.aes = list(size = 5))) 

# Display the plot
print(res) 

#@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@
#ndx for stratifying patients based on tissue biopsy
normal_tissue_idents <- rownames(dge$samples)[which(dge$samples$sample_type.samples == "Normal tissue")]
primary_tumor_idents <- rownames(dge$samples)[which(dge$samples$sample_type.samples == "Primary Tumor")]

#get expression matrix
expression_data <- cervicalData_tmm_normalize

#select "CCNA2" for statistical analysis
cervical_cancer_biomarker <- c("CCNA2")

#subset expression data for gene of interest
expression_data <- as.data.frame(expression_data[rownames(cervicalData_tmm_normalize) %in% cervical_cancer_biomarker,])
names(expression_data) <- "CCNA2"

#get tissue details from phenotype data
phenotype_data <- dge$samples %>% dplyr::select(sample_type.samples)
names(phenotype_data) <- "Tissue"

#match patient ident in expression and pData
phenotype_data <- phenotype_data[match(rownames(expression_data), rownames(phenotype_data)), ]

#Combine the data frames by column
expression_pData <- cbind(expression_data, phenotype_data)

#rename column
names(expression_pData)[2] <- "Tissue"

#####
#plot graph with statistical significance
#####
#Define your specific comparisons
comparisons_list <- list(c("Normal tissue", "Primary Tumor"))

#themes
themes <- theme(strip.text = element_text(face = "bold", size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1),
                legend.text = element_text(face = "bold", size = 15),legend.title = element_text(face = "bold", size = 15),
                plot.title = element_text(face = "bold", size = 15), axis.title = element_text(face = "bold", size = 15),
                axis.text.x = element_text(face = "bold", size = 15, colour = "black", angle = 40, hjust = 1), 
                axis.text.y = element_text(face = "bold", size = 15, colour = "black"),
                axis.ticks.x = element_blank(), legend.position = "None",
                panel.background = element_rect(fill = "white"), axis.line = element_line(linewidth = 0.8, colour = "black"))

expression_pData %>%
  ggplot(., aes(x=Tissue, y= CCNA2)) +
  geom_boxplot(aes(fill = Tissue), outlier.color = NA) + theme_classic() +
  #geom_jitter(size = 0.2, alpha = 0.5) +
  stat_compare_means(comparisons = comparisons_list, method = "wilcox.test", 
                     label = "p.signif", p.adjust.method = "BH", 
                     vjust = 0, 
                     map_signif_level = TRUE, step_increase = 0.1, 
                     textsize = 6, size = 7, 
                     tip_length = 5, na.rm = T) +
  labs(y= "Expression log(cpm)", x = " ", fill = "Tissue") +
  theme_classic() +
  themes +
  ggtitle("CCNA2 expression in primary tumor vs normal tissue")
