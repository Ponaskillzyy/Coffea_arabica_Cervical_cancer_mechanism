#NB: TCGA data can be assessed using tcgabiolink or https://xenabrowser.net

###################
#Set work directory
###################
setwd("/Users/prosperchukwuemeka/Movies/Ziglar/New_CC_arabica/Revision")

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
cervical.exp_mat.link <- "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-CESC.star_counts.tsv.gz"

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

#view
View(cervical.exp_mat)

#cbind Ensembl_ID with cervical.exp.mat
cervical.exp_mat  <- cbind(Ensembl_ID, cervical.exp_mat) %>% as.data.frame() %>% column_to_rownames(var = "Ensembl_ID")

#view
View(cervical.exp_mat)

#convert dataframe to matrix
cervical.exp_mat  <- as.matrix(cervical.exp_mat)

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

#convert gene expression data type from matrices to dataframe
cervical.exp_mat <- as.data.frame(cervical.exp_mat)

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
phenoData_link <- "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-CESC.clinical.tsv.gz"

#retrieve phenodata
download.file(phenoData_link, destfile= paste(getwd(), "cer_phenoData.tsv.gz", sep ="/"))
cer_phenodata <- fread(paste0(getwd(), "/cer_phenoData.tsv.gz"))

#subset only samples found in cervical.exp_mat
cer_phenodata <- cer_phenodata %>% filter(sample %in% colnames(cervical.exp_mat))

#########
#metadata to look into "disease_type", "age_at_index.demographic", "treatment_or_therapy.treatments.diagnoses"
#########
#subset "primary Tumor"
cer_phenodata <- cer_phenodata %>% 
  filter(sample_type.samples %in% c("Primary Tumor")) %>% as.data.frame() %>% 
  column_to_rownames(var = "sample")

#subset only samples that are "Primary Tumor" in cervical.exp_mat
cervical.exp_mat <- cervical.exp_mat[, colnames(cervical.exp_mat) %in% rownames(cer_phenodata)] %>% as.data.frame()

#check if desired sample type is appropriately retrieved
#table(cer_phenodata$sample_type.samples)

#View phenodata 
head(cer_phenodata)

#@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@
#retrieve survival data 
survData_link <- "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-CESC.survival.tsv.gz"

#retrieve SurvData
download.file(survData_link, destfile= paste(getwd(), "cer_survData.tsv.gz", sep ="/"))
cer_survdata <- fread(paste0(getwd(), "/cer_survData.tsv.gz"))

#convert overall survival time from days to month
cer_survdata <- cer_survdata %>% mutate(OS.time.year = OS.time / 365) %>%
  as.data.frame() %>% column_to_rownames(var = "sample")

#subset samples in "cervical.exp_mat" and "cer_phenodata" that is contained in "cer_survdata" (i.e samples with survival details)
cervical.exp_mat <- cervical.exp_mat[, colnames(cervical.exp_mat) %in% rownames(cer_survdata)] %>% as.data.frame()
cer_phenodata <- cer_phenodata %>% filter(rownames(cer_phenodata) %in% rownames(cer_survdata)) %>% as.data.frame()

#subset samples in "cer_survdata" that are contained in "cervical.exp_mat"
cer_survdata <- cer_survdata %>% filter(rownames(cer_survdata) %in% colnames(cervical.exp_mat)) %>% as.data.frame()

#check dimension of "cervical.exp_mat", "cer_phenodata" & "cer_survdata"; note that n_samples must be the same in all three dataframe
dim(cervical.exp_mat)
dim(cer_phenodata)
dim(cer_survdata)

#add a reference column for joining cer_phenodata and cer_survdata
cer_phenodata[["sample"]] <- rownames(cer_phenodata)
cer_survdata[["sample"]] <- rownames(cer_survdata)

#create a new dataframe containing phenotype & survival information
new_cer_phenodata <- left_join(cer_phenodata, cer_survdata, 
                               by= "sample") %>% column_to_rownames(var = "sample") 

#Reordering rownames of cer_phenodata to match cervical.exp_mat; NB: this should be done to prevent error creating a summarized experiment
reorder_ndx <- match(colnames(cervical.exp_mat), rownames(new_cer_phenodata))
new_cer_phenodata <- new_cer_phenodata[reorder_ndx, ]

#View new phenodata
#head(new_cer_phenodata)

##############################
#Normalization by TMM in edgeR
##############################
#Organize data in a summarized experiment (se) for downstream
cer_se <- SummarizedExperiment(assays = cervical.exp_mat,
                               rowData = rownames(cervical.exp_mat),
                               colData = new_cer_phenodata)

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
# 
# #read dge
# dge <- readRDS("cer_dge.RDS")

#get tmm normalized count
cervicalData_tmm_normalize <- edgeR::cpm(dge, log = FALSE) %>% as.matrix.default()

#View normalized count matrix
View(cervicalData_tmm_normalize)

###############################################
#subset Shared Compound-Cervical Cancer Targets for survival analysis
###############################################
#load compound-disease shared targets
compound_disease_targets <- fread("Common Targets.csv") %>% 
  pull(`Common Targets`) 

#exclude specific genes
compound_disease_targets <- setdiff(compound_disease_targets, c("AKR1C1", "AKR1C2", "NOS2"))

#subset genes
survival_genes <- cervicalData_tmm_normalize[rownames(cervicalData_tmm_normalize) %in% compound_disease_targets,]

#transpose survival_genes dataframe
survival_genes <- as.data.frame(t(survival_genes))

#subset overall survival (OS) and overall survival time (OS.time.year) from summarized experiment
new_survival_info <- as.data.frame(colData(cer_se)[, c("OS","OS.time.year", "initial_weight.samples")], row.names = rownames(colData(cer_se)))

#rename column name in new_survival_info to fit survminer package
colnames(new_survival_info) <- c(OS = "event", OS.time.year = "time", initial_weight.samples = "Initial_weight_of_samples")

#create the dataframe for survival analysis
survival_df <- cbind(survival_genes, new_survival_info)

#save survival_df as rds 
saveRDS(survival_df, "/Users/prosperchukwuemeka/Movies/Ziglar/New_CC_arabica/Revision/survival_df.rds")

#@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@
#read survival_df
survival_df <- readRDS("survival_df.rds")

#View survival genes
head(survival_df)

#survival analysis
#determine the optimal cut-point of genes
surv.cut <- surv_cutpoint(survival_df, time = "time", event = "event",
                          variables = colnames(survival_df))

#view cut-point summary
summary(surv.cut)

#categorize expression level of genes based on the optimal cut-point
surv.cat <- surv_categorize(surv.cut)
surv.cat$Initial_weight_of_samples <- ifelse(surv.cat$Initial_weight_of_samples == "high", "large", "small")

#view categories
surv.cat

#@@@@@@@@@@@
#@@@@@@@@@@@
#get shared compound-disease targets
genes <- compound_disease_targets

#Loop through each gene
for (gene in genes) {
  
  #Fit the survival curve for the current gene
  surv.fit <- survfit(Surv(time, event) ~ get(gene), data = surv.cat)
  
  #Generate the survival plot
  surv.plot <- ggsurvplot(surv.fit, data = surv.cat, risk.table = F, conf.int = T, surv.median.line = "hv", 
                          pval = T, legend.title = gene, legend.labs = c("High", "Low")) +
    labs(x = "Time (Years)", y = "Survival probability")
  
  #Extract the plot object from the ggsurvplot result
  surv.plot <- surv.plot$plot
  
  #Set plot parameters
  surv.plot <- ggpar(surv.plot,
                     font.main = c(16, "bold"),
                     font.x = c(16, "bold"),
                     font.y = c(16, "bold"),
                     font.caption = c(16, "bold"),
                     font.legend = c(16, "bold"),
                     font.tickslab = c(16, "bold"))
  
  #Save the plot with a dynamic filename
  ggsave(paste0(gene, "_survival_plot.png"), plot = surv.plot, width = 15, height = 12, units = "cm")
}

#@@@@@@@@@@@
#@@@@@@@@@@@
######
#multivariate cox hazard analysis
######
#make surv cat copy
surv.cat.copy <- surv.cat

#make controls for cox prop hazard
colnames(surv.cat.copy)[colnames(surv.cat.copy) == "LCK"] <- "LCK (Control 1)"
colnames(surv.cat.copy)[colnames(surv.cat.copy) == "STAT1"] <- "STAT1 (Control 2)"

#genes that correlated with poor outcome
poor.out.genes <- c("LCK (Control 1)", "STAT1 (Control 2)", "CA2", "MET", "MMP3", "CCNA2", "GART", "VDR", "MMP13",
                    "ADAM17", "FABP4", "MMP12", "PARP1", "MMP1", "PLAU", "MMP7")


#Loop through each gene to relevel it 
for (gene in poor.out.genes) {
  surv.cat.copy[[gene]] <- factor(surv.cat.copy[[gene]], levels = c("low", "high"))
}

# #relevel
# surv.cat$Initial_weight_of_samples <- factor(surv.cat$Initial_weight_of_samples,
#                                              levels = c("small", "large"))

#Combine gene names with "+"
gene_combined <- paste(poor.out.genes, collapse = " + ")

#Print the result
print(gene_combined)

#run model
multi.cox_model <- coxph(Surv(time, event) ~ `LCK (Control 1)` + `STAT1 (Control 2)` + CA2 + MET + MMP3 + 
                           CCNA2 + GART + VDR + MMP13 + ADAM17 + FABP4 + MMP12 + PARP1 + MMP1 + PLAU + MMP7,
                         data = surv.cat.copy)

# # Calculate hazard ratio
# cox_summary.multi <- summary(multi.cox_model)
# hazard_ratio.multi <- exp(cox_summary.multi$coefficients) %>% as.data.frame()
# hazard_ratio.multi <- hazard_ratio.multi$`exp(coef)`

#plot multivariate cox result
plot <- ggforest(multi.cox_model, fontsize = 1.3, refLabel = "reference") 

#save fig
# Save the plot
ggsave(filename = "/Users/prosperchukwuemeka/Movies/Ziglar/New_CC_arabica/cc_carabica_cox_plot.png", plot = plot, width = 16, height = 8, dpi = 600)
