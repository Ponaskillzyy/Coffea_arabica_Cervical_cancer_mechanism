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

#View(cervical.exp_mat)

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

#save transcript information
saveRDS(transcript.info, paste0(getwd(), "/transcript.info.rds"))

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

#subset "Metastatic" $ "primary Tumor"
cer_phenodata <- cer_phenodata %>% 
  filter(sample_type.samples %in% c("Primary Tumor")) %>% as.data.frame() %>% 
  column_to_rownames(var = "submitter_id.samples")

#subset only samples that are "Metastatic" & "Primary Tumor" in cervical.exp_mat
cervical.exp_mat <- cervical.exp_mat[, colnames(cervical.exp_mat) %in% rownames(cer_phenodata)] %>% as.data.frame()

#check if desired sample type is appropriately retrieved
#table(cer_phenodata$sample_type.samples)

#View phenodata 
head(cer_phenodata)

#@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@
#retrieve survival data 
survData_link <- "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-CESC.survival.tsv"

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

##############################
#Normalization by TMM in edgeR
##############################
#Organize data in a summarized experiment (se) for downstream
cer_se <- SummarizedExperiment(assays = cervical.exp_mat,
                               rowData = transcript,
                               colData = new_cer_phenodata)

#save summarized experiment as RDS for posterity
saveRDS(cer_se, "/Users/prosperchukwuemeka/Movies/Ziglar/cer_se.rds")

#read summarized experiment
cer_se <- readRDS("cer_se.rds")
#assay(cer_se)
#colData(cer_se)
#rowData(cer_se)

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

#save dge as RDS
saveRDS(object = dge,
        file = "cer_dge.RDS",
        compress = FALSE)

#read dge
dge <- readRDS("cer_dge.RDS")

#get tmm normalized count
cervicalData_tmm_normalize <- edgeR::cpm(dge, log = FALSE) %>% as.matrix.default()

#View normalized count matrix
#View(cervicalData_tmm_normalize)

##################################################################################################
### From our exploratory analysis both edgeR and vst normalization work best so either can be used
### NB: for the sake of this work shop we will proceed with vst transformed expression matrix
##################################################################################################

###############################################
#subset gene of interests for survival analysis
###############################################
#subset TP53, MDM2, & BCL6
survival_genes <- cervicalData_tmm_normalize[rownames(cervicalData_tmm_normalize) %in% c("TPX2", "AURKA", "CCNA2", "MDM2", "CDK2"),]

#transpose survival_genes dataframe
survival_genes <- as.data.frame(t(survival_genes))

#subset overall survival (OS) and overall survival time (OS.time.year) from summarized experiment
new_survival_info <- as.data.frame(colData(cer_se)[, c("OS","OS.time.year", "initial_weight.samples")], row.names = rownames(colData(cer_se)))

#rename column name in new_survival_info to fit survminer package
colnames(new_survival_info) <- c(OS = "event", OS.time.year = "time", initial_weight.samples = "Initial_weight_of_samples")

#create the dataframe for survival analysis
survival_df <- cbind(survival_genes, new_survival_info)

#save survival_df as rds 
saveRDS(survival_df, "/Users/prosperchukwuemeka/Movies/Ziglar/survival_df.rds")

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
#fit survival curves for "TPX2" and visualize
surv.fit <- survfit(Surv(time, event) ~ TPX2, data = surv.cat)
TPX2.plot <- ggsurvplot(surv.fit, data = surv.cat, risk.table = F, conf.int = T, surv.median.line = "hv",
           pval = T, legend.title = "TPX2", legend.labs = c("High", "Low")) + 
  labs(x = "Time (Years)", y = "Survival probability")

#set plot parameters
ggpar(TPX2.plot, 
      font.main = c(16, "bold"),
      font.x = c(16, "bold"),
      font.y = c(16, "bold"),
      font.caption = c(16, "bold"), 
      font.legend = c(16, "bold"), 
      font.tickslab = c(16, "bold"))

#save survival plot 
ggsave("TPX2.plot.png", width = 15, height = 12, units = "cm")

#@@@@@@@@@@@
#@@@@@@@@@@@
#fit survival curves for "AURKA" and visualize
surv.fit <- survfit(Surv(time, event) ~AURKA, data = surv.cat)
AURKA.plot <- ggsurvplot(surv.fit, data = surv.cat, risk.table = F, conf.int = T, surv.median.line = "hv",
           pval = T, legend.title = "AURKA", legend.labs = c("High", "Low")) + 
  labs(x = "Time (Years)", y = "Survival probability")

#set plot parameters
ggpar(AURKA.plot, 
      font.main = c(16, "bold"),
      font.x = c(16, "bold"),
      font.y = c(16, "bold"),
      font.caption = c(16, "bold"), 
      font.legend = c(16, "bold"), 
      font.tickslab = c(16, "bold"))

#save survival plot 
ggsave("AURKA.plot.png", width = 15, height = 12, units = "cm")

#@@@@@@@@@@@
#@@@@@@@@@@@
#fit survival curves for "CCNA2" and visualize
surv.fit <- survfit(Surv(time, event) ~CCNA2, data = surv.cat)
CCNA2.plot <- ggsurvplot(surv.fit, data = surv.cat, risk.table = F, conf.int = T, surv.median.line = "hv",
           pval = T, legend.title = "CCNA2", legend.labs = c("High", "Low")) + 
  labs(x = "Time (Years)", y = "Survival probability")

#set plot parameters
ggpar(CCNA2.plot, 
      font.main = c(16, "bold"),
      font.x = c(16, "bold"),
      font.y = c(16, "bold"),
      font.caption = c(16, "bold"), 
      font.legend = c(16, "bold"), 
      font.tickslab = c(16, "bold"))

#save survival plot 
ggsave("CCNA2.plot.png", width = 15, height = 12, units = "cm")

#@@@@@@@@@@@
#@@@@@@@@@@@
#fit survival curves for "MDM2" and visualize
surv.fit <- survfit(Surv(time, event) ~MDM2, data = surv.cat)
MDM2.plot <- ggsurvplot(surv.fit, data = surv.cat, risk.table = F, conf.int = T, surv.median.line = "hv",
                        pval = T, legend.title = "MDM2", legend.labs = c("High", "Low")) + 
  labs(x = "Time (Years)", y = "Survival probability")

#set plot parameters
ggpar(MDM2.plot, 
      font.main = c(16, "bold"),
      font.x = c(16, "bold"),
      font.y = c(16, "bold"),
      font.caption = c(16, "bold"), 
      font.legend = c(16, "bold"), 
      font.tickslab = c(16, "bold"))

#save survival plot
ggsave("MDM2.plot.png", width = 15, height = 12, units = "cm")

#@@@@@@@@@@@
#@@@@@@@@@@@
#fit survival curves for "CDK2" and visualize
surv.fit <- survfit(Surv(time, event) ~CDK2, data = surv.cat)
CDK2.plot <- ggsurvplot(surv.fit, data = surv.cat, risk.table = F, conf.int = T, surv.median.line = "hv",
                       pval = T, legend.title = "CDK2", legend.labs = c("High", "Low")) + 
  labs(x = "Time (Years)", y = "Survival probability")

#set plot parameters
ggpar(CDK2.plot, 
      font.main = c(16, "bold"),
      font.x = c(16, "bold"),
      font.y = c(16, "bold"),
      font.caption = c(16, "bold"), 
      font.legend = c(16, "bold"), 
      font.tickslab = c(16, "bold"))

#save survival plot 
ggsave("CDK2.plot.png", width = 15, height = 12, units = "cm")
