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
  library(stringr); library(edgeR); library(ggrepel);
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
#metadata to look into "disease_type (primary tumor, no metastatic disease)", "age_at_index.demographic (20 to 88 years)", "treatment_or_therapy.treatments.diagnoses (therapy naive)"
#########
#subset "primary Tumor"
cer_phenodata <- cer_phenodata %>% 
  filter(sample_type.samples %in% c("Primary Tumor")) %>% as.data.frame() %>% 
  column_to_rownames(var = "sample")

#subset only samples that are "Primary Tumor" in cervical.exp_mat
cervical.exp_mat <- cervical.exp_mat[, colnames(cervical.exp_mat) %in% rownames(cer_phenodata)] %>% as.data.frame()


###########
#Load GTEx normal cervix tissue expression data
###########
GTEx.normal.tissue.exp_mat.ecto <- fread(paste0("/Users/prosperchukwuemeka/Movies/Ziglar/New_CC_arabica/Cervical-cancer-healthy-tissue-ExpressionData", "/gene_reads_cervix_ectocervix.gct")) %>% dplyr::select(-Name) 
GTEx.normal.tissue.exp_mat.endo <- fread(paste0("/Users/prosperchukwuemeka/Movies/Ziglar/New_CC_arabica/Cervical-cancer-healthy-tissue-ExpressionData", "/gene_reads_cervix_endocervix.gct")) %>% dplyr::select(-Name) 

#check if description in both dataset matches
all(GTEx.normal.tissue.exp_mat.ecto$Description %in% GTEx.normal.tissue.exp_mat.endo$Description)

#cbind both GTEx dataset
GTEx.cervix.normal.tissue.exp_mat <- cbind(GTEx.normal.tissue.exp_mat.ecto, GTEx.normal.tissue.exp_mat.endo) %>% as.data.frame() 

#remove duplicate columns
names(GTEx.cervix.normal.tissue.exp_mat)[1] <- "Gene"
GTEx.cervix.normal.tissue.exp_mat <- GTEx.cervix.normal.tissue.exp_mat %>% dplyr::select(-Description)

#select only genes found in cervical cancer expression Data
GTEx.cervix.normal.tissue.exp_mat <- GTEx.cervix.normal.tissue.exp_mat %>% .[.$Gene %in% rownames(cervical.exp_mat),]

#check for duplicate rows based on gene column
any(duplicated(GTEx.cervix.normal.tissue.exp_mat$Gene))

#remove duplicates while retaining first occurrence
GTEx.cervix.normal.tissue.exp_mat <- GTEx.cervix.normal.tissue.exp_mat %>% .[!duplicated(.$Gene), ] 

#rename rownames with gene symbol
rownames(GTEx.cervix.normal.tissue.exp_mat) <- GTEx.cervix.normal.tissue.exp_mat$Gene
GTEx.cervix.normal.tissue.exp_mat <- GTEx.cervix.normal.tissue.exp_mat %>% dplyr::select(-Gene) %>% as.data.frame()

#filter out genes not present in healthy cervix tissue from cervical cancer expression data 
cervical.exp_mat <- cervical.exp_mat %>% .[rownames(.) %in% rownames(GTEx.cervix.normal.tissue.exp_mat),]

#reorder the gene names of the healthy tissue expression data by the cervical cancer gene order
GTEx.cervix.normal.tissue.exp_mat <- GTEx.cervix.normal.tissue.exp_mat[match(rownames(cervical.exp_mat), rownames(GTEx.cervix.normal.tissue.exp_mat)), ]

###########
#make phenotype data for GTEx normal cervix tissue expression data
###########
GTEx.cervix.phenoData <- data.frame("sample" = colnames(GTEx.cervix.normal.tissue.exp_mat), row.names = colnames(GTEx.cervix.normal.tissue.exp_mat))
GTEx.cervix.phenoData$sample_type.samples <- "Normal Cervix Tissue"

###########
#make new phenotype data concatenating TCGA tumor data and GTEx normal tissue data
###########
phenoData <-  bind_rows(cer_phenodata, GTEx.cervix.phenoData)

###########
#make new Expression data concatenating TCGA tumor data and GTEx normal tissue data
###########
ExpressionData <- bind_cols(cervical.exp_mat, GTEx.cervix.normal.tissue.exp_mat)

#check dimension of "ExpressionData" & "phenoData"; note that n_samples must be the same in all dataframe
dim(ExpressionData)
dim(phenoData)

#reorder samples in phenoData to match ExpressionData
phenoData <- phenoData %>% .[match(colnames(ExpressionData), rownames(.)), ]

##############################
#Normalization by TMM in edgeR
##############################
#Organize data in a summarized experiment (se) for downstream
cer_se <- SummarizedExperiment(assays = ExpressionData,
                               rowData = rownames(ExpressionData),
                               colData = phenoData)

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
#order gene expression
res <- res %>% 
  group_by(Expression) %>%
  arrange(desc(logFC)) %>% column_to_rownames("Gene")
res$Gene <- rownames(res)

# #save DEG result
# saveRDS(res, "cc_deg_result.rds")

#readRDS
res <- readRDS("cc_deg_result.rds")

#filter upregulated genes
upreg.gene <- res %>% filter(Expression == "Up-regulated in Tumor")

#write csv
write_csv(upreg.gene, "upreg.gene.csv")

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

# #Apply a small offset to p-value
# offset <- 0 + 1e-50 #prevent having a pvalue of exactly zero since log10(0) is undefined
# 
# #index genes with padj with exactly 0
# ndx.0.padj <- which(res$FDR==0)
# res[ndx.0.padj,] <- offset

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

#@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@
#load compound-disease shared targets
compound_disease_targets <- fread("Common Targets.csv") %>% pull(`Common Targets`)

#ndx for stratifying patients based on tissue biopsy
normal_tissue_idents <- rownames(dge$samples)[which(dge$samples$group == "Normal Cervix Tissue")]
primary_tumor_idents <- rownames(dge$samples)[which(dge$samples$group == "Primary Tumor")]

#get expression matrix
expression_data <- cervicalData_tmm_normalize

#subset expression data for gene of interest
expression_data <- as.data.frame(expression_data[rownames(cervicalData_tmm_normalize) %in% compound_disease_targets,]) %>% t()

#get tissue details from phenotype data
phenotype_data <- dge$samples %>% dplyr::select(group)
names(phenotype_data) <- "Tissue"

#Combine the dataframes by column
expression_pData <- bind_cols(expression_data, phenotype_data)

#convert data to long format
expression_pData <- expression_pData %>% rownames_to_column(var = "Sample") %>%
  pivot_longer(names_to = "Gene", values_to = "Expression", 2:47)

#####
#plot graph with statistical significance
#####
#themes
themes <- theme(strip.text = element_text(face = "bold", size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1),
                legend.text = element_text(face = "bold", size = 15),legend.title = element_text(face = "bold", size = 15),
                plot.title = element_text(face = "bold", size = 15), axis.title = element_text(face = "bold", size = 15),
                axis.text.x = element_text(face = "bold", size = 15, colour = "black", angle = 40, hjust = 1), 
                axis.text.y = element_text(face = "bold", size = 15, colour = "black"),
                axis.ticks.x = element_blank(), legend.position = "None",
                panel.background = element_rect(fill = "white"), axis.line = element_line(linewidth = 0.8, colour = "black"))

p <- expression_pData %>%
  ggplot(., aes(x=Tissue, y= Expression)) +
  geom_boxplot(outlier.shape = NA, fill = "grey80", show.legend = FALSE) +
  geom_point(aes(color = Tissue), size = 4, alpha = 1, position = position_jitter(width = 0.03, height = 0)) +
  labs(title = "Expression of Shared Compound-Cervical Cancer Targets in Normal (n=47) vs. Primary Tumor (n=304)", 
       x = "Gene", 
       y = "Expression (cpm)", 
       color = "Tissue") + 
  theme_bw() + 
  themes +
  theme(legend.position = "right") + # Keep the legend at the top
  facet_wrap(~Gene, scales = "free", nrow = 7) 

#Add statistical analysis
p + stat_compare_means(method = "wilcox.test", aes(group = Tissue), 
                       fontface = "bold", label = "p.signif", hide.ns = T, 
                       size = 8, vjust = 0.7)













