########
#load libraries
########
suppressWarnings(suppressPackageStartupMessages({library(dplyr); library(ggExtra);
  library(tidyr); library(ggplotify); library(missForest); library(visdat); 
  library(ggrastr); library(stringr); library(enrichplot); library(data.table); 
  library(tibble); library(stats)}))

#setwd 
setwd("/Users/prosperchukwuemeka/Movies/Ziglar/New_CC_arabica/Revision/Descriptor_Prediction")

#@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@
#load dataset
BindingDB_QSAR <- fread("BindingDB_MMP7_inhibitors_for_QSAR.csv") %>%
  as.data.frame()

#filter dataset remove compounds with IC50 outliers 
BindingDB_QSAR <- BindingDB_QSAR %>% 
  filter(!str_detect(`IC50 (nM)`, ">"))

#save data
saveRDS(BindingDB_QSAR, "BindingDB_QSAR.rds")

#select compound name and smile strings
BindingDB_QSAR.padel <- BindingDB_QSAR %>% 
  dplyr::select(`BindingDB Reactant_set_id`, `Ligand SMILES`)

# #write csv file
# write_csv(BindingDB_QSAR.padel, "BindingDB_QSAR.padel.csv")

###########
#Load mmp7 inh PaDel Descriptors
###########
#load
mmp7_inh_descriptors <- fread("MMP7_inh_decriptor_results.csv") 

#add inh IC50 to dataframe
mmp7_inh_qsar_data <- merge(mmp7_inh_descriptors, BindingDB_QSAR %>% 
                              dplyr::select(`BindingDB Reactant_set_id`, `IC50 (nM)`),
                            by = "BindingDB Reactant_set_id", sort = FALSE)

#convert IC50 to pIC50
mmp7_inh_qsar_data$pIC50 <- -log10(as.numeric(mmp7_inh_qsar_data$`IC50 (nM)`) * 1e-9)

#remove IC50 and smiles column
mmp7_inh_qsar_data <- mmp7_inh_qsar_data %>% 
  dplyr::select(-`IC50 (nM)`, - `Ligand SMILES`)

##############
#impute data to fill missing data
##############
#set rownames
mmp7_inh_qsar_data <- mmp7_inh_qsar_data %>% 
  column_to_rownames("BindingDB Reactant_set_id") 

#check data for missingness
vis_miss(mmp7_inh_qsar_data, warn_large_data = FALSE)

#use random forest model to handle missingness
imputed_data <- missForest(mmp7_inh_qsar_data)$ximp

#check data for missingness
vis_miss(imputed_data, warn_large_data = FALSE)

# #saveRDS
# saveRDS(imputed_data, "mmp7_inh_qsar_imputed_data.rds")
# 
# #write csv file
# write.csv(imputed_data_qsar, "mmp7_inh_qsar_imputed_data.csv")

###########
#Load c_arabica docked compounds PaDel Descriptors
###########
#load
c_arabica_docked_compounds_desc <- fread("Coffea_arabica_Docked_Compounds_desc.csv") %>%
  dplyr::select(-SMILES) %>%
  column_to_rownames("Compound") %>%
  as.data.frame()

##############
#impute data to fill missing data
##############
#check data for missingness
vis_miss(c_arabica_docked_compounds_desc, warn_large_data = FALSE)

#use random forest model to handle missingness
c_arabica_docked_compounds_desc_imputed_data <- missForest(c_arabica_docked_compounds_desc)$ximp

#check data for missingness
vis_miss(c_arabica_docked_compounds_desc_imputed_data, warn_large_data = FALSE)

# #saveRDS
# saveRDS(c_arabica_docked_compounds_desc_imputed_data, "docked_compounds_desc_imputed_data.rds")
# 
# #write csv file
# write.csv(c_arabica_docked_compounds_desc_imputed_data, "docked_compounds_desc_imputed_data.csv")
# 




