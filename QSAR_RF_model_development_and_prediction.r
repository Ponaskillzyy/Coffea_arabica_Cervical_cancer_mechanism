########
#load libraries
########
suppressWarnings(suppressPackageStartupMessages({library(dplyr); library(ggExtra); library(pROC);
  library(tidyr); library(ggplotify); library(missForest); library(visdat); library(caret); library(corrplot);
  library(ggrastr); library(stringr); library(enrichplot); library(data.table); library(ggpubr);
  library(tibble); library(stats); library(randomForest)}))

#setwd 
setwd("/Users/prosperchukwuemeka/Movies/Ziglar/New_CC_arabica/Revision/Descriptor_Prediction")

#load data
mmp7_inh_qsar_imputed_data <- readRDS("mmp7_inh_qsar_imputed_data.rds")

#check data for missingness
vis_miss(mmp7_inh_qsar_imputed_data, warn_large_data = FALSE)

############
#Data preprocessing
############
#remove near-zero variance predictors
nzv <- nearZeroVar(mmp7_inh_qsar_imputed_data, saveMetrics = TRUE)
rf_data <- mmp7_inh_qsar_imputed_data[, !nzv$nzv]

#split into features (X) and target (y)
X <- rf_data[, -which(names(rf_data) == "pIC50")]  # Features
y <- rf_data$pIC50  # Target variable

#normalize data (Scaling)
preProcValues <- preProcess(X, method = c("center", "scale"))
X <- predict(preProcValues, X)

#combine processed features with target variable
rf_data_preprocessed <- data.frame(X, pIC50 = y)

################
#split data into training (80%) and testing (20%)
################
#set seed
set.seed(123)

#split data
trainIndex <- createDataPartition(rf_data_preprocessed$pIC50, p = 0.8, list = FALSE)
trainData <- rf_data_preprocessed[trainIndex, ]
testData <- rf_data_preprocessed[-trainIndex, ]

#################
#Feature selection using Recursive Feature Elimination (RFE)
#################
#set seed
set.seed(123)

#set rfe hyperparameter 
control <- rfeControl(functions = rfFuncs, method = "cv", number = 10)

#run rfe model
rfeModel <- rfe(trainData[, -ncol(trainData)], trainData$pIC50, sizes = c(5, 10, 15, 20),
                rfeControl = control)

# #saveRDS
# saveRDS(rfeModel, "rfeModel.rds")
# 
# #load rfe model
# rfeModel <- readRDS("rfeModel.rds")

#selected features
selected_features <- predictors(rfeModel)
trainData <- trainData[, c(selected_features, "pIC50")]
testData <- testData[, c(selected_features, "pIC50")]

################
#Train Random Forest model with Cross-Validation
################
#set seed
set.seed(123)

#define the fit control
fitControl <- trainControl(method = "cv", number = 10)  # 10-fold cross-validation

#train model
rf_model <- train(pIC50 ~ ., data = trainData, method = "rf",
                  trControl = fitControl, importance = TRUE)

#train set prediction
train_pred <- predict(rf_model, trainData)

#model performance Metrics on training set
train_r2 <- cor(train_pred, trainData$pIC50)^2
train_rmse <- sqrt(mean((trainData$pIC50 - train_pred)^2))
cat("R² Train: ", train_r2, "\nRMSE Test:", train_rmse, "\n")

#########
#predict on test set
#########
#test set prediction
test_pred <- predict(rf_model, newdata = testData)

#model performance Metrics on test set
test_r2 <- cor(testData$pIC50, test_pred)^2
test_rmse <- sqrt(mean((testData$pIC50 - test_pred)^2))
cat("R² Test:", test_r2, "\nRMSE Test:", test_rmse, "\n")

########
#create a dataframe for ggplot
########
train_df <- data.frame(Observed = trainData$pIC50, Predicted = train_pred, Set = "Training")
test_df <- data.frame(Observed = testData$pIC50, Predicted = test_pred, Set = "Test")

#combine both dataframes
plot_data <- rbind(train_df, test_df)

#re-level set
plot_data$Set <- factor(plot_data$Set, levels = c("Training", "Test"))

############
#plot observed vs predicted
############
#create annotation text
annotation_text <- paste0(
  "Training: R² = ", round(train_r2, 3), ", RMSE = ", round(train_rmse, 3), "; ",
  "Test: R² = ", round(r2_test, 3), ", RMSE = ", round(rmse_test, 3)
)

#themes
themes <- theme(strip.text = element_text(face = "bold", size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1),
                legend.text = element_text(face = "bold", size = 15),legend.title = element_text(face = "bold", size = 15),
                plot.title = element_text(face = "bold", size = 15), axis.title = element_text(face = "bold", size = 15),
                axis.text.x = element_text(face = "bold", size = 15, colour = "black", angle = 52, hjust = 1), 
                axis.text.y = element_text(face = "bold", size = 15, colour = "black"), legend.position = "none",
                panel.background = element_rect(fill = "white"), axis.line = element_line(linewidth = 0.8, colour = "black"))

#plot
ggplot(plot_data, aes(x = Observed, y = Predicted, color = Set)) +
  geom_point(aes(color = Set), alpha = 0.6, size = 4) +  # Scatter points
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", size = 1.2) +  # Regression line
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "solid", size = 1.2) +  # Perfect fit
  labs(title = "Random Forest Regression Model (Observed vs Predicted pIC50)", x = "Observed pIC50", y = "Predicted pIC50") +
  theme_minimal() + themes +
  theme(legend.position = "top") +
  scale_color_manual(values = c("Training" = "blue", "Test" = "red")) +
  coord_cartesian(xlim = c(0, max(plot_data$Observed)), ylim = c(0, max(plot_data$Predicted))) +
  annotate("text", x = max(plot_data$Observed) * 1, y = min(plot_data$Predicted), 
           label = annotation_text, hjust = 1, size = 5, color = "black") 

#############
#View import features contributing to the model performance
#############
#Extract importance and standard deviation
importance_df <- as.data.frame(importance(rf_model$finalModel, type = 1, scale = TRUE))  # type=1 gives %IncMSE
importance_df$Feature <- rownames(importance_df)

#compute standard deviation of %IncMSE
sd_incMSE <- sd(importance_df$`%IncMSE`)

#compute standard error
n <- nrow(importance_df)
SE <- sd_incMSE / sqrt(n)

#compute 95% confidence interval
importance_df <- importance_df %>%
  mutate(
    Lower_CI = `%IncMSE` - 1.96 * SE,
    Upper_CI = `%IncMSE` + 1.96 * SE
  )

#restructure data into a long format
importance_df <- importance_df %>%
  pivot_longer(cols = c(`%IncMSE`, `Lower_CI`, `Upper_CI`), 
               names_to = "Importance Metric", values_to = "Value")

#rename item
importance_df$Feature <- ifelse(as.character(importance_df$Feature) == "MDEC.33", "MDEC-33", as.character(importance_df$Feature))

#reorder the Feature factor by the median of Value, in descending order
importance_df$Feature <- factor(importance_df$Feature, 
                                levels = importance_df %>%
                                  group_by(Feature) %>%
                                  summarise(median_value = median(Value)) %>%
                                  arrange(desc(median_value)) %>%
                                  pull(Feature))

#print results
print(importance_df)

#Plot %IncMSE
importance_df %>%
  ggplot(aes(x = Feature, y = Value)) +  # Use Feature directly here instead of reorder
  geom_boxplot(outlier.shape = NA, fill = "grey80", show.legend = FALSE) + 
  geom_point(aes(color = `Importance Metric`), size = 4, alpha = 1, position = position_jitter(width = 0.03, height = 0)) + 
  labs(
    title = "Feature Importance Plot", 
    x = "Features", 
    y = "% Increase in MSE", 
    color = "Importance Metric"
  ) + 
  theme_bw() + 
  themes + 
  theme(legend.position = "top") + 
  #coord_flip() + 
  stat_compare_means(
    aes(group = Feature),  # Group by Feature
    method = "t.test",  # Use the t-test for comparisons
    label = "p.signif",  # Show significance levels
    ref.group = "ATSC4v",  # Use "ATSC4v" as the reference group
    size = 6,  # Set text size for better visibility
    fontface = "bold",  # Make the significance text bold
    hide.ns = T
  )

##########
#Extract the 10-Fold CV Results and plot
##########
#add an ID to each fold to track them
cv_long <- rf_model$resample %>%
  mutate(Fold = row_number()) %>%
  pivot_longer(cols = c(RMSE, Rsquared, MAE), names_to = "Metric", values_to = "Value")

# Paired line plot
ggplot(cv_long, aes(x = Metric, y = Value, group = Fold)) +
  geom_line(aes(color = as.factor(Fold)), linewidth = 1.2, alpha = 0.6) +  # Connect RMSE and R² for each fold
  geom_point(aes(color = as.factor(Fold)), size = 7) +  # Show points
  scale_color_viridis_d(option = "B") +  # Optional: make it look nice
  theme_bw() +
  labs(title = " Dynamics of Random Forest Regression Model Performance across Folds (CV = 10)",
       x = "Metric",
       y = "Value",
       color = "Fold") + 
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  ) + themes


##########
#Model descriptors orthogonality test
##########
#get model trainData
orth_data <- rf_model$trainingData %>%
  dplyr::select(-.outcome)

#rename column MDEC.33 to MDEC-33
names(orth_data) <- ifelse(as.character(names(orth_data)) == "MDEC.33", 
                           "MDEC-33", as.character(names(orth_data)))

#get descriptor order from %IncMSE plot and relevel descriptor orthogonality data
orth_data <- orth_data %>% .[, print(levels(importance_df$Feature))] 

#plot correlation between descriptors
corrplot.mixed(cor(orth_data, method = "pearson"), 
               lower = "square", 
               upper = "number",
               tl.pos = "lt",
               diag = "l",
               tl.col = "black",
               tl.cex = 1.4) 

# #add plot title
# title("Descriptor Correlation Matrix", line = 1, cex.main = 1.5, font.main = 2)

#@@@@@@@@@@@@
#############
#Load c_arabica docked compounds imputed PaDel Descriptors
#############
#@@@@@@@@@@@@
c_arabica_docked_compounds_desc_imputed_data <- readRDS("docked_compounds_desc_imputed_data.rds")

#check data for missingness
vis_miss(c_arabica_docked_compounds_desc_imputed_data, warn_large_data = FALSE)

#################
#filter dataset to contain only features used during the QSAR model development
#################
#use selected features from rfe model
c_arabica_docked_compounds_desc_imputed_data <- c_arabica_docked_compounds_desc_imputed_data %>% 
  data.frame() %>%
  .[, names(.) %in% selected_features]

#########
#normalize data (Scaling) using same approach as in model development
#########
c_arabica_preProcValues <- preProcess(c_arabica_docked_compounds_desc_imputed_data, 
                                      method = c("center", "scale"))
c_arabica_data_preprocessed <- predict(c_arabica_preProcValues, 
                                       c_arabica_docked_compounds_desc_imputed_data)

#predict c_arabica compound bioactivity 
c_arabica_data_preprocessed$pIC50 <- predict(rf_model, c_arabica_data_preprocessed)

#view data
View(c_arabica_data_preprocessed)

#rename item
names(c_arabica_data_preprocessed) <- ifelse(as.character(names(c_arabica_data_preprocessed)) == "MDEC.33", 
                                                              "MDEC-33", as.character(names(c_arabica_data_preprocessed)))

# #saveRSD
# saveRDS(c_arabica_data_preprocessed, "c_arabica_docked_compound_pIC50_pred.rds")
# 
# #write csv file
# write.csv(c_arabica_data_preprocessed, "c_arabica_docked_compound_pIC50_pred.csv")




