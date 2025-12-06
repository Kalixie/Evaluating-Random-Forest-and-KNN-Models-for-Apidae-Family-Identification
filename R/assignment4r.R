### Evaluating Random Forest and KNN Models for Apidae Family Identification, Assignment 4 BINF 6210 ----
# Nadira Robertson
# 05/12/25

# https://github.com/Kalixie/Evaluating-Random-Forest-and-KNN-Models-for-Apidae-Family-Identification

## Packages Used ----

library(rentrez)
library(seqinr)
library(Biostrings)
library(tidyverse)
library(caret)
library(randomForest)
library(ggplot2)

## Data searching and collection from NCBI ----

df_Apidae_search <- entrez_search(db = "nuccore", term = "Apidae [Organism] AND COI[Gene]", use_history = T)

## Load and Process data ----

# Load FASTA

ss_apidae <- readDNAStringSet("../data/apidae.fasta")

# Create data frame with titles and sequences

df_apidae <- data.frame(title = names(ss_apidae), sequence = paste(ss_apidae))

seq_total <- nrow(df_apidae)
print(seq_total)

rm(ss_apidae)

# Create nucleotides2 column in DNAStringSet format

df_apidae$nucleotides2 <- DNAStringSet(df_apidae$sequence)

# Generate dinucleotide frequencies for analysis

df_apidae <- cbind(df_apidae, as.data.frame(dinucleotideFrequency(df_apidae$nucleotides2, as.prob = TRUE)))

# Adding trinucleotide frequency, k-mers = 3

df_apidae <- cbind(df_apidae, as.data.frame(trinucleotideFrequency(df_apidae$nucleotides2, as.prob = TRUE)))

# Remove nucleotides2 column

df_apidae <- df_apidae %>%
  select(-nucleotides2)

# Extract genus and unique ID from FASTA headers

df_apidae <- df_apidae %>%
  mutate(genus = word(title, 2L), unique_id = word(title, 1L))

unique(df_apidae$genus)

# Remove improper/unverified sequences

remove <- c("UNVERIFIED:", "Apidae", "Hymenoptera", "Apinae", "Meliponini")

df_apidae <- df_apidae %>%
  filter(!genus %in% remove)

unique(df_apidae$genus)

rm(remove)

# Calculate sequence lengths

df_apidae$seq_length <- nchar(df_apidae$sequence)

# Calculate summary of sequences

summary(df_apidae$seq_length)

# Basic plot histogram of sequence lengths

hist(df_apidae$seq_length)

# Filter outlier sequences

df_apidae <- df_apidae %>%
  filter(between(seq_length, 350, 700))

# Better histogram with ggplot

ggplot(df_apidae, aes(x = seq_length)) +
  geom_histogram(binwidth = 70, fill = "plum3", color = "black") +
  scale_x_continuous(breaks = seq(0, max(df_apidae$seq_length), by = 200)) +
  theme_minimal() +
  labs(title = "Histogram of Apidae Sequence Lengths", x = "Sequence Length (bp)", y = "Count")

summary(df_apidae$seq_length)

unique(df_apidae$genus)

## Prepare Data for Random Forest ----

# Keep only data with over 150 per genus

ls_vectors <- df_apidae %>%
  group_by(genus) %>%
  dplyr::count() %>%
  filter(n > 150) %>%
  pull(genus)

df_apidae_vector <- df_apidae %>%
  filter(genus %in% ls_vectors)

unique(df_apidae_vector$genus)

rm(ls_vectors)

# For equal sampling setting min sample size

smaller_sample <- min(table(df_apidae_vector$genus))

# Split into validation (20%) and training (80%)

set.seed(10)
df_apidae_validation <- df_apidae_vector %>%
  group_by(genus) %>%
  sample_n(floor(0.2 * smaller_sample))

set.seed(12)
df_apidae_training <- df_apidae_vector %>%
  filter(!unique_id %in% df_apidae_validation$unique_id) %>%
  group_by(genus) %>%
  sample_n(ceiling(0.8 * smaller_sample))

rm(df_apidae_vector)

# Output of training set

str(df_apidae_training)

## Random forest model on dinucleotide frequencies (k-mer = 2) ----

apidae_class <- randomForest(
  x = df_apidae_training[, 3:18],
  y = as.factor(df_apidae_training$genus),
  ntree = 500, importance = TRUE
)

# Plot results and importance

apidae_class
plot(apidae_class)
varImpPlot(apidae_class)

# Predict on validation set

predict_validation <- predict(
  apidae_class,
  df_apidae_validation[, 3:18]
)

# Compare predicted vs observed

table(
  observed = df_apidae_validation$genus,
  predicted = predict_validation
)

## Random forest model on trinucleotide frequencies (k-mer = 3) ----

# Random Forest using trinucleotide frequencies

apidae_class_tri <- randomForest(
  x = df_apidae_training[, 19:82],
  y = as.factor(df_apidae_training$genus),
  ntree = 500, importance = TRUE
)

# Plot results

apidae_class_tri
plot(apidae_class_tri)
varImpPlot(apidae_class_tri)

# Predict on validation set

predict_validation_tri <- predict(
  apidae_class_tri,
  df_apidae_validation[, 19:82]
)

# Compare predicted vs observed
table(
  observed = df_apidae_validation$genus,
  predicted = predict_validation_tri
)

## Prepare alternate machine learning method using K-Nearest Neighbors (KNN) ----

# Set safe k values (smaller than smallest class)

max_k <- smaller_sample - 1

# Make a list of potential k values

k_values <- seq(1, max_k, by = 1)

rm(smaller_sample)

# Prepare target

target_train <- as.factor(df_apidae_training$genus)
target_val <- as.factor(df_apidae_validation$genus)

## Run KNN for dinucleotide frequencies (k-mer = 2)----

# Prepare features

features_train <- as.data.frame(df_apidae_training[, 3:18])
features_val <- as.data.frame(df_apidae_validation[, 3:18])

# Train k-NN with caret

set.seed(12)
control <- trainControl(method = "cv", number = 10)

model_knn <- train(
  x = features_train,
  y = target_train,
  method = "knn",
  trControl = control,
  preProcess = c("center", "scale"),
  tuneGrid = data.frame(k = k_values)
)

#  View results for k value selection

print(model_knn)

plot(model_knn)

# Predict on validation set

pred_val <- predict(model_knn, features_val)

rm(features_train, features_val)

# Compare predicted vs observed

table(observed = target_val, predicted = pred_val)

## Run KNN for trinucleotide frequencies (k-mer = 3)----

features_train_tri <- as.data.frame(df_apidae_training[, 19:82])
features_val_tri <- as.data.frame(df_apidae_validation[, 19:82])

set.seed(12)
model_knn_tri <- train(
  x = features_train_tri,
  y = target_train,
  method = "knn",
  trControl = control,
  preProcess = c("center", "scale"),
  tuneGrid = data.frame(k = k_values)
)

#  View results for k value selection

print(model_knn_tri)

plot(model_knn_tri)

# Predict

pred_val_tri <- predict(model_knn_tri, features_val_tri)

rm(features_train_tri, features_val_tri)

# Compare predicted vs observed

table(observed = target_val, predicted = pred_val_tri)

## Compare Overall model accuracy (Barplot) ----

dfaccurancy <- data.frame(
  model_type = rep(c("RF", "KNN"), each = 2),
  kmer_type = rep(c("Dinucleotide", "Trinucleotide"), 2),
  accuracy = c(
    # Calculate mean accuracy values by comparing
    
    mean(predict_validation == df_apidae_validation$genus),
    mean(predict_validation_tri == df_apidae_validation$genus),
    mean(pred_val == target_val),
    mean(pred_val_tri == target_val)
  )
)

ggplot(dfaccurancy, aes(x = model_type, y = accuracy, fill = kmer_type)) +
  geom_col(position = "dodge") +
  ylim(0, 1) +
  theme_minimal() +
  scale_fill_manual(values = c("Dinucleotide" = "lightblue", "Trinucleotide" = "plum3")) +
  labs(title = "Overall Model Accuracy Comparison", y = "Accuracy", x = "Model", fill = "K-mer Value") +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 14))

## Confusion Matrix Generation ----

# Make factor to prevent errors

true_val <- factor(df_apidae_validation$genus)
pred_rf_di <- factor(predict_validation)
pred_rf_tri <- factor(predict_validation_tri)
pred_knn_di <- factor(pred_val)
pred_knn_tri <- factor(pred_val_tri)

# Function for creating data frames

make_cm_df <- function(tr, pr) {
  as.data.frame(table(True = tr, Pred = pr))
}

cm_rf_di_df <- make_cm_df(true_val, pred_rf_di)
cm_rf_tri_df <- make_cm_df(true_val, pred_rf_tri)
cm_knn_di_df <- make_cm_df(true_val, pred_knn_di)
cm_knn_tri_df <- make_cm_df(true_val, pred_knn_tri)

# Plot with ggplot function

# Define colour palettes

cols_di  <- c("lightblue", "steelblue3")
cols_tri <- c("plum", "mediumpurple1")


plot_matrix <- function(true, pred, title, colours) {
  cm <- as.data.frame(table(True = true, Pred = pred))
  
  ggplot(cm, aes(Pred, True, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = Freq)) +
    scale_fill_gradient(low = colours[1], high = colours[2]) +
    labs(title = title, x = "Predicted genus", y = "Observed genus") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Use plot function

plot_matrix(true_val, pred_rf_di,  "RF — Dinucleotide",  cols_di)
plot_matrix(true_val, pred_rf_tri, "RF — Trinucleotide", cols_tri)
plot_matrix(true_val, pred_knn_di, "KNN — Dinucleotide", cols_di)
plot_matrix(true_val, pred_knn_tri,"KNN — Trinucleotide",cols_tri)


