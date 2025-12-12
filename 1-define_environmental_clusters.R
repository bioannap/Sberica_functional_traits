#############
# K-means clustering analysis and PCA
############

# Load required libraries
library(readxl)
library(cluster)   # to compute silhouette scores
library(corrplot)     # to visualize correlation matrices
library(factoextra)
library(tidyverse)
library(tidyplots)

setwd("C:/Users/Saxifraga_berica/Statistics/")
# Import dataset from Excel file
Saxydata_per_analisi <- read_excel("Saxydata_per_analisi_1.xlsx", 
                                   col_types = c("text", "text", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric",
                                                 "numeric"))
View(Saxydata_per_analisi)

# Convert tibble to data.frame
Saxydata_per_analisi <- as.data.frame(Saxydata_per_analisi)

# Use first column as rownames (IDs), then drop it from the dataset
rownames(Saxydata_per_analisi) <- Saxydata_per_analisi[[1]]
Saxydata_per_analisi <- Saxydata_per_analisi[,-1]
#select only environmental data
env_data <- Saxydata_per_analisi[, 2:16]
env_data<-na.omit(env_data)
# Copy to df (scaled version for k-means input)
df <- scale(env_data)
# Also keep a non-scaled version to append clusters later
data <- env_data

# Compute correlation matrix between variables
corrmatrix <- cor(env_data)

# Visualize correlation matrix with correlation coefficients as numbers
corrplot(corrmatrix, method = 'number')

#use uncorrelated variables
# Select relevant columns for clustering (columns 1–4 and 30–41 after removing col 1)
vars_to_use <- c("Light_Intensity(PAR)",
                 "Ambient_Humidity",
                 "Ambient_Temperature",
                 "PPFD(400-700nm)",
                 "R:FR",
                 "R:B")
env_data <- env_data[, vars_to_use]
# Copy to df (scaled version for k-means input)
df <- scale(env_data)
# Also keep a non-scaled version to append clusters later
data <- env_data
# Compute correlation matrix between variables
corrmatrix <- cor(env_data)
# Visualize correlation matrix with correlation coefficients as numbers
corrplot(corrmatrix, method = 'number')

# Define function to compute average silhouette score for a given number of clusters (k)
silhouette_score <- function(k){
  km <- kmeans(df, centers = k, nstart = 25)        # run K-means with 'k' clusters
  ss <- silhouette(km$cluster, dist(df))            # compute silhouette values
  mean(ss[, 3])                                     # return average silhouette width
}

# Test different cluster numbers (from 2 to 10)
k <- 2:10
avg_sil <- sapply(k, silhouette_score)

#plot optimal number of clusters
fviz_nbclust(df, kmeans, method='silhouette')


# Perform final clustering with chosen number of clusters (e.g., k=2)
km.final <- kmeans(df, 2)

# Print total within-cluster sum of squares (measure of compactness)
km.final$tot.withinss
# Print cluster sizes (number of elements per cluster)
km.final$size

# Add cluster assignments as a new column in the dataset
data$cluster <- km.final$cluster
data$site_tag <- Saxydata_per_analisi$site_tag[match(rownames(data), rownames(Saxydata_per_analisi))]

# Display first 6 rows of dataset with new 'cluster' column
head(data, 6)


#PCA
# Run PCA on the scaled data
pca <- prcomp(df, center = TRUE, scale. = TRUE)

# Extract scores (PC coordinates of samples)
scores <- data.frame(pca$x[,1:2])
scores$cluster <- factor(km.final$cluster)# add cluster info
scores$site_tag <- factor(data$site_tag)
scores$sample <- rownames(df)

#SEE CONTRIBUTION OF EACH VARIABLE
# Extract loadings (directions of variables)
loadings <- data.frame(pca$rotation[,1:2])
loadings$variable <- rownames(loadings)

# Scale arrows to fit on the same plot as scores
arrow_scale <- 5  # adjust as needed for visibility
loadings_scaled <- loadings
loadings_scaled$PC1 <- loadings_scaled$PC1 * arrow_scale
loadings_scaled$PC2 <- loadings_scaled$PC2 * arrow_scale

#decide cluster colors for points (bright/saturated)
point_colors <- c("#0072B2","#F0E442", "#56B4E9", "#E69F00",  "#009E73",   "#D55E00", "#CC79A7")

# Create muted/pastel versions for density contours
density_colors <- c("#B3D9F2", "#F7F1A1", "#ABCEF4", "#F2CF80", "#80CEB9", "#EAAF80", "#E6BCD3")

##### IMPROVED KERNEL DENSITY PLOT
# Biplot with better distinction between datapoints and densities
ggplot(scores, aes(x = PC1, y = PC2)) +
  # Add density contours per cluster with muted colors and lower alpha
  stat_density_2d_filled(geom = "polygon",
                         contour_var = "ndensity", 
                         breaks = c(0.2, 0.5, 0.8, 0.95, 1),
                         aes(x = PC1, y = PC2, fill = cluster),
                         alpha = 0.3,  # Lower alpha for subtlety
                         inherit.aes = FALSE) +
  
  # Add density contour lines for better definition
  stat_density_2d(aes(x = PC1, y = PC2, color = cluster), 
                  contour_var = "ndensity",
                  breaks = c(0.2, 0.5, 0.8, 0.95),
                  alpha = 0.6,
                  linewidth = 0.5,
                  inherit.aes = FALSE) +
  
  # Add points with bright colors, white outline for pop
  geom_point(aes(x = PC1, y = PC2, color = cluster, shape = site_tag), 
             size = 2.5, 
             stroke = 1, 
             fill = "white") +  # White fill for hollow shapes
  
  # Add arrows for loadings
  geom_segment(data = loadings_scaled,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               inherit.aes = FALSE,
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "black", 
               linewidth = 0.8) +
  
  # Add variable labels with white background for readability
  geom_text(data = loadings_scaled,
            aes(x = PC1, y = PC2, label = variable),
            inherit.aes = F,
            size = 3, color = "black") +
  
  labs(title = "PCA Biplot with K-means clusters and Kernel Density") +
  
  # Use bright colors for points
  scale_color_manual(values = point_colors, name = "Cluster") +
  # Use muted colors for density fills
  scale_fill_manual(values = density_colors, name = "Cluster") +
  scale_shape_manual(values = c(1, 2, 3, 4, 5, 6, 7, 8), name = "Site") +
  
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) +
  coord_equal()

#SAVE PCA COORDINATES
pca_coord<-scores[,c(1,2,4)]
write.table(pca_coord, "C:/Users/annap/Documents/Ricerca/Plant_genomics/Saxifraga_berica/Statistics/PCA_coord_2clusters_env.tsv", 
            sep="\t", row.names=FALSE, quote=FALSE)
metadata_2clusters<-scores[,c(3,5)]
write.table(metadata_2clusters, "C:/Users/annap/Documents/Ricerca/Plant_genomics/Saxifraga_berica/Statistics/metadata_2clusters_env.tsv", 
            sep="\t", row.names=FALSE, quote=FALSE)

####CHECK LIGHT/DARK

#Tidy the dataset
dataset_boxplot<- data %>%
  pivot_longer( names_to = "Env_var", values_to = "value", cols = c("Light_Intensity(PAR)",
                                                                    "Ambient_Humidity",
                                                                    "Ambient_Temperature",
                                                                    "PPFD(400-700nm)",
                                                                    "R:FR",
                                                                    "R:B")) %>%
  arrange(value, Env_var)
View(dataset_boxplot)
dataset_boxplot <- dataset_boxplot %>%
  mutate(
    cluster = as.factor(cluster),
  )
#PLOT with wilocox test
color <- c("#0072B2","#F0E442")
dataset_boxplot |>
  tidyplot(x = cluster, y = value, color = cluster) |>
  add_boxplot(alpha = 0.4) |>
  add_test_asterisks(method = "wilcox_test", hide_info = TRUE) |>
  adjust_colors(new_colors = color,saturation = 1) |>
  split_plot(by = Env_var)


## Get pairwise Wilcoxon results
library(rstatix)
# Initialize results list
results <- list()

# Loop through each trait variable
for (var in env_vars) {
  test <- wilcox.test(data[[var]] ~ data$cluster)
  results[[var]] <- c(
    p_value = test$p.value,
    statistic = test$statistic
  )
}

# Convert to dataframe
results_table <- as.data.frame(do.call(rbind, results))
results_table$variable <- rownames(results_table)
rownames(results_table) <- NULL

# Reorder columns
results_table <- results_table[, c("variable", "p_value", "statistic.W")]

# View/save results
print(results_table)
write.csv(results_table, "pairwise_wilcox_results_environment.csv", row.names = FALSE)
