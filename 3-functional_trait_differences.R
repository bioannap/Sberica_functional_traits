# Load required packages
library(tidyverse)
library(readxl)
library(vegan)
library(ggfortify)
library(broom)
library(tidyplots)
library(corrplot)

setwd("C:/Users/Saxifraga_berica/Statistics/trait-env_correlations")
# -------------------------------
# 1. Import and merge data
# -------------------------------
metadata <- read_tsv("metadata_2clusters_env.tsv",col_types = "cc")
traits_data <- read_excel("Saxydata_per_analisi_1.xlsx", 
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

df <- traits_data %>%
  inner_join(metadata, by = "sample") %>%
  mutate(cluster = factor(cluster, levels = c("1", "2")))

df <- df %>%
  rename(
    Light_Intensity_PAR = `Light_Intensity(PAR)`,
    PPFD_400_700nm      = `PPFD(400-700nm)`,
    R_FR                = `R:FR`,
    R_B                 = `R:B`
  )
# -------------------------------
# 2. Define variable groups
# -------------------------------
trait_vars <- c("LEF","NPQt","Phi2","PhiNO","PhiNPQ","FmPrime",
                "FoPrime","Fs","FvP/FmP","qL","Leaf_temperature_differential",
                "Chl_a[mg/g]","Chl_b[mg/g]","Car[mg/g]","Chl_tot[mg/g]", "Chl_a+b/Car",
                "Chl_a/Chl b","LDMC","SLA")

env_vars <- c("Light_Intensity_PAR",
              "Ambient_Humidity",
              "Ambient_Temperature",
              "PPFD_400_700nm",
              "R_FR",
              "R_B")

df<-df %>% drop_na("LEF","NPQt","Phi2","PhiNO","PhiNPQ","FmPrime",
               "FoPrime","Fs","FvP/FmP","qL","Leaf_temperature_differential",
               "Chl_a[mg/g]","Chl_b[mg/g]","Car[mg/g]","Chl_tot[mg/g]","Chl_a+b/Car",
               "Chl_a/Chl b","LDMC","SLA","Light_Intensity_PAR",
               "Ambient_Humidity",
               "Ambient_Temperature",
               "PPFD_400_700nm",
               "R_FR",
               "R_B")

df_traits <- df %>% select(all_of(trait_vars))
df_env    <- df %>% select(cluster, all_of(env_vars))

#check for autocorrelating traits
corrmatrix <- cor(df_traits)
# Visualize correlation matrix with correlation coefficients as numbers
corrplot(corrmatrix, method = 'number')

# -------------------------------
# 3. Global multivariate test
# -------------------------------
# PERMANOVA with categorical predictor
permanova1 <- adonis2(df_traits ~ cluster, 
                     data = df, method = "euclidean")
print(permanova1)

# Test for homogeneity of group dispersions (important for PERMANOVA validity)
disper <- betadisper(vegdist(df_traits, method = "euclidean"), df_env$cluster)
anova(disper)
#test is not significant, the assumption holds and the PERMANOVA result is valid

# -------------------------------
# 4. Ordination (PCA)
# -------------------------------
# PCA of traits
pca_res <- prcomp(df_traits, scale. = TRUE)

# Extract scores (PC coordinates of samples)
scores <- data.frame(pca_res$x[,1:2])
scores$cluster <- factor(df$cluster)# add cluster info
scores$site_tag <- factor(df$site_tag)
scores$sample <- df$sample
# Extract loadings (directions of variables)
loadings <- data.frame(pca_res$rotation[,1:2])
loadings$variable <- rownames(loadings)

# Scale arrows to fit on the same plot as scores
arrow_scale <- 5  # adjust as needed for visibility
loadings_scaled <- loadings
loadings_scaled$PC1 <- loadings_scaled$PC1 * arrow_scale
loadings_scaled$PC2 <- loadings_scaled$PC2 * arrow_scale
#calculate total variance explained by each principal component
pca_res$sdev^2 / sum(pca_res$sdev^2)

#decide cluster colors
color <- c("#0072B2","#F0E442", "#56B4E9", "#E69F00",  "#009E73",   "#D55E00", "#CC79A7")
##### ADD KERNEL DENSITY PLOT
# Biplot with kernel density contours
ggplot(scores, aes(x = PC1, y = PC2, color = cluster)) +
  # add density contours per cluster
  stat_density_2d_filled(geom = "polygon",
                         contour_var = "ndensity", 
                         breaks = c(0.2, 0.5, 0.8, 0.95, 1),
                         aes(x = PC1, y = PC2, fill= cluster, alpha = after_stat(level)),
                         inherit.aes = F) +
  # add points
  geom_point(size = 2
             ) +
  # add arrows for loadings
  #geom_segment(data = loadings_scaled,
  #             aes(x = 0, y = 0, xend = PC1*1.2, yend = PC2*1.2),
  #             inherit.aes = F,
  #             color = "black") +
  # add variable labels at arrow tips
  geom_text(data = loadings_scaled,
            aes(x = PC1*1.5, y = PC2*1.5, label = variable),
            inherit.aes = F,
            size = 3, color = "black") +
  labs(title = "PCA Biplot with K-means clusters and Kernel Density") +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) +
  theme_minimal() +
  coord_equal()

#SAVE PCA COORDINATES
pca_coord<-scores[,c(1,2,4)]
write.table(pca_coord, "C:/Users/Saxifraga_berica/Statistics/trait-env_correlations/pca_coord_traits.tsv", 
            sep="\t", row.names=FALSE, quote=FALSE)
plot(rda_model, scaling = 2)

# -------------------------------
# 5. Quantitative differences
# -------------------------------
## Get pairwise Wilcoxon results
library(rstatix)
# Initialize results list
results <- list()

# Loop through each trait variable
for (var in trait_vars) {
  test <- wilcox.test(df[[var]] ~ df$cluster)
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
write.csv(results_table, "pairwise_wilcox_results_functional-traits.csv", row.names = FALSE)

#Tidy the dataset
dataset_boxplot <- df %>%
  pivot_longer( names_to = "Character", values_to = "value", cols = trait_vars) %>%
  arrange(value, Character)
View(dataset_boxplot)

#PLOT with wilocox test
color <- c("#0072B2","#F0E442")
dataset_boxplot |> 
  tidyplot(x = cluster, y = value, color = cluster) |> 
  add_boxplot(alpha = 0.4) |> 
  add_test_asterisks(method = "wilcox_test", hide_info = TRUE) |> 
  adjust_colors(new_colors = color,saturation = 1) |> 
  add_data_points_beeswarm() |>
  split_plot(by = Character)

#------------------------------
# 6. Trait-by-trait regressions (GLM)
# -------------------------------

trait_tests <- lapply(trait_vars, function(tr) {
  m <- glm(df[[tr]] ~ Light_Intensity_PAR + Ambient_Humidity + Ambient_Temperature+ PPFD_400_700nm + R_FR + R_B,
           data = df,
           family = gaussian())  # GLM with Gaussian family
  
  tidy <- broom::tidy(m)
  data.frame(Trait = tr, tidy)
}) %>% bind_rows()

# Extract effects for Environment + each covariate
trait_effects <- trait_tests %>%
  filter(term != "(Intercept)") %>%
  group_by(term) %>%
  mutate(p_adj = p.adjust(p.value, method = "fdr")) %>%
  ungroup()

# Top hits per predictor
Trait_predictors<- trait_effects %>%
  filter(p_adj < 0.05) %>%
  arrange(p_adj)

write.table(Trait_predictors, "C:/Users/Plant_genomics/Saxifraga_berica/Statistics/trait-env_correlations/Trait_predictors.tsv", 
            sep="\t", row.names=FALSE, quote=FALSE)

#------------------------------
# 7. Seed differences
# -------------------------------
## Get pairwise Wilcoxon results
# Initialize results list
seed_vars <- c("G%", "MGT")
results <- list()

# Loop through each trait variable
for (var in seed_vars) {
  test <- wilcox.test(df[[var]] ~ df$cluster)
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
write.csv(results_table, "pairwise_wilcox_results_germination-traits.csv", row.names = FALSE)

#Tidy the dataset
dataset_seeds<- df %>%
  pivot_longer( names_to = "Character", values_to = "value", cols = c("G%", "MGT")) %>%
  arrange(value, Character)
View(dataset_seeds)

#PLOT with wilocox test
color <- c("#0072B2","#F0E442")
dataset_seeds |> 
  tidyplot(x = cluster, y = value, color = cluster) |> 
  add_boxplot(alpha = 0.4) |> 
  add_test_asterisks(method = "t_test", hide_info = TRUE) |> 
  adjust_colors(new_colors = color,saturation = 1) |> 
  add_data_points_beeswarm() |>
  split_plot(by = Character)

#------------------------------
# 8. Verify homogeneous seed weight and size
# ------------------------------
#Weight
anova_weight <- aov(`Average_weight_of_50_seeds[g]` ~ factor(cluster), data = df)
summary(anova)
#Size
anova_size <- aov(seed_axis_ratio ~ factor(cluster), data = df)
summary(anova)

seeds<- read_excel("dataset_semi.xlsx", 
                                 col_types = c("text", "text", "numeric", 
                                               "numeric"))
anova_seeds <- aov(`Weight of 50 seeds  [g]` ~ Individual, data = seeds)
summary(anova_seeds)
kw_seeds <- kruskal.test(`Weight of 50 seeds  [g]` ~ Individual, data = seeds)

