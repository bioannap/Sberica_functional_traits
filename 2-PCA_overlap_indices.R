# PCA Contour Overlap Analysis (Adapted for Environmental Data)


# ==========================
# 1️ Setup
# ==========================
setwd("C:/Users/Saxifraga_berica/Statistics/k-means_2clusters")

# Load required libraries
required_packages <- c("sf", "ks", "MASS", "ggplot2", "dplyr", 
                       "grDevices", "alphahull", "geometry")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Function to prepare PCA data
prepare_pca_data <- function(pca_data_file,
                             metadata_file,
                             coord_cols = c("PC1", "PC2"),
                             cluster_col = "cluster",
                             sample_col = "sample") {
  # Read PCA data safely
  pca_data <- read.table(pca_data_file,
                         header = TRUE,
                         sep = "\t",
                         stringsAsFactors = FALSE,
                         comment.char = "",   # keep # in IDs
                         quote = "")
  
  # Check coordinate columns exist
  missing_coords <- setdiff(coord_cols, colnames(pca_data))
  if (length(missing_coords) > 0) {
    stop("Missing coordinate columns in PCA file: ",
         paste(missing_coords, collapse = ", "))
  }
  
  # Keep only relevant columns
  pca_data <- pca_data[, c(sample_col, coord_cols), drop = FALSE]
  
  # Read metadata safely
  metadata <- read.table(metadata_file,
                         header = TRUE,
                         sep = "\t",
                         stringsAsFactors = FALSE,
                         comment.char = "",
                         quote = "")
  
  # Convert sample and cluster columns to character to ensure proper handling
  metadata[[sample_col]] <- as.character(metadata[[sample_col]])
  metadata[[cluster_col]] <- as.character(metadata[[cluster_col]])
  pca_data[[sample_col]] <- as.character(pca_data[[sample_col]])
  
  # Merge by sample ID
  merged <- merge(pca_data, metadata,
                  by.x = sample_col,
                  by.y = sample_col)
  
  # Check cluster column exists
  if (!(cluster_col %in% colnames(merged))) {
    stop("Cluster column '", cluster_col, "' not found in metadata.")
  }
  
  return(merged)
}

# Function 1: Convex Hull Overlap Analysis (updated for sf)
calculate_convex_hull_overlap <- function(coords1, coords2, group1_name, group2_name) {
  
  cat(paste("Calculating convex hull overlap:", group1_name, "vs", group2_name, "\n"))
  
  tryCatch({
    # Create convex hulls
    hull1_indices <- chull(coords1)
    hull2_indices <- chull(coords2)
    
    hull1_coords <- coords1[hull1_indices, ]
    hull2_coords <- coords2[hull2_indices, ]
    
    # Close the polygons by adding first point at the end
    hull1_coords <- rbind(hull1_coords, hull1_coords[1, ])
    hull2_coords <- rbind(hull2_coords, hull2_coords[1, ])
    
    # Convert to sf polygons
    hull1_matrix <- as.matrix(hull1_coords)
    hull2_matrix <- as.matrix(hull2_coords)
    
    # Create sf polygons
    hull1_poly <- st_polygon(list(hull1_matrix))
    hull2_poly <- st_polygon(list(hull2_matrix))
    
    hull1_sf <- st_sfc(hull1_poly)
    hull2_sf <- st_sfc(hull2_poly)
    
    # Calculate areas
    area1 <- as.numeric(st_area(hull1_sf))
    area2 <- as.numeric(st_area(hull2_sf))
    
    # Calculate intersection
    intersection <- st_intersection(hull1_sf, hull2_sf)
    
    # Handle case where there's no intersection
    if (length(intersection) == 0) {
      overlap_area <- 0
    } else {
      overlap_area <- as.numeric(st_area(intersection))
    }
    
    # Calculate metrics
    union_area <- area1 + area2 - overlap_area
    jaccard_index <- overlap_area / union_area
    overlap_percent_group1 <- overlap_area / area1
    overlap_percent_group2 <- overlap_area / area2
    
    results <- list(
      method = "convex_hull",
      group1 = group1_name,
      group2 = group2_name,
      area1 = area1,
      area2 = area2,
      overlap_area = overlap_area,
      union_area = union_area,
      jaccard_index = jaccard_index,
      overlap_percent_group1 = overlap_percent_group1,
      overlap_percent_group2 = overlap_percent_group2,
      hull1_coords = hull1_coords[1:(nrow(hull1_coords)-1), ],  # Remove duplicate last point for plotting
      hull2_coords = hull2_coords[1:(nrow(hull2_coords)-1), ]
    )
    
    return(results)
    
  }, error = function(e) {
    cat("Error in convex hull calculation:", e$message, "\n")
    cat("Trying fallback method...\n")
    
    # Fallback: simple bounding box overlap
    tryCatch({
      bbox1 <- apply(coords1, 2, range)
      bbox2 <- apply(coords2, 2, range)
      
      # Calculate bounding box areas
      area1 <- (bbox1[2,1] - bbox1[1,1]) * (bbox1[2,2] - bbox1[1,2])
      area2 <- (bbox2[2,1] - bbox2[1,1]) * (bbox2[2,2] - bbox2[1,2])
      
      # Calculate overlap
      x_overlap <- max(0, min(bbox1[2,1], bbox2[2,1]) - max(bbox1[1,1], bbox2[1,1]))
      y_overlap <- max(0, min(bbox1[2,2], bbox2[2,2]) - max(bbox1[1,2], bbox2[1,2]))
      overlap_area <- x_overlap * y_overlap
      
      union_area <- area1 + area2 - overlap_area
      jaccard_index <- overlap_area / union_area
      
      results <- list(
        method = "bounding_box_fallback",
        group1 = group1_name,
        group2 = group2_name,
        area1 = area1,
        area2 = area2,
        overlap_area = overlap_area,
        union_area = union_area,
        jaccard_index = jaccard_index,
        overlap_percent_group1 = overlap_area / area1,
        overlap_percent_group2 = overlap_area / area2,
        hull1_coords = NULL,
        hull2_coords = NULL
      )
      
      return(results)
      
    }, error = function(e2) {
      cat("Both convex hull methods failed:", e2$message, "\n")
      return(NULL)
    })
  })
}

# Function 2: Kernel Density Overlap Analysis
calculate_kernel_density_overlap <- function(coords1, coords2, group1_name, group2_name, 
                                             grid_size = 100, bandwidth = NULL) {
  
  cat(paste("Calculating kernel density overlap:", group1_name, "vs", group2_name, "\n"))
  
  tryCatch({
    # Determine bandwidth if not provided
    if (is.null(bandwidth)) {
      bandwidth1 <- Hpi(coords1)
      bandwidth2 <- Hpi(coords2)
      bandwidth <- (bandwidth1 + bandwidth2) / 2
    }
    
    # Create common grid
    x_range <- range(c(coords1[,1], coords2[,1]))
    y_range <- range(c(coords1[,2], coords2[,2]))
    
    x_seq <- seq(x_range[1], x_range[2], length.out = grid_size)
    y_seq <- seq(y_range[1], y_range[2], length.out = grid_size)
    grid_points <- expand.grid(x = x_seq, y = y_seq)
    
    # Calculate kernel densities
    kde1 <- kde(coords1, H = bandwidth, eval.points = grid_points)
    kde2 <- kde(coords2, H = bandwidth, eval.points = grid_points)
    
    # Normalize densities to sum to 1
    kde1_norm <- kde1$estimate / sum(kde1$estimate)
    kde2_norm <- kde2$estimate / sum(kde2$estimate)
    
    # Calculate overlap metrics
    # Bhattacharyya coefficient (measure of overlap between probability distributions)
    bhattacharyya_coeff <- sum(sqrt(kde1_norm * kde2_norm))
    
    # Hellinger distance (1 - Bhattacharyya coefficient)
    hellinger_distance <- sqrt(1 - bhattacharyya_coeff)
    
    # Overlap area (intersection of normalized densities)
    overlap_density <- pmin(kde1_norm, kde2_norm)
    overlap_area <- sum(overlap_density)
    
    results <- list(
      method = "kernel_density",
      group1 = group1_name,
      group2 = group2_name,
      bhattacharyya_coefficient = bhattacharyya_coeff,
      hellinger_distance = hellinger_distance,
      overlap_area = overlap_area,
      kde1 = kde1,
      kde2 = kde2,
      grid_points = grid_points,
      bandwidth = bandwidth
    )
    
    return(results)
    
  }, error = function(e) {
    cat("Error in kernel density calculation:", e$message, "\n")
    
    # Fallback to simpler 2D density estimation
    tryCatch({
      cat("Trying alternative density estimation...\n")
      
      # Use MASS::kde2d as fallback
      x_range <- range(c(coords1[,1], coords2[,1]))
      y_range <- range(c(coords1[,2], coords2[,2]))
      
      kde1_alt <- MASS::kde2d(coords1[,1], coords1[,2], 
                              lims = c(x_range, y_range), n = grid_size)
      kde2_alt <- MASS::kde2d(coords2[,1], coords2[,2], 
                              lims = c(x_range, y_range), n = grid_size)
      
      # Normalize
      kde1_norm <- kde1_alt$z / sum(kde1_alt$z)
      kde2_norm <- kde2_alt$z / sum(kde2_alt$z)
      
      # Calculate overlap
      bhattacharyya_coeff <- sum(sqrt(kde1_norm * kde2_norm))
      hellinger_distance <- sqrt(1 - bhattacharyya_coeff)
      overlap_area <- sum(pmin(kde1_norm, kde2_norm))
      
      results <- list(
        method = "kernel_density_fallback",
        group1 = group1_name,
        group2 = group2_name,
        bhattacharyya_coefficient = bhattacharyya_coeff,
        hellinger_distance = hellinger_distance,
        overlap_area = overlap_area,
        kde1_alt = kde1_alt,
        kde2_alt = kde2_alt
      )
      
      return(results)
      
    }, error = function(e2) {
      cat("Both density methods failed:", e2$message, "\n")
      return(NULL)
    })
  })
}

# Function 3: Distance-based Separation Analysis
calculate_distance_separation <- function(coords1, coords2, group1_name, group2_name) {
  
  cat(paste("Calculating distance-based separation:", group1_name, "vs", group2_name, "\n"))
  
  # Calculate centroids
  centroid1 <- colMeans(coords1)
  centroid2 <- colMeans(coords2)
  centroid_distance <- sqrt(sum((centroid1 - centroid2)^2))
  
  # Calculate within-group distances
  dist1 <- as.matrix(dist(coords1))
  dist2 <- as.matrix(dist(coords2))
  
  # Remove diagonal (distance to self = 0)
  dist1_values <- dist1[upper.tri(dist1)]
  dist2_values <- dist2[upper.tri(dist2)]
  
  mean_within_dist1 <- mean(dist1_values)
  mean_within_dist2 <- mean(dist2_values)
  mean_within_dist_combined <- mean(c(dist1_values, dist2_values))
  
  # Calculate between-group distances
  between_distances <- as.vector(as.matrix(dist(rbind(coords1, coords2)))[1:nrow(coords1), (nrow(coords1)+1):(nrow(coords1)+nrow(coords2))])
  mean_between_dist <- mean(between_distances)
  
  # Separation metrics
  separation_index <- centroid_distance / mean_within_dist_combined
  hopkins_statistic <- mean_between_dist / mean_within_dist_combined
  
  # Silhouette-like index
  silhouette_index <- (mean_between_dist - mean_within_dist_combined) / 
    max(mean_between_dist, mean_within_dist_combined)
  
  results <- list(
    method = "distance_separation",
    group1 = group1_name,
    group2 = group2_name,
    centroid_distance = centroid_distance,
    mean_within_dist1 = mean_within_dist1,
    mean_within_dist2 = mean_within_dist2,
    mean_within_dist_combined = mean_within_dist_combined,
    mean_between_dist = mean_between_dist,
    separation_index = separation_index,
    hopkins_statistic = hopkins_statistic,
    silhouette_index = silhouette_index,
    centroid1 = centroid1,
    centroid2 = centroid2
  )
  
  return(results)
}

# Function 4: Comprehensive Overlap Analysis (updated for PCA data)
perform_comprehensive_overlap_analysis <- function(pca_data, comparisons = NULL) {
  
  cat("=== COMPREHENSIVE PCA OVERLAP ANALYSIS ===\n\n")
  
  # If no comparisons specified, do all pairwise comparisons
  if (is.null(comparisons)) {
    clusters <- unique(pca_data$cluster)
    comparisons <- list()
    for (i in 1:(length(clusters)-1)) {
      for (j in (i+1):length(clusters)) {
        comparisons[[paste("Cluster", clusters[i], "vs Cluster", clusters[j])]] <- c(clusters[i], clusters[j])
      }
    }
  }
  
  all_results <- list()
  
  for (comp_name in names(comparisons)) {
    cat(paste("Analyzing:", comp_name, "\n"))
    
    group1_name <- comparisons[[comp_name]][1]
    group2_name <- comparisons[[comp_name]][2]
    
    # Extract coordinates
    coords1 <- pca_data[pca_data$cluster == group1_name, c("PC1", "PC2")]
    coords2 <- pca_data[pca_data$cluster == group2_name, c("PC1", "PC2")]
    
    if (nrow(coords1) < 3 || nrow(coords2) < 3) {
      cat(paste("Skipping", comp_name, "- insufficient data points\n"))
      next
    }
    
    # Perform all three analyses
    hull_results <- calculate_convex_hull_overlap(coords1, coords2, group1_name, group2_name)
    density_results <- calculate_kernel_density_overlap(coords1, coords2, group1_name, group2_name)
    distance_results <- calculate_distance_separation(coords1, coords2, group1_name, group2_name)
    
    all_results[[comp_name]] <- list(
      convex_hull = hull_results,
      kernel_density = density_results,
      distance_separation = distance_results
    )
    
    cat("\n")
  }
  
  return(all_results)
}

# Function 5: Generate Summary Report (updated for clusters)
generate_overlap_summary <- function(overlap_results, output_file = "pca_overlap_results.txt") {
  
  cat("Generating overlap analysis summary...\n")
  
  sink(output_file)
  
  cat("=== PCA FUNCTIONAL SPACE OVERLAP ANALYSIS ===\n")
  cat("Analysis Date:", format(Sys.Date(), "%Y-%m-%d"), "\n\n")
  
  for (comp_name in names(overlap_results)) {
    cat(paste("===", toupper(comp_name), "===\n"))
    
    results <- overlap_results[[comp_name]]
    
    # Convex Hull Results
    if (!is.null(results$convex_hull)) {
      cat("\nCONVEX HULL OVERLAP:\n")
      cat(sprintf("- Jaccard Index: %.3f (0 = no overlap, 1 = complete overlap)\n", 
                  results$convex_hull$jaccard_index))
      cat(sprintf("- Cluster %s area overlap: %.1f%%\n", 
                  results$convex_hull$group1, 
                  results$convex_hull$overlap_percent_group1 * 100))
      cat(sprintf("- Cluster %s area overlap: %.1f%%\n", 
                  results$convex_hull$group2, 
                  results$convex_hull$overlap_percent_group2 * 100))
    }
    
    # Kernel Density Results
    if (!is.null(results$kernel_density)) {
      cat("\nKERNEL DENSITY OVERLAP:\n")
      cat(sprintf("- Bhattacharyya Coefficient: %.3f (0 = no overlap, 1 = identical)\n", 
                  results$kernel_density$bhattacharyya_coefficient))
      cat(sprintf("- Hellinger Distance: %.3f (0 = identical, 1 = no overlap)\n", 
                  results$kernel_density$hellinger_distance))
      cat(sprintf("- Density Overlap Area: %.3f\n", 
                  results$kernel_density$overlap_area))
    }
    
    # Distance Separation Results
    if (!is.null(results$distance_separation)) {
      cat("\nDISTANCE-BASED SEPARATION:\n")
      cat(sprintf("- Separation Index: %.3f (>1 = well separated)\n", 
                  results$distance_separation$separation_index))
      cat(sprintf("- Hopkins Statistic: %.3f (>1 = well separated)\n", 
                  results$distance_separation$hopkins_statistic))
      cat(sprintf("- Silhouette Index: %.3f (-1 to 1, >0 = well separated)\n", 
                  results$distance_separation$silhouette_index))
      cat(sprintf("- Centroid Distance: %.3f\n", 
                  results$distance_separation$centroid_distance))
    }
    
    # Overall interpretation
    cat("\nINTERPRETATION:\n")
    
    if (!is.null(results$convex_hull)) {
      jaccard <- results$convex_hull$jaccard_index
      if (jaccard < 0.1) {
        cat("- Convex hull: Minimal overlap - functionally distinct clusters\n")
      } else if (jaccard < 0.3) {
        cat("- Convex hull: Low overlap - moderately distinct clusters\n")
      } else if (jaccard < 0.6) {
        cat("- Convex hull: Moderate overlap - some functional similarity\n")
      } else {
        cat("- Convex hull: High overlap - functionally similar clusters\n")
      }
    }
    
    if (!is.null(results$distance_separation)) {
      sep_index <- results$distance_separation$separation_index
      if (sep_index > 2) {
        cat("- Distance separation: Well separated clusters\n")
      } else if (sep_index > 1) {
        cat("- Distance separation: Moderately separated clusters\n")
      } else {
        cat("- Distance separation: Overlapping clusters\n")
      }
    }
    
    cat("\n", paste(rep("=", 60), collapse = ""), "\n\n")
  }
  
  cat("INTERPRETATION GUIDELINES:\n")
  cat("- Jaccard Index: Proportion of shared functional space\n")
  cat("- Bhattacharyya Coefficient: Similarity of probability distributions\n")
  cat("- Separation Index: Centroid distance relative to within-group spread\n")
  cat("- Higher separation = more functionally distinct clusters\n")
  cat("- Lower overlap = more specialized environmental profiles\n")
  
  sink()
  
  cat(paste("Overlap analysis summary saved to:", output_file, "\n"))
}

# Function 6: Create Overlap Visualizations (updated for PCA data)
create_overlap_visualizations <- function(pca_data, overlap_results, output_prefix = "pca_overlap") {
  
  cat("Creating overlap visualizations...\n")
  
  plots <- list()
  
  for (comp_name in names(overlap_results)) {
    
    cat(paste("Creating plot for:", comp_name, "\n"))
    
    group1_name <- overlap_results[[comp_name]]$convex_hull$group1
    group2_name <- overlap_results[[comp_name]]$convex_hull$group2
    
    # Extract data for this comparison
    plot_data <- pca_data[pca_data$cluster %in% c(group1_name, group2_name), ]
    
    # Base plot
    p <- ggplot(plot_data, aes(x = PC1, y = PC2, color = cluster, fill = cluster)) +
      geom_point(size = 2, alpha = 0.7) +
      theme_minimal() +
      labs(title = paste("Environmental Space Overlap:", comp_name),
           subtitle = "PCA Analysis") +
      theme(legend.position = "bottom")
    
    # Add convex hulls if available
    if (!is.null(overlap_results[[comp_name]]$convex_hull)) {
      hull1_coords <- overlap_results[[comp_name]]$convex_hull$hull1_coords
      hull2_coords <- overlap_results[[comp_name]]$convex_hull$hull2_coords
      
      # Add hulls to plot
      hull1_df <- data.frame(hull1_coords, cluster = group1_name)
      hull2_df <- data.frame(hull2_coords, cluster = group2_name)
      
      p <- p + 
        geom_polygon(data = hull1_df, alpha = 0.2) +
        geom_polygon(data = hull2_df, alpha = 0.2)
    }
    
    # Add density contours if available
    if (!is.null(overlap_results[[comp_name]]$kernel_density)) {
      p <- p + stat_density_2d(alpha = 0.3, bins = 5)
    }
    
    plots[[comp_name]] <- p
    
    # Save individual plot
    safe_comp_name <- gsub("[^A-Za-z0-9_]", "_", comp_name)
    ggsave(paste0(output_prefix, "_", safe_comp_name, ".png"), p, 
           width = 8, height = 6, dpi = 300)
  }
  
  return(plots)
}

# Main function to run pairwise overlap analysis (updated for environmental clusters)
run_pairwise_pca_overlap_analysis <- function(pca_data_file, metadata_file, 
                                              cluster1_name, cluster2_name,
                                              coord_cols = c("PC1", "PC2"),
                                              cluster_col = "cluster",
                                              sample_col = "sample",
                                              output_prefix = NULL) {
  
  cat("=== PAIRWISE PCA OVERLAP ANALYSIS ===\n")
  cat(paste("Comparing: Cluster", cluster1_name, "vs Cluster", cluster2_name, "\n\n"))
  
  # Set default output prefix if not provided
  if (is.null(output_prefix)) {
    output_prefix <- paste0("pca_overlap_cluster", cluster1_name, "_vs_cluster", cluster2_name)
  }
  
  # Prepare data
  full_pca_data <- prepare_pca_data(pca_data_file, metadata_file, coord_cols, 
                                    cluster_col, sample_col)
  
  # Filter to only the two clusters of interest
  pca_data <- full_pca_data[full_pca_data$cluster %in% c(cluster1_name, cluster2_name), ]
  
  if (nrow(pca_data) == 0) {
    stop(paste("No data found for clusters:", cluster1_name, "and", cluster2_name))
  }
  
  cat("Data summary for this comparison:\n")
  print(table(pca_data$cluster))
  cat("\n")
  
  # Set up single comparison
  comparisons <- list()
  comp_name <- paste("Cluster", cluster1_name, "vs Cluster", cluster2_name)
  comparisons[[comp_name]] <- c(cluster1_name, cluster2_name)
  
  # Perform overlap analysis
  overlap_results <- perform_comprehensive_overlap_analysis(pca_data, comparisons)
  
  # Generate summary report
  generate_overlap_summary(overlap_results, paste0(output_prefix, "_results.txt"))
  
  # Create visualizations
  plots <- create_overlap_visualizations(pca_data, overlap_results, output_prefix)
  
  cat("\nPairwise overlap analysis completed successfully!\n")
  cat("Files created:\n")
  cat("- Results:", paste0(output_prefix, "_results.txt"), "\n")
  cat("- Plot:", paste0(output_prefix, "_", gsub("[^A-Za-z0-9_]", "_", comp_name), ".png"), "\n")
  
  return(list(
    pca_data = pca_data,
    overlap_results = overlap_results,
    plots = plots,
    comparison = comp_name
  ))
}

# Function to generate combined summary of all comparisons (updated for clusters)
generate_combined_summary <- function(all_results, output_file = "combined_pca_overlap_summary.txt") {
  
  cat("Generating combined summary of all pairwise comparisons...\n")
  
  sink(output_file)
  
  cat("=== COMBINED PCA OVERLAP ANALYSIS SUMMARY ===\n")
  cat("Analysis Date:", format(Sys.Date(), "%Y-%m-%d"), "\n")
  cat("Total Comparisons:", length(all_results), "\n\n")
  
  # Summary table
  cat("OVERLAP SUMMARY TABLE:\n")
  cat(sprintf("%-25s %-15s %-15s %-15s %-15s\n", 
              "Comparison", "Jaccard_Index", "Bhattacharyya", "Separation", "Interpretation"))
  cat(paste(rep("-", 85), collapse = ""), "\n")
  
  for (comp_name in names(all_results)) {
    result <- all_results[[comp_name]]
    
    if (!is.null(result$overlap_results)) {
      overlap_data <- result$overlap_results[[comp_name]]
      
      jaccard <- ifelse(!is.null(overlap_data$convex_hull), 
                        sprintf("%.3f", overlap_data$convex_hull$jaccard_index), "N/A")
      
      bhatta <- ifelse(!is.null(overlap_data$kernel_density), 
                       sprintf("%.3f", overlap_data$kernel_density$bhattacharyya_coefficient), "N/A")
      
      separation <- ifelse(!is.null(overlap_data$distance_separation), 
                           sprintf("%.3f", overlap_data$distance_separation$separation_index), "N/A")
      
      # Determine interpretation
      if (!is.null(overlap_data$convex_hull)) {
        j_val <- overlap_data$convex_hull$jaccard_index
        if (j_val < 0.2) interpretation <- "Distinct"
        else if (j_val < 0.5) interpretation <- "Moderate"
        else interpretation <- "Similar"
      } else {
        interpretation <- "Unknown"
      }
      
      cat(sprintf("%-25s %-15s %-15s %-15s %-15s\n", 
                  comp_name, jaccard, bhatta, separation, interpretation))
    }
  }
  
  cat("\n", paste(rep("=", 85), collapse = ""), "\n\n")
  
  cat("BIOLOGICAL/ENVIRONMENTAL INTERPRETATION:\n")
  cat("Based on environmental space overlap patterns:\n\n")
  
  for (comp_name in names(all_results)) {
    result <- all_results[[comp_name]]
    if (!is.null(result$overlap_results)) {
      overlap_data <- result$overlap_results[[comp_name]]
      
      cat(paste("•", toupper(comp_name), ":\n"))
      
      if (!is.null(overlap_data$convex_hull)) {
        jaccard <- overlap_data$convex_hull$jaccard_index
        if (jaccard < 0.2) {
          cat("  - Environmentally DISTINCT: Low overlap suggests different environmental niches\n")
        } else if (jaccard < 0.5) {
          cat("  - Environmentally INTERMEDIATE: Moderate overlap with some niche specialization\n")
        } else {
          cat("  - Environmentally SIMILAR: High overlap suggests similar environmental preferences\n")
        }
      }
      
      if (!is.null(overlap_data$distance_separation)) {
        sep_index <- overlap_data$distance_separation$separation_index
        if (sep_index > 2) {
          cat("  - Well-separated environmental clusters\n")
        } else if (sep_index > 1) {
          cat("  - Moderately separated environmental profiles\n") 
        } else {
          cat("  - Overlapping environmental space\n")
        }
      }
      cat("\n")
    }
  }
  
  sink()
  
  cat(paste("Combined summary saved to:", output_file, "\n"))
}

# ==========================
# 2️ USAGE
# ==========================

# Usage for environmental data:
# To run the analysis between cluster 1 and cluster 2:

result <- run_pairwise_pca_overlap_analysis(
   pca_data_file = "PCA_coord_2clusters_env.tsv",
   metadata_file = "metadata_2clusters_env.tsv", 
   cluster1_name = "1",
   cluster2_name = "2",
   coord_cols = c("PC1", "PC2"),
   cluster_col = "cluster",
   sample_col = "sample",
   output_prefix = "environmental_clusters_overlap"
 )

