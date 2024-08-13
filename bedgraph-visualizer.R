# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)  # Alternatively, you could use patchwork
library(scales)   # For formatting axis labels

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) < 4) {
  stop("Please provide the mode ('genome_plot' or 'region_plot'), followed by the required arguments: <Mode> <BAF file> <LRR file> <Region BED file (only for region_plot)> <Output Directory>")
}

# Assign command-line arguments to variables
mode <- args[1]
baf_file <- args[2]
lrr_file <- args[3]
output_dir <- args[length(args)]

# Ensure the output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Load the BAF and LRR data using the provided file paths
baf_data <- read.table(gzfile(baf_file), header = TRUE)
lrr_data <- read.table(gzfile(lrr_file), header = TRUE)

# Ensure columns are named correctly
colnames(baf_data) <- c("Chr", "Start", "End", "BAF")
colnames(lrr_data) <- c("Chr", "Start", "End", "LRR")

# Function to plot the entire genome
genome_plot <- function() {
  # Unique chromosomes
  chromosomes <- unique(lrr_data$Chr)

  # Create a list to hold plots
  lrr_plots <- list()
  baf_plots <- list()

  # Generate plots for each chromosome
  for (chr in chromosomes) {
    # Filter data for the current chromosome
    lrr_chr_data <- lrr_data %>% filter(Chr == chr)
    baf_chr_data <- baf_data %>% filter(Chr == chr)

    # Create LRR plot
    lrr_plot <- ggplot(lrr_chr_data) +
      geom_point(aes(x = Start, y = LRR), color = "#2a2a2a") +
      geom_smooth(aes(x = Start, y = LRR), se = FALSE, method = "lm", color = "#ff3333") +  # Add smooth line
      theme_minimal() +
      labs(x = NULL, y = "LRR", title = paste("chr", chr, sep="")) +
      geom_hline(yintercept = 0, linetype = "dotted", color = "#555555") +  # Guide line for LRR at y = 0
      coord_cartesian(ylim = c(-1, 1)) +  # Set the y-axis limits for LRR
      scale_x_continuous(labels = label_number()) +  # Format x-axis labels as numeric
      scale_y_continuous(expand = c(0, 0))  # Avoid extra space at the plot edges

    # Create BAF plot
    baf_plot <- ggplot(baf_chr_data) +
      geom_point(aes(x = Start, y = BAF), color = "#2a2a2a") +
      theme_minimal() +
      labs(x = "Position", y = "BAF") +
      geom_hline(yintercept = 0.5, linetype = "dotted", color = "#555555") +  # Guide line for BAF at y = 0.5
      coord_cartesian(ylim = c(0, 1)) +
      scale_x_continuous(labels = label_number()) +  # Format x-axis labels as numeric
      scale_y_continuous(expand = c(0, 0))  # Avoid extra space at the plot edges

    # Store plots in lists
    lrr_plots[[chr]] <- lrr_plot
    baf_plots[[chr]] <- baf_plot
  }

  # Combine plots for each chromosome
  all_plots <- lapply(chromosomes, function(chr) {
    plot_grid(lrr_plots[[chr]], baf_plots[[chr]], align = "v", ncol = 1)
  })
  # Split all plots into chunks of 12
  plot_chunks <- split(all_plots, ceiling(seq_along(all_plots) / 12))

  # Combine all chunks into a single plot
  combined_plot <- plot_grid(plotlist = unlist(plot_chunks, recursive = FALSE), ncol = 4)

  output_file <- file.path(output_dir, "plot_genome.png")
  ggsave(output_file, plot = combined_plot, width = 30, height = 40, dpi = 300)
}

# Function to plot specific regions
region_plot <- function(region_file) {
  # Load region data
  regions <- read.table(region_file, header = FALSE)
  colnames(regions) <- c("Chr", "RegionStart", "RegionEnd", "Type")

  # Iterate over each region and create plots
  for (i in 1:nrow(regions)) {
    region <- regions[i, ]

    # Filter BAF and LRR data for the current region
    filtered_baf <- baf_data %>%
      filter(Chr == region$Chr & Start >= region$RegionStart & End <= region$RegionEnd)

    filtered_lrr <- lrr_data %>%
      filter(Chr == region$Chr & Start >= region$RegionStart & End <= region$RegionEnd)

    # Create the LRR plot with smooth line
    lrr_plot <- ggplot(filtered_lrr) +
      geom_point(aes(x = Start, y = LRR), color = "#2a2a2a") +
      geom_smooth(aes(x = Start, y = LRR), se = FALSE, method = "lm", color = "#ff3333") +  # Add smooth line
      theme_minimal() +
      labs(x = NULL, y = "LRR", title = paste(region$Chr, ":", region$RegionStart, "-", region$RegionEnd, sep="")) +
      geom_hline(yintercept = 0, linetype = "dotted", color = "#555555") +  # Guide line for LRR at y = 0
      coord_cartesian(ylim = c(-1, 1)) +  # Set the y-axis limits for LRR
      scale_x_continuous(labels = label_number())  # Format x-axis labels as numeric

    # Create the BAF plot
    baf_plot <- ggplot(filtered_baf) +
      geom_point(aes(x = Start, y = BAF), color = "#2a2a2a") +
      theme_minimal() +
      labs(x = "Position", y = "BAF") +
      geom_hline(yintercept = 0.5, linetype = "dotted", color = "#555555") +  # Guide line for BAF at y = 0.5
      coord_cartesian(ylim = c(0, 1)) +
      scale_x_continuous(labels = label_number())  # Format x-axis labels as numeric

    # Combine the plots vertically
    combined_plot <- plot_grid(lrr_plot, baf_plot, align = "v", ncol = 1, rel_heights = c(1, 1))

    # Create a filename for the current region
    output_file <- file.path(output_dir, paste0("plot_", region$Chr, "_", region$RegionStart, "_", region$RegionEnd, ".png"))

    # Save the combined plot to a PNG file
    ggsave(output_file, plot = combined_plot, width = 10, height = 8, dpi = 300)
  }
}

# Execute the appropriate function based on the mode
if (mode == "genome_plot") {
  genome_plot()
} else if (mode == "region_plot") {
  region_file <- args[4]
  region_plot(region_file)
} else {
  stop("Invalid mode. Use 'genome_plot' or 'region_plot'.")
}
