#############################################
#### Taxonomic Binning Data Visualization
#############################################

# ============= USER INPUTS (MODIFY THESE) =============
# File paths
tax_binning_data_path <- "./4.Taxonomic-Binning/Taxonomic-Binning/5.Genera.all.tab"  # Path to taxonomic binning data
sample_data_path <- "./0.Original-Data/Template-Data/mapping_file.tab"  # Path to sample metadata

# ID Column in sample metadata/mapping file
id_col <- "#SampleID"

# Plot customization
top_n_bins <- 10  # Number of top taxonomic bins to show individually in the barplot. Currently up to 15 supported
facet_variables <- c("Diet", "Facility")  # Variables to use for faceting (grouping)
# Base width per sample
width_per_sample <- 1  
# Base width for legend and margins
base_width <- 5
# Cap of plot width
max_plot_width <- 100
# Barplot output height
barplot_height <- 3


# =================== SCRIPT START ===================
# Load required packages ---------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(ggprism)
  library(RColorBrewer)
  library(scales)  # For log transformation
})

# Options -----------------------------------------------------------------
facet_label_width <- 15  # Width for wrapping facet labels
barplot_color_palette <- c("Dark2", "Set2", "Set3")  # RColorBrewer palettes to use
heatmap_color_palette <- "RdYlBu"  # RColorBrewer palette for heatmap

# Prepare Output Path -----------------------------------------------------
stacked_barplot_output <- paste0(dirname(tax_binning_data_path),
                                 "/stacked_barplot-",
                                 basename(tax_binning_data_path),
                                 ".pdf")

heatmap_output <- paste0(dirname(tax_binning_data_path),
                         "/heatmap-",
                         basename(tax_binning_data_path),
                         ".pdf")

# Load and prepare data ---------------------------------------------------
# Load taxonomic binning data
message("Loading taxonomic binning data from: ", tax_binning_data_path)
tax_binning_dat <- read_delim(tax_binning_data_path, show_col_types = FALSE)
colnames(tax_binning_dat)[1] <- "Taxonomic_Bin"

# Clean taxonomic bin names (remove prefix)
tax_binning_dat$Taxonomic_Bin <- gsub(pattern = ".__", 
                                      replacement = "",
                                      x = tax_binning_dat$Taxonomic_Bin)

# Load sample metadata
message("Loading sample metadata from: ", sample_data_path)
sample_dat <- read_delim(sample_data_path, show_col_types = FALSE)

# Calculate dynamic width based on number of samples -----------------------
# Get number of samples
n_samples <- ncol(tax_binning_dat) - 1  # Subtract 1 for the Taxonomic_Bin column

# Get number of groups (for width calculation)
# Extract the facet variables from the sample data
facet_data <- sample_dat %>% 
  dplyr::select(all_of(c(id_col, facet_variables))) %>%
  rename(Sample = id_col)

# Count unique combinations of facet variables
n_groups <- facet_data %>%
  dplyr::select(-Sample) %>%
  distinct() %>%
  nrow()

# Calculate width based on number of samples and groups

# Calculate the total width 
plot_width <- base_width + (n_samples / n_groups * width_per_sample)
# Ensure minimum width
plot_width <- max(plot_width, 8)
# Cap maximum width to reasonable value
plot_width <- min(plot_width, max_plot_width)

message("Calculated dynamic plot width: ", round(plot_width, 1), 
        " inches (", n_samples, " samples across ", n_groups, " groups)")

# Generate stacked barplot -------------------------------------------------
message("Generating stacked barplot with top ", top_n_bins, " taxonomic bins")

# Get names of most abundant taxonomic bins
top_bins <- tax_binning_dat %>% 
  pivot_longer(!Taxonomic_Bin) %>% 
  group_by(Taxonomic_Bin) %>% 
  summarise(sum = sum(value)) %>% 
  arrange(desc(sum)) %>% 
  slice_head(n = top_n_bins) %>% 
  pull(Taxonomic_Bin)

# Create color palette for top bins and "other"
palette1 <- brewer.pal(min(8, top_n_bins), barplot_color_palette[1])
palette2 <- brewer.pal(min(7, max(0, top_n_bins - 8)), barplot_color_palette[2])
palette3 <- brewer.pal(min(7, max(0, top_n_bins - 15)), barplot_color_palette[3])
color_palette <- c(palette1, palette2, palette3)[1:top_n_bins]
names(color_palette) <- top_bins
color_palette <- c(color_palette, "other" = "grey")

# Create a facet formula from user-specified variables
facet_formula <- as.formula(paste("~", paste(facet_variables, collapse = " + ")))

# Create stacked barplot
barplot <- tax_binning_dat %>% 
  pivot_longer(!Taxonomic_Bin, names_to = "Sample", values_to = "Abundance") %>% 
  left_join(sample_dat, by = c("Sample" = id_col)) %>% 
  mutate(Bin = if_else(Taxonomic_Bin %in% top_bins, Taxonomic_Bin, "other")) %>% 
  ggplot(aes(x = Sample, y = Abundance, fill = reorder(Bin, Abundance))) +
  geom_bar(position = "stack", stat = "identity") +
  facet_grid(facet_formula, 
             labeller = label_wrap_gen(width = facet_label_width),
             scale = "free_x", space = "free_x") +
  scale_y_continuous(guide = "prism_offset") +
  scale_fill_manual(values = color_palette, name = "Taxonomic Bin") +
  labs(y = "Relative abundance (%)") +
  theme_prism() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  rremove("x.axis")

# Display and save the stacked barplot
message("Saving stacked barplot to: ", stacked_barplot_output)
# plot(barplot)
ggsave(plot = barplot, filename = stacked_barplot_output, width = plot_width, height = barplot_height)

# Generate heatmap --------------------------------------------------------
message("Generating heatmap of taxonomic abundance")

# Create heatmap with log-transformed color scale
heatmap <- tax_binning_dat %>% 
  pivot_longer(!Taxonomic_Bin, names_to = "Sample", values_to = "Abundance") %>% 
  left_join(sample_dat, by = c("Sample" = id_col)) %>% 
  # Add small value to avoid log(0) issues
  mutate(Abundance = Abundance + 0.01) %>% 
  ggplot(aes(y = reorder(Taxonomic_Bin, Abundance), x = Sample, fill = Abundance)) +
  geom_raster() +
  facet_grid(facet_formula,
             labeller = label_wrap_gen(width = facet_label_width),
             scale = "free_x") +
  # Apply log transformation to color scale
  scale_fill_distiller(palette = heatmap_color_palette, 
                       trans = "log10",
                       name = "% Relative\nAbundance\n(log scale)",
                       breaks = c(0.01, 0.1, 1, 10, 100),
                       labels = c("0.01", "0.10", "1.00", "10.00", "100")) +
  labs(y = "") +
  theme_prism() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text()) +
  rremove("x.axis")

# Calculate dynamic height for heatmap based on number of taxonomic bins
n_tax_bins <- nrow(tax_binning_dat)
height_per_bin <- 0.15  # Height per taxonomic bin
base_height <- 3  # Base height for margins and legend
heatmap_height <- base_height + (n_tax_bins * height_per_bin)
# Ensure minimum height
heatmap_height <- max(heatmap_height, 6)
# Cap maximum height to reasonable value
heatmap_height <- min(heatmap_height, 30)

message("Calculated dynamic heatmap height: ", round(heatmap_height, 1), 
        " inches (", n_tax_bins, " taxonomic bins)")

# Display and save the heatmap
message("Saving heatmap to: ", heatmap_output)
# plot(heatmap)
ggsave(plot = heatmap, filename = heatmap_output, width = plot_width, height = heatmap_height)

message("Script execution complete!")