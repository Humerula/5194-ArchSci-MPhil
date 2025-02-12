##################################Creating graphs for sample matching######################################
setwd("C:/Users/Desktop/DISS/graphs")

#Loading libraries
library(ggplot2)
library(viridis)
library(ggthemes) # for better colorblind-friendly themes
library(dplyr)
library(hrbrthemes)

#Loading dataset for creating graphs - tables with matching samples, individual time and distances plus combined distances
plague_graphs <- read.csv("C:/Users/Desktop/DISS/graphs/plague-human_forgraphs.csv")
leprosy_graphs <- read.csv("C:/Users/Desktop/DISS/graphs/leprosy-human_forgraphs.csv")
pathogen_graphs <- read.csv("C:/Users/Desktop/DISS/graphs/bothpathogenscombined_forgraphs.csv")


################ PLOTTING REAL DISTANCE VS EUCLIDIAN DIST ONTO SCATTER PLOT##############################
# Function to plot time difference vs. geographic distance against Euclidean distance
plot_distances <- function(closest_samples, title) {
    ggplot(closest_samples, aes(x = temp_dist, y = geo_dist / 1000, color = combined_dist)) +
        geom_point(size = 2, alpha = 0.7) +  # Increase point size and transparency
        scale_color_viridis_c(option = "plasma") +  # Better color scale
        labs(
            title = title,
            x = "Temporal Distance (years)",
            y = "Geographic Distance (km)",
            color = "Euclidean Distance"
        ) +
        theme_bw(base_size = 12) +  # Improve readability with larger text
}

# Save high-quality plots
png(
    file = "C:/Users/Desktop/DISS/graphs/leprosy_realdistance_vs_euclidean_feb.png",
    width = 6000, height = 4800, res = 800
)
plot_distances(leprosy_graphs, "Leprosy: Temporal vs. Geographic Distance")
dev.off()

png(
    file = "C:/Users/Desktop/DISS/graphs/plague_realdistance_vs_euclidean_feb.png",
    width = 6000, height = 4800, res = 800
)
plot_distances(plague_graphs, "Plague: Temporal vs. Geographic Distance")
dev.off()

# Save as SVG to combine
svg(file = "C:/Users/Desktop/DISS/graphs/leprosy_realdistance_vs_euclidean_feb2.svg")
plot_distances(leprosy_graphs, "Leprosy: Temporal vs. Geographic Distance")
dev.off()


svg(file = "C:/Users/Desktop/DISS/graphs/plague_realdistance_vs_euclidean_feb2.svg")
plot_distances(plague_graphs, "Plague: Temporal vs. Geographic Distance")
dev.off()

####################CREATING HISTOGRAMS FOR BOTH PATHOGENS TOGETHER#########################

######Histogram of time differences for plague and for leprosy, plotted on a single axis########

# Define colorblind-friendly palette
cbPalette <- c("Plague" = "#0072B2","Leprosy" = "#D55E00")  # Orange for leprosy, blue for plague
# Save high-res PNG
png(
    filename = "C:/Users/Desktop/DISS/graphs/time_diff_density_plot.png",
    width = 6500, height = 4800, res = 800  # High resolution
)
# Density plot for time differences
ggplot(pathogen_graphs, aes(x = temp_dist, fill = pathogen_type, color = pathogen_type)) +
    geom_density(alpha = 0.3, linewidth = 1.2) +  # Density curve with transparency
    scale_fill_manual(values = cbPalette) + 
    scale_color_manual(values = cbPalette) + 
    theme_minimal(base_size = 16) +
    theme(legend.title = element_blank(), legend.position = "right") +
    labs(
        title = "Density Plot of Time Differences for Plague snd Leprosy",
        x = "Time Difference (years)", 
        y = "Density"
    )

dev.off()  # Close the graphics device

# Define custom colors (Plague = Blue, Leprosy = Orange)
cbPalette <- c("Plague" = "#0072B2", "Leprosy" = "#D55E00")  

# Create the histogram with dodge position (side-by-side bars)
ggplot(pathogen_graphs, aes(x = temp_dist, fill = pathogen_type)) +
    geom_histogram(binwidth = 10, color = "black", alpha = 0.7, 
                   position = position_dodge(width = 10)) + 
    scale_fill_manual(values = cbPalette) +  # Set custom colors
    theme_minimal(base_size = 18) +  # Increase base font size
    labs(title = "Histogram of Time Differences for Plague and Leprosy", 
         x = "Time Difference (years)", 
         y = "Frequency", 
         fill = "Pathogen Type") +
    theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 22, margin = margin(b = 15)),  # Center title & move it higher
        axis.title.x = element_text(hjust = 0.5, size = 20, face = "bold"),  # Center & enlarge x-label
        axis.title.y = element_text(hjust = 0.5, size = 20, face = "bold"),  # Center & enlarge y-label
        axis.text = element_text(size = 16)  # Enlarge axis numbers
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.02)))  # Add some spacing

ggsave("C:/Users/Desktop/DISS/graphs/time_diff_histogram2.png", width = 10, height = 6, dpi = 300)


######Histogram of spatial differences for plague and for leprosy, plotted on a single axis########

# Define custom colors (Plague = Blue, Leprosy = Orange)
cbPalette <- c("Plague" = "#0072B2", "Leprosy" = "#D55E00")  

# Create the histogram with dodge position (side-by-side bars)
ggplot(pathogen_graphs, aes(x = geo_dist_km, fill = pathogen_type)) +
    geom_histogram(binwidth = 10, color = "black", alpha = 0.7, 
                   position = position_dodge(width = 10)) + 
    scale_fill_manual(values = cbPalette) +  # Set custom colors
    theme_minimal(base_size = 18) +  # Increase base font size
    labs(title = "Histogram of Spatial Differences for Plague and Leprosy", 
         x = "Spatial Distance (km)", 
         y = "Frequency", 
         fill = "Pathogen Type") +
    theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 22, margin = margin(b = 15)),  # Center title & move it higher
        axis.title.x = element_text(hjust = 0.5, size = 20, face = "bold"),  # Center & enlarge x-label
        axis.title.y = element_text(hjust = 0.5, size = 20, face = "bold"),  # Center & enlarge y-label
        axis.text = element_text(size = 16)  # Enlarge axis numbers
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.02)))  # Add some spacing

ggsave("C:/Users/Desktop/DISS/graphs/geo_diff_histogram.png", width = 20, height = 12, dpi = 300)


##Define colour palette
cbPalette <- c("Plague" = "#0072B2","Leprosy" = "#D55E00")  # Orange for leprosy, blue for plague

# Save high-res PNG
png(
    filename = "C:/Users/Desktop/DISS/graphs/geog_diff_density_plot.png",
    width = 6500, height = 4800, res = 800  # High resolution
)
# Density plot for time differences
ggplot(pathogen_graphs, aes(x = geo_dist_km, fill = pathogen_type, color = pathogen_type)) +
    geom_density(alpha = 0.3, linewidth = 1.2) +  # Density curve with transparency
    scale_fill_manual(values = cbPalette) + 
    scale_color_manual(values = cbPalette) + 
    theme_minimal(base_size = 16) +
    theme(legend.title = element_blank(), legend.position = "right") +
    labs(
        title = "Density Plot of Spatial Differences for Plague and Leprosy",
        x = "Spatial Difference (km)", 
        y = "Density"
    )

dev.off()  # Close the graphics device


################## Histogram of Euclidean distance metric #####################
summary(pathogen_graphs$combined_dist)

# Define custom colors (Plague = Blue, Leprosy = Orange)
cbPalette <- c("Plague" = "#0072B2", "Leprosy" = "#D55E00")  

# Create the histogram with dodge position (side-by-side bars)
ggplot(pathogen_graphs, aes(x = combined_dist, fill = pathogen_type)) +
    geom_histogram(binwidth = 0.005, color = "black", alpha = 0.7, 
                   position = position_dodge(width = 0.005)) + 
    scale_fill_manual(values = cbPalette) +  # Set custom colors
    theme_minimal(base_size = 18) +  # Increase base font size
    labs(title = "Histogram of Euclidean Distances for Plague and Leprosy", 
         x = "Combined Eucledian Distance", 
         y = "Frequency", 
         fill = "Pathogen Type") +
    theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 22, margin = margin(b = 15)),  # Center title & move it higher
        axis.title.x = element_text(hjust = 0.5, size = 20, face = "bold"),  # Center & enlarge x-label
        axis.title.y = element_text(hjust = 0.5, size = 20, face = "bold"),  # Center & enlarge y-label
        axis.text = element_text(size = 16)  # Enlarge axis numbers
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.02)))  # Add some spacing

ggsave("C:/Users/Desktop/DISS/graphs/euc_diff_histogram.png", width = 30, height = 12, dpi = 300)


##Define colour palette
cbPalette <- c("Plague" = "#0072B2","Leprosy" = "#D55E00")  # Orange for leprosy, blue for plague

# Save high-res PNG
png(
    filename = "C:/Users/Desktop/DISS/graphs/euc_diff_density_plot.png",
    width = 6500, height = 4800, res = 800  # High resolution
)
# Density plot for combined differences
ggplot(pathogen_graphs, aes(x = combined_dist, fill = pathogen_type, color = pathogen_type)) +
    geom_density(alpha = 0.3, linewidth = 1.2) +  # Density curve with transparency
    scale_fill_manual(values = cbPalette) + 
    scale_color_manual(values = cbPalette) + 
    theme_minimal(base_size = 16) +
    theme(legend.title = element_blank(), legend.position = "right") +
    labs(
        title = "Density Plot of Euclidean Distances for Plague and Leprosy",
        x = "Combined Euclidean Distance", 
        y = "Density"
    )

dev.off()  # Close the graphics device

############# Scatter plot: Sample time vs. Time difference with correlation ###########
png(
    filename = "C:/Users/Desktop/DISS/graphs/scatter_sampletime_vs_time_diff3.png",
    width = 6500, height = 4800, res = 800  # High resolution
)

svg(filename = "C:/Users/Desktop/DISS/graphs/scatter_sampletime_vs_time_diff3.svg")

ggplot(pathogen_graphs, aes(x = age_pathogen, y = temp_dist, color = pathogen_type)) +
    geom_point(alpha = 0.7, size = 3) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
    scale_color_manual(values = cbPalette) +
    labs(title = "Scatter Plot: Pathogen Age vs. Time Difference",
         x = "Pathogen Age (years BP)", y = "Time Difference (years)") +
    theme_minimal(base_size = 16) +
    theme(legend.title = element_blank())

dev.off()

cor.test(pathogen_graphs$age_pathogen, pathogen_graphs$temp_dist)
##Pearson's product-moment correlation
#data:  pathogen_graphs$age_pathogen and pathogen_graphs$temp_dist
#t = -3.1651, df = 73, p-value = 0.002262
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.5323667 -0.1307219
#sample estimates:
#       cor
#-0.3473774
