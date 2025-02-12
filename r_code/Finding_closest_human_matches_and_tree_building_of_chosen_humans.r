library(readr)
library(readxl)
library(geosphere)
library(tidyverse)
library(dplyr)
library(stringi)
library(ggplot2)
library(janitor)
library(scatterplot3d)

setwd("C:/Users/Desktop/DISS/Human_closest")

## Reading Allen data
anno_data <- read_delim("C:/Users/Desktop/DISS/v54.1.p1_1240K_public.anno", delim = "\t")

head(anno_data)
## selecting variables of interest
anno_var <- c("Master ID", "Lat", "Long", "Date_mean_B1950")
subset_anno <- anno_data[anno_var]

## convert the lat and long into floating pts
subset_anno$Lat <- as.numeric(subset_anno$Lat)
subset_anno$Long <- as.numeric(subset_anno$Long)

class(subset_anno$Long)

# Handle invalid numbers (e.g., NAs)
subset_anno_clean <- subset_anno[complete.cases(subset_anno$Lat), ]
subset_anno_clean <- subset_anno_clean[complete.cases(subset_anno_clean$Long), ]
subset_anno_clean <- subset_anno_clean[complete.cases(subset_anno_clean$Date_mean_B1950), ]

################################# FOR ANCIENT PATHOGENS ONLY############################################

# Loading pathogen data
leprosy_metadata <- read.csv("C:/Users/Desktop/DISS/Human_closest/lep_metadata_for_closest.csv")
plague_metadata <- read.csv("C:/Users/Desktop/DISS/Human_closest/plague_metadata_for_closest.csv")

# Ensuring dates are numeric (representing ages in years ago from 1950)
leprosy_metadata$Age_bp_1950 <- as.numeric(leprosy_metadata$Age_bp_1950)
plague_metadata$Age_bp_1950 <- as.numeric(plague_metadata$Age_bp_1950)
subset_anno_clean$Date_mean_B1950 <- as.numeric(subset_anno_clean$Date_mean_B1950)

# renaming human date column to match others
subset_anno_clean <- subset_anno_clean %>%
    rename(Age_bp_1950 = Date_mean_B1950)

#############################################################################

# Defining a function to filter human samples based on pathogen metadata
filter_human_samples <- function(pathogen_metadata, human_metadata, time_buffer = 200, lat_buffer = 3, long_buffer = 3) {
    # Calculate time bounds
    min_date <- min(pathogen_metadata$Age_bp_1950) - time_buffer
    max_date <- max(pathogen_metadata$Age_bp_1950) + time_buffer

    # Calculate geographic bounds
    min_lat <- min(pathogen_metadata$Lat) - lat_buffer
    max_lat <- max(pathogen_metadata$Lat) + lat_buffer
    min_long <- min(pathogen_metadata$Long) - long_buffer
    max_long <- max(pathogen_metadata$Long) + long_buffer

    # Filter human samples
    human_metadata %>%
        filter(
            Age_bp_1950 >= min_date & Age_bp_1950 <= max_date,
            Lat >= min_lat & Lat <= max_lat,
            Long >= min_long & Long <= max_long
        )
}

# Filter human samples for both leprosy and plague
filtered_human_leprosy <- filter_human_samples(leprosy_metadata, subset_anno_clean)
filtered_human_plague <- filter_human_samples(plague_metadata, subset_anno_clean)

#############################################################################
# Function to calculate distances and find the closest human sample
find_closest_human_sample <- function(pathogen_metadata, human_metadata, temporal_weight = 0.4, geo_weight = 0.6) {
    # Initialize an empty dataframe to store results
    closest_samples <- data.frame()

    # Calculating maximum distances for normalization
    max_geo_dist <- 0
    max_temp_dist <- 0

    for (i in 1:nrow(pathogen_metadata)) {
        for (j in 1:nrow(human_metadata)) {
            geo_dist <- distHaversine(
                c(pathogen_metadata$Long[i], pathogen_metadata$Lat[i]),
                c(human_metadata$Long[j], human_metadata$Lat[j])
            )
            temp_dist <- abs(pathogen_metadata$Age_bp_1950[i] - human_metadata$Age_bp_1950[j])
            max_geo_dist <- max(max_geo_dist, geo_dist)
            max_temp_dist <- max(max_temp_dist, temp_dist)
        }
    }

    # Calculating normalized distances and combined Euclidean distance
    for (i in 1:nrow(pathogen_metadata)) {
        min_distance <- Inf
        closest_sample <- NULL

        for (j in 1:nrow(human_metadata)) {
            # Calculate normalised geo distance
            geo_dist <- distHaversine(
                c(pathogen_metadata$Long[i], pathogen_metadata$Lat[i]),
                c(human_metadata$Long[j], human_metadata$Lat[j])
            ) / max_geo_dist

            # Calculate normalised temporal distance
            temp_dist <- abs(pathogen_metadata$Age_bp_1950[i] - human_metadata$Age_bp_1950[j]) / max_temp_dist

            # Calculate combined euclidean dist with weights
            # could account for human migration like paper X
            combined_dist <- sqrt((geo_weight * geo_dist)^2 + (temporal_weight * temp_dist)^2)

            if (combined_dist < min_distance) {
                min_distance <- combined_dist
                closest_sample <- human_metadata[j, ]
                closest_sample$combined_dist <- combined_dist
                closest_sample$geo_dist <- geo_dist * max_geo_dist
                closest_sample$temp_dist <- temp_dist * max_temp_dist
            }
        }

        closest_samples <- rbind(closest_samples, cbind(pathogen_metadata[i, ], closest_sample))
    }

    return(closest_samples)
}

# Find closest human samples for both leprosy and plague
closest_human_leprosy <- find_closest_human_sample(leprosy_metadata, filtered_human_leprosy, temporal_weight = 0.4, geo_weight = 0.6)
closest_human_plague <- find_closest_human_sample(plague_metadata, filtered_human_plague, temporal_weight = 0.4, geo_weight = 0.6)

# Ensure unique column names in the resulting dataframes

cleaned_closest_human_leprosy <- janitor::clean_names(closest_human_leprosy)
cleaned_closest_human_plague <- janitor::clean_names(closest_human_plague)

cleaned_closest_human_leprosy <- cleaned_closest_human_leprosy %>%
    rename(
        sample_leprosy = sample,
        lat_leprosy = lat,
        long_leprosy = long,
        age_leprosy = age_bp_1950,
        sample_human_leprosy = master_id,
        lat_human_leprosy = lat_2,
        long_human_leprosy = long_2,
        age_human_leprosy = age_bp_1950_2
    )

cleaned_closest_human_plague <- cleaned_closest_human_plague %>%
    rename(
        sample_plague = sample,
        lat_plague = lat,
        long_plague = long,
        age_plague = age_bp_1950,
        sample_human_plague = master_id,
        lat_human_plague = lat_2,
        long_human_plague = long_2,
        age_human_plague = age_bp_1950_2
    )

# Write clean dataframe to csv
write.csv(cleaned_closest_human_leprosy, "C:/Users/Desktop/DISS/Human_closest/leprosy_closest_human_NEW.csv", row.names = FALSE)
write.csv(cleaned_closest_human_plague, "C:/Users/Desktop/DISS/Human_closest/plague_closest_human_febfinal.csv", row.names = FALSE)

###################### CREATING TREES FOR HUMANS CLOSE TO LEPROSY#######################################

## create trees from closest samples for leprosy
library(Rcpp)
library(tidyverse)
library(igraph)
library(plotly)
library(devtools)
library(admixtools)
library(snpStats)
library(ape)
library(genio)

# is working only if .ind file has same names, check carefully and add .ind file type names instead of those generated for closest match ID
humanindfile <- read.table("C:/Users/Desktop/DISS/v54.1.p1_1240K_public.ind", header = FALSE)

human_leprosy_names <- list(cleaned_closest_human_leprosy$sample_human_leprosy)

# list of inds matched
inds_leprosy <- c(
    "JK2888",
    "I3044",
    "VK256_noUDG.SG",
    "I3044",
    "SWG002.A",
    "SWG002.A",
    "SWG002.A",
    "SWG001.A",
    "SWG001.A",
    "SWG002.A",
    "SWG002.A",
    "vik_urm161_noUDG.SG",
    "VK330_noUDG.SG",
    "I3044",
    "VK256_noUDG.SG",
    "VK256_noUDG.SG",
    "VK252_noUDG.SG",
    "VK414_noUDG.SG",
    "I3044",
    "I3044",
    "I16510",
    "I3044",
    "VK207_noUDG.SG",
    "I8145",
    "S_Estonian-2.DG",
    "HGDP01364.DG",
    "I12515",
    "I12515",
    "I10853",
    "I10897",
    "I10897",
    "I10853",
    "I0774",
    "VEN005",
    "FU-215.SG",
    "SZAK-4.SG",
    "VK256_noUDG.SG"
)

# Check if `inds` are all present in `indfile$X1`
if (!all(inds_leprosy %in% humanindfile$V1)) {
    missing_individuals <- setdiff(inds_leprosy, humanindfile$V1)
    stop(paste("The following individuals are missing in the genotype file:", paste(missing_individuals, collapse = ", ")))
}

eigenstrat_to_plink(
    inpref = "C:/Users/Desktop/DISS/v54.1.p1_1240K_public",
    outpref = "C:/Users/Desktop/DISS/Human_closest/v54.1.p1_1240K_public_plink_leprosymatch2907",
    inds = c(
        "JK2888",
        "I3044",
        "VK256_noUDG.SG",
        "I3044",
        "SWG002.A",
        "SWG002.A",
        "SWG002.A",
        "SWG001.A",
        "SWG001.A",
        "SWG002.A",
        "SWG002.A",
        "vik_urm161_noUDG.SG",
        "VK330_noUDG.SG",
        "I3044",
        "VK256_noUDG.SG",
        "VK256_noUDG.SG",
        "VK252_noUDG.SG",
        "VK414_noUDG.SG",
        "I3044",
        "I3044",
        "I16510",
        "I3044",
        "VK207_noUDG.SG",
        "I8145",
        "S_Estonian-2.DG",
        "HGDP01364.DG",
        "I12515",
        "I12515",
        "I10853",
        "I10897",
        "I10897",
        "I10853",
        "I0774",
        "VEN005",
        "FU-215.SG",
        "SZAK-4.SG",
        "VK256_noUDG.SG"
    ),
    verbose = TRUE
)


# To make distance matrices for APE
plink_path <- "C:/Users/Desktop/DISS/plink_win64_20231211/plink.exe"
input_humans <- "C:/Users/Desktop/DISS/Human_closest/v54.1.p1_1240K_public_plink_leprosymatch2907"
output_dist_mat_humans <- "C:/Users/Desktop/DISS/Human_closest/human_dist_mat_leprosymatch2907"


# Define the PLINK command
command_distmat_hum <- paste0(plink_path, " --bfile ", input_humans, " --distance square 1-ibs --out ", output_dist_mat_humans)

# Run the command using system
system(command_distmat_hum)


# Read the distance matrix into R
dist_matrix_human_leprosy <- as.matrix(read.table("C:/Users/Desktop/DISS/Human_closest/human_dist_mat_leprosymatch2907.mdist"))


# Find the largest value in the distance matrix
# Since there are values not larger than 100 which bionj does not like so no division needed
largest_value <- max(dist_matrix_human_leprosy, na.rm = TRUE)

print(paste("The largest value in the distance matrix is:", largest_value))


# Read the .fam file into R
human_fam_for_samplenames_leprosy <- read.table("C:/Users/Desktop/DISS/Human_closest/v54.1.p1_1240K_public_plink_leprosymatch2907.fam")

# Extract the sample IDs
sample_ids <- human_fam_for_samplenames_leprosy$V2

# Replace the row names and column names of the distance matrix
rownames(dist_matrix_human_leprosy) <- sample_ids
colnames(dist_matrix_human_leprosy) <- sample_ids


# Compute the BIONJ tree
tree_human_leprosy <- bionj(dist_matrix_human_leprosy)

# Print the BIONJ tree
plot(tree_human_leprosy)

# Export the BIONJ tree
png(
    file = "C:/Users/Desktop/DISS/Human_closest/human_tree_leprosy_ancientNEW.png",
    width = 2000, height = 4500
)

plot(tree_human_leprosy, main = "Human Tree for best Ancient Leprosy matches")

dev.off()

write.tree(tree_human_leprosy, file = "C:/Users/Desktop/DISS/Human_closest/human_tree_leprosy_ancientNEW.txt")

read.tree("C:/Users/Desktop/DISS/Human_closest/human_tree_leprosy_ancientNEW.txt")

############### CREATING TREES FOR HUMANS CLOSE TO PLAGUE##########################
# creating tree for closest samples to plague
library(Rcpp)
library(tidyverse)
library(igraph)
library(plotly)
library(devtools)
library(admixtools)
library(snpStats)
library(ape)
library(genio)


# is working only if .ind file has same names, check .anno and carefully and update to match name to .anno/.ind and not that generated as close match ID
humanindfile <- read.table("C:/Users/Desktop/DISS/v54.1.p1_1240K_public.ind", header = FALSE)

human_plague_names <- list(cleaned_closest_human_plague$sample_human_plague)


# updated ind names
inds_plague <- c(
    "I7205",
    "RK4001",
    "I0117",
    "VLI047",
    "VLI047",
    "NV3001",
    "MIB034",
    "OTTM_81_d",
    "MIB034",
    "OTTM_156_noUDG",
    "OTTM_91_d",
    "ARS006.A0101",
    "KZL003",
    "I5037",
    "DA145_noUDG.SG",
    "HGDP00530.DG",
    "HGDP00530.DG",
    "HGDP00530.DG",
    "HGDP00530.DG",
    "HGDP00530.DG",
    "HG00126.SG",
    "HG00126.SG",
    "HG00126.SG",
    "HG00126.SG",
    "R2066.SG",
    "R2066.SG",
    "I17313",
    "Sunghir6_noUDG.SG",
    "HG00126.SG",
    "R2066.SG",
    "I17313",
    "R2066.SG",
    "R2066.SG",
    "R2066.SG",
    "I3044",
    "ber1M_noUDG.SG",
    "RK1007",
    "RISE684_noUDG.SG",
    "RISE512_noUDG.SG",
    "I10897",
    "I2741",
    "Spiginas2",
    "OTTM_152_d",
    "WEHR_1192SkB",
    "RK4001",
    "EKA1_noUDG.SG",
    "I0359",
    "GLZ001",
    "GLZ002",
    "Sunghir6_noUDG.SG",
    "AED125b_noUDG.SG",
    "R2066.SG",
    "R2066.SG",
    "R2066.SG",
    "I0773"
)

# Check if `inds` are all present in `indfile$X1`
if (!all(inds_plague %in% humanindfile$V1)) {
    missing_individuals <- setdiff(inds_plague, humanindfile$V1)
    stop(paste("The following individuals are missing in the genotype file:", paste(missing_individuals, collapse = ", ")))
}

eigenstrat_to_plink(
    inpref = "C:/Users/Desktop/DISS/v54.1.p1_1240K_public",
    outpref = "C:/Users/Desktop/DISS/Human_closest/v54.1.p1_1240K_public_plink_plaguematchFEB",
    inds = c(
        "I7205",
        "RK4001",
        "I0117",
        "VLI047",
        "VLI047",
        "NV3001",
        "MIB034",
        "OTTM_81_d",
        "MIB034",
        "OTTM_156_noUDG",
        "OTTM_91_d",
        "ARS006.A0101",
        "KZL003",
        "I5037",
        "DA145_noUDG.SG",
        "HGDP00530.DG",
        "HGDP00530.DG",
        "HGDP00530.DG",
        "HGDP00530.DG",
        "HGDP00530.DG",
        "HG00126.SG",
        "HG00126.SG",
        "HG00126.SG",
        "HG00126.SG",
        "R2066.SG",
        "R2066.SG",
        "I17313",
        "Sunghir6_noUDG.SG",
        "HG00126.SG",
        "R2066.SG",
        "I17313",
        "R2066.SG",
        "R2066.SG",
        "R2066.SG",
        "I3044",
        "ber1M_noUDG.SG",
        "RK1007",
        "RISE684_noUDG.SG",
        "RISE512_noUDG.SG",
        "I10897",
        "I2741",
        "Spiginas2",
        "OTTM_152_d",
        "WEHR_1192SkB",
        "RK4001",
        "EKA1_noUDG.SG",
        "I0359",
        "GLZ001",
        "GLZ002",
        "Sunghir6_noUDG.SG",
        "AED125b_noUDG.SG",
        "R2066.SG",
        "R2066.SG",
        "R2066.SG",
        "I0773"
    ),
    verbose = TRUE
)


# To make distance matrices for APE
plink_path <- "C:/Users/Desktop/DISS/plink_win64_20231211/plink.exe"
input_humans <- "C:/Users/Desktop/DISS/Human_closest/v54.1.p1_1240K_public_plink_plaguematchFEB"
output_dist_mat_humans <- "C:/Users/Desktop/DISS/Human_closest/human_dist_mat_plaguematchFEB"


# Define the PLINK command
command_distmat_hum <- paste0(plink_path, " --bfile ", input_humans, " --distance square 1-ibs --out ", output_dist_mat_humans)

# Run the command using system
system(command_distmat_hum)


# Read the distance matrix into R
dist_matrix_human_plague <- as.matrix(read.table("C:/Users/Desktop/DISS/Human_closest/human_dist_mat_plaguematchFEB.mdist"))


# Read the .fam file into R
human_fam_for_samplenames_plague <- read.table("C:/Users/Desktop/DISS/Human_closest/v54.1.p1_1240K_public_plink_plaguematchFEB.fam")

# Extract the sample IDs
sample_ids <- human_fam_for_samplenames_plague$V2

# Replace the row names and column names of the distance matrix
rownames(dist_matrix_human_plague) <- sample_ids
colnames(dist_matrix_human_plague) <- sample_ids

dist_matrix_human_plague

## clean rows with invalid names
dist_matrix_human_plague_cleanrow <- dist_matrix_human_plague[!rownames(dist_matrix_human_plague) %in% c("I17313", "R2066.SG", "HG00126.SG", "I2741"), ]
dist_matrix_human_plague_cleanall <- dist_matrix_human_plague_cleanrow[, !colnames(dist_matrix_human_plague_cleanrow) %in% c("I17313", "R2066.SG", "HG00126.SG", "I2741")]

# Find the largest value in the distance matrix
# if there are values larger than 1 which bionj does not like ONLY THEN divide
largest_value <- max(dist_matrix_human_plague_cleanall, na.rm = TRUE)

print(paste("The largest value in the distance matrix is:", largest_value))

# diving by largest value as bionj does not like resultant values above 100 and reading as matrix
dist_matrix_human_plague <- dist_matrix_human_plague / 2.59869
largest_value <- max(dist_matrix_human_plague, na.rm = TRUE)
print(paste("The largest value in the distance matrix is:", largest_value))

write.table(dist_matrix_human_plague_cleanall, file = "C:/Users/Desktop/DISS/Human_closest/human_dist_mat_plaguematchFEB_clean.mdist", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

# Compute the BIONJ tree#
tree_human_plague <- bionj(dist_matrix_human_plague_cleanall)

# Check for negative branch lengths and set them to zero
negative_branches <- tree_human_plague$edge.length < 0
if (any(negative_branches)) {
    cat("Negative branch lengths found:\n")
    print(tree_human_plague$edge.length[negative_branches])
    cat("Setting negative branch lengths to zero.\n")
    tree_human_plague$edge.length[negative_branches] <- 0
} else {
    cat("No negative branch lengths found.\n")
}


# Print the BIONJ tree
plot(tree_human_plague)

# Export the BIONJ tree
png(
    file = "C:/Users/Desktop/DISS/Human_closest/human_tree_plague_ancientNEW.png",
    width = 2000, height = 4500
)

plot(tree_human_plague, main = "Human Tree for best Ancient Plague matches")

dev.off()

write.tree(tree_human_plague, file = "C:/Users/Desktop/DISS/Human_closest/human_tree_plague_ancientFEB.txt")

read.tree("C:/Users/Desktop/DISS/Human_closest/human_tree_plague_ancientFEB.txt")
#####################################################
