library(readr)
library(readxl)
library(geosphere)
library(tidyverse)
library(dplyr)
library(stringi)
library(ape)
library(vcfR)
library(admixtools)
library(snpStats)
library(genio)
library(data.table)
library(phytools)

############################################## URBAN ET AL DATA#######################################################
# input URBAN ET AL 2024 DATA from SUPERVISOR- haploids only with single alleles and convert to Plink

## Changing sample IDs for problematic ones and removing samples first

# Read the VCF file
vcf_path <- "C:/Users/Desktop/DISS/squirrel_leprosy/Urban2024_leprosy_snp_table.haploid.biallelicOnly.vcf"
haploid_Urban_lep <- read.vcfR(vcf_path)

# Extract the sample IDs from the VCF file
sample_ids <- colnames(haploid_Urban_lep@gt)[-1] # Exclude the FORMAT column
cat("Original sample IDs:", sample_ids, "\n")

# Create a data table with the current and new IDs
id_data <- data.table(
    current_id = c("column_13", "column_14", "column_15", "column_16"), # The current IDs you want to update
    new_id = c("sample_2936", "sample_3077", "sample_85054", "sample_97016") # The new IDs
)

# Create a named vector for mapping current IDs to new IDs
id_map <- setNames(id_data$new_id, id_data$current_id)

# Print the mapping for debugging
cat("ID Mapping:\n")
print(id_map)

# Update the sample IDs in the VCF file
updated_ids <- sample_ids # Copy the original IDs
for (i in seq_along(sample_ids)) {
    if (sample_ids[i] %in% names(id_map)) {
        updated_ids[i] <- id_map[sample_ids[i]]
    }
}
cat("Updated sample IDs:", updated_ids, "\n")

# Apply the updated IDs to the VCF object
colnames(haploid_Urban_lep@gt)[-1] <- updated_ids

# Define the non-human samples to remove
samples_to_remove <- c("Brw15-10M2", "Brw15-5E", "Brw15-25E", "Brw15-20M2", "Brw15-1E", "Brw15-12M2", "I30_W09", "Ch4", "CM1", "SM1", "GB64", "TNP418")

# Remove the specified samples
keep_samples <- !colnames(haploid_Urban_lep@gt)[-1] %in% samples_to_remove
haploid_Urban_lep@gt <- haploid_Urban_lep@gt[, c(TRUE, keep_samples)] # Keep the FORMAT column and specified samples

# Write the updated VCF file to a new file
updated_vcf_path <- "C:/Users/Desktop/DISS/squirrel_leprosy/Urban2024_leprosy_snp_table.haploid.biallelicOnly.updated.vcf"
write.vcf(haploid_Urban_lep, updated_vcf_path)

# Verify the updated VCF file
updated_vcf <- read.vcfR(updated_vcf_path)
cat("Updated sample IDs in the new VCF:", colnames(updated_vcf@gt)[-1], "\n")

# Define the path to the PLINK executable
plink_path <- "C:/Users/Desktop/DISS/plink_win64_20231211/plink.exe"

# Define the base name for the output files
out_path <- "C:/Users/Desktop/DISS/squirrel_leprosy/Urban2024_leprosy_snp_table.haploid.biallelicOnly.plink"

# Construct the PLINK command
plink_command <- sprintf(
    "%s --vcf %s --double-id --make-bed --out %s --chr-set -1",
    plink_path, updated_vcf_path, out_path
)

# Print the command to check it (optional)
cat("Running command:", plink_command, "\n")

# Execute the PLINK command
system(plink_command)


# To make distance matrices for APE
plink_path <- "C:/Users/Desktop/DISS/plink_win64_20231211/plink.exe"
input_filtered <- "C:/Users/Desktop/DISS/squirrel_leprosy/Urban2024_leprosy_snp_table.haploid.biallelicOnly.plink"
output_dist_mat <- "C:/Users/Desktop/DISS/squirrel_leprosy/Urban2024_leprosy_distance_matrix_filtered"


# Define the distance matrix command
command_distmat <- paste0(plink_path, " --bfile ", input_filtered, " --distance square 1-ibs --out ", output_dist_mat)

# Run the command using system
system(command_distmat)


# Read the distance matrix into R
dist_matrix_urban <- as.matrix(read.table("C:/Users/Desktop/DISS/squirrel_leprosy/Urban2024_leprosy_distance_matrix_filtered.MDIST"))

# Find the largest value in the distance matrix
# Since there are values larger than 100 which bionj does not like
largest_value <- max(dist_matrix_urban, na.rm = TRUE)

print(paste("The largest value in the distance matrix is:", largest_value))

## not diving by largest value since it is all under 1 and while reading as matrix values above 100 are not bionj compatible


# Read the .fam file into R
fam_for_samplenames <- read.table("C:/Users/Desktop/DISS/squirrel_leprosy/Urban2024_leprosy_snp_table.haploid.biallelicOnly.plink.fam")

# Extract the sample IDs
sample_ids <- fam_for_samplenames$V2

# Replace the row names and column names of the distance matrix
rownames(dist_matrix_urban) <- sample_ids
colnames(dist_matrix_urban) <- sample_ids


# Compute the BIONJ tree
tree <- bionj(dist_matrix_urban)

# Check for negative branch lengths and set them to zero
negative_branches <- tree$edge.length < 0
if (any(negative_branches)) {
    cat("Negative branch lengths found:\n")
    print(tree$edge.length[negative_branches])
    cat("Setting negative branch lengths to zero.\n")
    tree$edge.length[negative_branches] <- 0
} else {
    cat("No negative branch lengths found.\n")
}

# exporting for visualisation using FigTree as R code below does not allow for 'vertical expansion' between branches without skewing tree
write.tree(tree, file = "C:/Users/Desktop/DISS/squirrel_leprosy/Urban2024_leprosy_distance_matrix_filtered.tree")


# Improve tree visualization
tree <- ladderize(tree) # Arrange branches for readability
par(mar = c(4, 4, 4, 4)) # Adjust margins


# Save high-resolution PNG
png("C:/Users/Desktop/DISS/squirrel_leprosy/improved_tree.png",
    width = 3508, height = 4961, res = 300
) # Larger A3 size for better spacing
plot(tree, main = "Leprosy Phylogeny Tree from URBAN 2024", cex = 0.5) # Reduce label size
dev.off()

# Save as SVG
svg("C:/Users/Desktop/DISS/squirrel_leprosy/improved_tree.svg",
    width = 11.69, height = 16.53
) # A3 size in inches
plot(tree, main = "Leprosy Phylogeny Tree from URBAN 2024", cex = 0.5)
dev.off()

########################## CREATING TREE FOR ONLY ANCIENT SAMPLES###############################################################

# Define the path to the updated VCF file
updated_vcf_path <- "C:/Users/Desktop/DISS/squirrel_leprosy/Urban2024_leprosy_snp_table.haploid.biallelicOnly.updated.vcf"

# Read the updated VCF file
updated_vcf <- read.vcfR(updated_vcf_path)

# Define the subset of samples to include
subset_samples <- c(
    "Abusir1630", "MMW_H94_1", "MMW_H80_1", "MMW_H50_1", "Jorgen_625", "Jorgen_404", "Jorgen_427",
    "Jorgen_507", "Jorgen_533", "Jorgen_722", "Jorgen_749", "sample_3077", "Refshale_16", "SK2", "SK8", "SK14", "BEL024",
    "Bergen", "CHRY023", "CHRY044", "EDI006", "JDS097", "KirkHill", "PAVd09_1.5", "R7546-671", "UF11", "UF21", "UF25",
    "UF101", "UF700", "UF703", "UF803", "GC96CU", "T18", "SK11", "Body_188", "SK27"
)

# Filter to keep only the subset of samples
keep_samples <- colnames(updated_vcf@gt)[-1] %in% subset_samples
filtered_vcf <- updated_vcf
filtered_vcf@gt <- filtered_vcf@gt[, c(TRUE, keep_samples)] # Keep the FORMAT column and subset samples

# Write the filtered VCF file to a new file
filtered_vcf_path <- "C:/Users/Desktop/DISS/squirrel_leprosy/Urban2024_leprosy_snp_table.haploid.biallelicOnly.filtered.vcf"
write.vcf(filtered_vcf, filtered_vcf_path)

# Verify the filtered VCF file
filtered_vcf_verify <- read.vcfR(filtered_vcf_path)
cat("Filtered sample IDs in the new VCF:", colnames(filtered_vcf_verify@gt)[-1], "\n")

# Define the path to the PLINK executable
plink_path <- "C:/Users/Desktop/DISS/plink_win64_20231211/plink.exe"

# Define the base name for the output files
filtered_out_path <- "C:/Users/Desktop/DISS/squirrel_leprosy/Urban2024_leprosy_snp_table.haploid.biallelicOnly.filtered.plink"

# Construct the PLINK command
plink_command_filtered <- sprintf(
    "%s --vcf %s --double-id --make-bed --out %s --chr-set -1",
    plink_path, filtered_vcf_path, filtered_out_path
)

# Print the command to check it (optional)
cat("Running command:", plink_command_filtered, "\n")

# Execute the PLINK command
system(plink_command_filtered)


# To make distance matrices for APE
plink_path <- "C:/Users/Desktop/DISS/plink_win64_20231211/plink.exe"
input_filtered <- "C:/Users/Desktop/DISS/squirrel_leprosy/Urban2024_leprosy_snp_table.haploid.biallelicOnly.filtered.plink"
output_dist_mat <- "C:/Users/Desktop/DISS/squirrel_leprosy/leprosy_onlyancient_distmat"


# Define the PLINK command
command_distmat <- paste0(plink_path, " --bfile ", input_filtered, " --distance square 1-ibs --out ", output_dist_mat)

# Run the command using system
system(command_distmat)


# Read the distance matrix into R
dist_matrix_ancient_leprosy <- as.matrix(read.table("C:/Users/Desktop/DISS/squirrel_leprosy/leprosy_onlyancient_distmat.mdist"))

# Find the largest value in the distance matrix
# Since there are values larger than 100 which bionj does not like
largest_value <- max(dist_matrix_ancient_leprosy, na.rm = TRUE)

print(paste("The largest value in the distance matrix is:", largest_value))

## NOT diving by largest value since it is all under 1 due to 1-IBS


# Read the .fam file into R
fam_for_samplenames <- read.table("C:/Users/Desktop/DISS/squirrel_leprosy/Urban2024_leprosy_snp_table.haploid.biallelicOnly.filtered.plink.fam")

# Extract the sample IDs
sample_ids <- fam_for_samplenames$V2

# Replace the row names and column names of the distance matrix
rownames(dist_matrix_ancient_leprosy) <- sample_ids
colnames(dist_matrix_ancient_leprosy) <- sample_ids


# Compute the BIONJ tree
tree <- bionj(dist_matrix_ancient_leprosy)

# Check for negative branch lengths and set them to zero
negative_branches <- tree$edge.length < 0
if (any(negative_branches)) {
    cat("Negative branch lengths found:\n")
    print(tree$edge.length[negative_branches])
    cat("Setting negative branch lengths to zero.\n")
    tree$edge.length[negative_branches] <- 0
} else {
    cat("No negative branch lengths found.\n")
}

# Print the BIONJ tree
plot(tree)

# Export the BIONJ tree
png(
    file = "C:/Users/Desktop/DISS/Human_closest/Leprosy_onlyancient_negativesas0_newfeb.png"
)

plot(tree, main = "Ancient Leprosy Phylogeny Tree")

dev.off()

write.tree(tree, file = "C:/Users/Desktop/DISS/Human_closest/ancient_leprosy_tree_newfeb.tree")
write.tree(tree, file = "C:/Users/Desktop/DISS/Human_closest/modern_leprosy_tree_new2907.txt")

###### Subset of only modern samples##################################

# Define the path to the updated VCF file
updated_vcf_path <- "C:/Users/Desktop/DISS/squirrel_leprosy/Urban2024_leprosy_snp_table.haploid.biallelicOnly.updated.vcf"

# Read the updated VCF file
updated_vcf <- read.vcfR(updated_vcf_path)

# Define the subset of samples to include
subset_samples_modern <- c(
    "TN", "2188-2007", "3208-2015", "3208-2007", "sample_2936", "sample_85054", "sample_97016", "1126-2007", "LRC-1A", "Airaku-3", "Amami", "ARLP_07", "ARLP_08", "ARLP_10", "ARLP_11", "ARLP_12", "ARLP_13", "ARLP_14", "ARLP_20", "ARLP_23", "ARLP_25", "ARLP_27", "ARLP_29", "ARLP_30", "ARLP_32", "ARLP_37", "ARLP_40", "ARLP_46", "ARLP_48", "ARLP_49", "ARLP_52", "ARLP_57", "ARLP_62", "ARLP_63", "ARLP_65", "ARLP_68", "ARLP_73", "ARLP_74", "Bn7-39", "Bn7-41", "Bn8-46", "Bn8-47", "Bn8-51", "Bn8-52", "BP", "Br14-1", "Br14-2", "Br14-3", "Br14-4", "Br14-5", "Br2016-14", "Br2016-15", "Br2016-16", "Br2016-17", "Br2016-18", "Br2016-19", "Br2016-20", "Br2016-21", "Br2016-24", "Br2016-26", "Br2016-27", "Br2016-45", "Br2016-46", "Br2016-47", "1126-2011", "2188-2014", "EGG", "Fio3", "Gu4-17", "Gu4-19L", "Gu5-23", "Izumi", "Kanazawa", "Kitasato", "Kusatsu-6", "Kyoto-1", "Ml2-10", "Ml6-50", "Ml6-55", "S15_92041", "Tsukuba-1", "US57", "Zensho-2", "Zensho-4", "Zensho-5", "Zensho-9", "Oku-4", "Ryukyu-2", "Kyoto-2", "Md1638", "Comore_3", "Md01046D", "Md02018_2C", "Md03031D", "Md04032", "Md04033", "Md05036_1D", "Md05036_2B", "Md05036_3BL", "Md06008", "Md07027B", "Md07037_1", "Md07037_2B", "Md07037_4", "Md08011", "Md09054B", "Md09055D", "Md10013", "Md10044", "Md10044C", "Md11060B", "MdEL_10", "MdEL_13", "MdEL_4", "MdEL_5", "MdEL_6", "MdEL_7", "MdEL_8", "MdEL_9", "MM23897", "MM24219", "MM29149", "MM29754", "Br4923", "1262-16", "2DDS", "Br1", "BrMM1", "BrMM2", "BrMM4", "BrMM5", "Indonesia-1", "Korea-3-2", "Thai-53", "Ml10-91", "Ml10-93", "Ml10-94", "Ml10-95", "Ml10-96", "Ml10-97", "Ml10-98", "Ml10-99", "ML2-5", "Ml9-79", "Ml9-80", "Ml9-81", "Ml9-82", "Ml9-83", "Ml9-84", "Ml9-86", "Ml9-87", "Ng12-33", "Ng13-32", "Ng13-33", "Ng14-35", "Ng15-36", "Ng15-37", "Ng16-38", "Ng17-39", "NHDP-55", "NHDP-63", "NHDP-98", "Pak", "S10_Ch-04", "S11_Inde_2", "S13_Ml-3-28", "S14_Ml-2-07", "S2_95034", "S9_96008", "Thai-237", "Thai-311", "Ye2-3", "Ye3s2", "Ye4-10", "Ye4-11", "Ye4-12", "Ye4-8", "B14096", "B14632", "B14668", "B191", "B204", "Br2016_68", "Br2016-141", "Comore_1", "Comore_2", "Congo1", "MGI-1", "MM137", "MM149", "MM29890", "Np14624", "Np14689", "Np14958", "Ph15", "Ph17", "RB022_2", "RB041", "RB048", "RB053", "RB065", "RB066", "RB067", "RB180", "RB181", "RB184", "RB187", "RB188", "VB-04", "VB-06", "VB-11", "VB-16", "VB-17", "VB-21", "VB-29", "VB-40", "VB-43", "VB-52", "VB-67", "VB-82", "RB001", "RB003", "RB063", "RB069", "RB073", "RB074", "RB182", "RB185", "RB186", "RN001", "RN022", "RN059", "RN084", "RN095", "RN165", "Bn7-37", "Bn7-38", "Bn9-59", "Bn9-66", "Bn9-67", "Bn10-71", "CI-7", "Ml2-7", "Ml3-17", "Ml6-52", "Ml11-101", "Ml11-102", "Ml11-103", "Ml12-109", "Ml12-110", "Ng19-42", "Ng22-45", "Ng27-58", "Ng27-60", "Ng29-64", "Sen1-1", "KFMC-1"
)

# Filter to keep only the subset of samples
keep_samples <- colnames(updated_vcf@gt)[-1] %in% subset_samples_modern
filtered_vcf <- updated_vcf
filtered_vcf@gt <- filtered_vcf@gt[, c(TRUE, keep_samples)] # Keep the FORMAT column and subset samples

# Write the filtered VCF file to a new file
filtered_vcf_path <- "C:/Users/Desktop/DISS/squirrel_leprosy/Urban2024_leprosy_snp_table.haploid.biallelicOnly.filteredmodern.vcf"
write.vcf(filtered_vcf, filtered_vcf_path)

# Verify the filtered VCF file
filtered_vcf_verify <- read.vcfR(filtered_vcf_path)
cat("Filtered sample IDs in the new VCF:", colnames(filtered_vcf_verify@gt)[-1], "\n")

# Define the path to the PLINK executable
plink_path <- "C:/Users/Desktop/DISS/plink_win64_20231211/plink.exe"

# Define the base name for the output files
filtered_out_path <- "C:/Users/Desktop/DISS/squirrel_leprosy/Urban2024_leprosy_snp_table.haploid.biallelicOnly.filteredmodern.plink"

# Construct the PLINK command
plink_command_filtered <- sprintf(
    "%s --vcf %s --double-id --make-bed --out %s --chr-set -1",
    plink_path, filtered_vcf_path, filtered_out_path
)

# Print the command to check it (optional)
cat("Running command:", plink_command_filtered, "\n")

# Execute the PLINK command
system(plink_command_filtered)

######################## creating modern only tree##################################

# To make distance matrices for APE
plink_path <- "C:/Users/Desktop/DISS/plink_win64_20231211/plink.exe"
input_filtered <- "C:/Users/Desktop/DISS/squirrel_leprosy/Urban2024_leprosy_snp_table.haploid.biallelicOnly.filteredmodern.plink"
output_dist_mat <- "C:/Users/Desktop/DISS/squirrel_leprosy/leprosy_onlymodern_distmat"


# Define the PLINK command
command_distmat <- paste0(plink_path, " --bfile ", input_filtered, " --distance square 1-ibs --out ", output_dist_mat)

# Run the command using system
system(command_distmat)


# Read the distance matrix into R
dist_matrix_modern_leprosy <- as.matrix(read.table("C:/Users/Desktop/DISS/squirrel_leprosy/leprosy_onlymodern_distmat.mdist"))

# Find the largest value in the distance matrix
# Since there are values larger than 100 which bionj does not like
largest_value <- max(dist_matrix_modern_leprosy, na.rm = TRUE)

print(paste("The largest value in the distance matrix is:", largest_value))

## NOT diving by largest value since it is all under 1 due to 1-IBS


# Read the .fam file into R
fam_for_samplenames <- read.table("C:/Users/Desktop/DISS/squirrel_leprosy/Urban2024_leprosy_snp_table.haploid.biallelicOnly.filteredmodern.plink.fam")

# Extract the sample IDs
sample_ids <- fam_for_samplenames$V2

# Replace the row names and column names of the distance matrix
rownames(dist_matrix_modern_leprosy) <- sample_ids
colnames(dist_matrix_modern_leprosy) <- sample_ids


# Compute the BIONJ tree
tree <- bionj(dist_matrix_modern_leprosy)

# Check for negative branch lengths and set them to zero
negative_branches <- tree$edge.length < 0
if (any(negative_branches)) {
    cat("Negative branch lengths found:\n")
    print(tree$edge.length[negative_branches])
    cat("Setting negative branch lengths to zero.\n")
    tree$edge.length[negative_branches] <- 0
} else {
    cat("No negative branch lengths found.\n")
}

# Print the BIONJ tree
plot(tree)

# Export the BIONJ tree
png(
    file = "C:/Users/Desktop/DISS/Human_closest/Leprosy_onlymodern_negativesas0_new2907.png",
    width = 2000, height = 4500
)

plot(tree, main = "Modern Leprosy Phylogeny Tree")

dev.off()

write.tree(tree, file = "C:/Users/Desktop/DISS/Human_closest/modern_leprosy_tree_new2907.txt")

read.tree("C:/Users/Desktop/DISS/Human_closest/modern_leprosy_tree_new2907.txt")
