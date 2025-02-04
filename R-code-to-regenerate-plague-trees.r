setwd("C:/Users/Desktop/DISS/yersinia_pestis")

library(ape)


## on bash cp /path/to/LNBAdataset2021_gok2_RV2039_SNPExcluded_genotyped.raxml.bestTree /mnt/c/path/to/destination/


# Define the path to the Newick tree file
plague_tree_file <- "C:/Users/Desktop/DISS/yersinia_pestis/raxml/LNBAdataset2021_gok2_RV2039_SNPExcluded_genotyped--threads4.raxml.bestTree"

# Read the tree
plague_tree <- read.tree(plague_tree_file)


# Export the tree
png(
    file = "C:/Users/Desktop/DISS/yersinia_pestis/LNBA_plague_tree_ver3.png",
    width = 2000, height = 4500
)
# Plot the tree
plot(plague_tree, main = "Plague Tree from LNBA raxML")

dev.off()

# Define the path to the Newick tree file
plague_tree_file2 <- "C:/Users/Desktop/DISS/yersinia_pestis/raxml/LNBAdataset2021_gok2_RV2039_SNPExcluded_genotyped--threads4.raxml.bestTreeCollapsed"

# Read the tree
plague_tree2 <- read.tree(plague_tree_file2)


# Export the tree
png(
    file = "C:/Users/Desktop/DISS/yersinia_pestis/LNBA_plague_tree_ver4_collapsed.png",
    width = 2000, height = 4500
)
# Plot the tree
plot(plague_tree2, main = "Plague Tree from LNBA raxML collapsed")

dev.off()
