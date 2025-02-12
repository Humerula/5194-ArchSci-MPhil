setwd("C:/Users/Desktop/DISS/yersinia_pestis")

library(ape)


## on bash cp /path/to/LNBAdataset2021_gok2_RV2039_SNPExcluded_genotyped.raxml.bestTree /mnt/c/path/to/destination/

## Reading tree from Valtuena et al 2022

# Define the path to the Newick tree file
plague_tree_file <- "C:/Users/Desktop/DISS/yersinia_pestis/raxml/LNBAdataset2021_gok2_RV2039_SNPExcluded_genotyped--threads4.raxml.bestTree"

# Read the tree
plague_tree <- read.tree(plague_tree_file)
print(plague_tree$tip.label)



# Export the tree
png(
    file = "C:/Users/Desktop/DISS/yersinia_pestis/LNBA_plague_tree_ver9.png",
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
    file = "C:/Users/Desktop/DISS/yersinia_pestis/LNBA_plague_tree_ver4_collapsed_9.png",
    width = 2000, height = 4500
)
# Plot the tree
plot(plague_tree2, main = "Plague Tree from LNBA raxML collapsed")

dev.off()

########################## CREATING TREE FOR ANCIENT SAMPLES TO USE FOR COPHYLO###############################################################

## check if all ancient subset names match with the tree tip labels

# creating a subset of samples from plague_tree to create distance matrix of relevant samples
### (I5884 and I2470 skipped as missing from Valtuena et al 2022 raxml)

plague_tree_ancient_trim <- keep.tip(plague_tree, tip = c(
    "VLI092", "KNK001_c", "GRS004", "HOP001", "HOP004",
    "KLZ001_c", "CHC004", "KLE031", "MIB054", "KLE048",
    "OOH003_c", "ARS007_c", "KZL002", "VEL003_c", "GRH001_c",
    "OBS124", "OBS107", "OBS116", "OBS137", "OBS110",
    "BED028.A0102_8103_PE_SE", "NMS002.A0101.combined", "BED030.A0102_8127_PE_SE", "BED024.A0102_8052_PE_SE", "STN020.A0101_SE",
    "STN021.A0101_SE", "STA001.A0101_PE_SE", "LAI009.A0101_RT88_PE_SE", "BED034.A0102_8198_PE_SE", "STN014.A0101_SE",
    "BRA001.A0101_PE_SE", "STN008.A0101_SE", "STN007.A0101_SE", "MAN008.B0101_SE", "London_EastSmithfield_8124_8291_11972",
    "Gok2", "RV2039", "RISE509", "RISE505", "Barcelona",
    "GEN72", "Gyvakarai1", "PST006.A_Post6", "1343UNTAL85", "RK1001.C",
    "KunilaII", "MIK005.A_RT5", "GZL001.A0101_02.YP2.1", "GZL002.A0101_02.YP2.1", "Bolgar",
    "Altenerding2018", "STN019.A0101_SE", "STN002.A0101_SE", "2.MED1_M-549", "EDI001.A_shotgun"
))

print(plague_tree_ancient_trim$tip.label)
length(plague_tree_ancient_trim$tip.label)


# trying to get distance matrix from the uncollapsed tree for finding close matching samples
distmat_a <- cophenetic(plague_tree_ancient_trim)

View(distmat_a)

write.csv(distmat_a, file = "C:/Users/Desktop/DISS/yersinia_pestis/ancient_plague_Valtuena_distmat.mdist")

# Exporting to make ancient tree only using FigTree
write.tree(plague_tree_ancient_trim, file = "C:/Users/Desktop/DISS/yersinia_pestis/ancient_plague_Valtuena.tree")
write.tree(plague_tree_ancient_trim, file = "C:/Users/Desktop/DISS/yersinia_pestis/ancient_plague_Valtuena.txt")
