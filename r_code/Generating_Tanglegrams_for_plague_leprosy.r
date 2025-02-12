setwd("C:/Users/Desktop/DISS/Human_closest")


library(ape)
library(phytools)
library(stringi)
library(parallel)
library(Rtapas)
library(data.table)
library(paco)
library(distory)
library(parallel)
library(xlsx)

######### CREATING TANGLEGRAMS for LEPROSY#############
# path to phylogenetic tree files
leprosy_tree_file <- "C:/Users/Desktop/DISS/Human_closest/ancient_leprosy_tree_new2907.txt"
human_tree_file <- "C:/Users/Desktop/DISS/Human_closest/human_tree_leprosy_ancientNEW.txt"

# Read the trees
leprosy_tree <- read.tree(leprosy_tree_file)
human_tree <- read.tree(human_tree_file)
leprosy_tree$tip.label
human_tree$tip.label

#### TESTING TANGLEGRAMS

# Create the host-symbiont association matrix (HS)

associations <- read.csv("C:/Users/Desktop/DISS/Human_closest/leprosy_closest_human_NEW.csv")

# Initialize the HS matrix
HS <- matrix(0,
    nrow = length(leprosy_tree$tip.label), ncol = length(human_tree$tip.label),
    dimnames = list(leprosy_tree$tip.label, human_tree$tip.label)
)


# Populate the HS matrix based on the CSV associations
for (i in 1:nrow(associations)) {
    pathogen_tip <- associations$sample_leprosy[i]
    human_tip <- associations$sample_human_leprosy[i]

    if (pathogen_tip %in% leprosy_tree$tip.label && human_tip %in% human_tree$tip.label) {
        HS[pathogen_tip, human_tip] <- 1
    } else {
        if (!pathogen_tip %in% leprosy_tree$tip.label) {
            print(paste("Pathogen tip not found:", pathogen_tip))
        }
        if (!human_tip %in% human_tree$tip.label) {
            print(paste("Human tip not found:", human_tip))
        }
    }
}

HS_matrix_df <- as.data.frame(HS)

HS_matrix_flipped <- t(HS_matrix_df)


############################### Compute max_cong using fqtab, human_tree, leprosy_tree, and HS
n <- 4
N <- 10000

fqtab <- max_cong(HS_matrix_flipped, human_tree, leprosy_tree, n, N,
    method = "paco", symmetric = TRUE, ei.correct = "sqrt.D",
    percentile = 0.01, res.fq = TRUE, strat = "parallel", cl = 4
)


## PAPER EXAMPLE TEST
data(nuc_cp)
N <- 10 # for the example, we recommend 1e+4 value
n <- 8
NPc <- max_cong(np_matrix, NUCtr, CPtr, n, N,
    method = "paco",
    symmetric = TRUE, ei.correct = "sqrt.D",
    percentile = 0.01, res.fq = FALSE,
    strat = "parallel", cl = 4
)
col <- c("darkorchid4", "gold")
tangle_gram(NUCtr, CPtr, np_matrix, NPc,
    colscale = "sequential",
    colgrad = col, nbreaks = 50, node.tag = TRUE
)
#### END OF EXAMPLE#####

# Define color gradient
colgrad <- c("#2C7BB6", "#FDCB58", "#D7191C")

# Call the tangle_gram function
treeH <- human_tree
treeS <- leprosy_tree

distmat_lep <- as.matrix(read.table("C:/Users/Desktop/DISS/squirrel_leprosy/leprosy_onlyancient_distmat.mdist"))
distmat_human <- as.matrix(read.table("C:/Users/Desktop/DISS/Human_closest/human_dist_mat_leprosymatch2907.mdist"))

plot(tangle_gram(treeH, treeS, HS_matrix_flipped, fqtab, colscale = "diverging", colgrad = colgrad, nbreaks = 50, node.tag = TRUE))


# Save a high-resolution PNG
png(
    file = "C:/Users/Desktop/DISS/Human_closest/leprosy_tanglegram_highres_SYMFALSE.png",
    width = 2000, height = 2000, res = 300
) # Increase resolution (300 dpi) and size

tangle_gram(treeH, treeS, HS_matrix_flipped, fqtab,
    colscale = "diverging", colgrad = colgrad,
    nbreaks = 50, node.tag = TRUE, cexpt = 1
) # Increase text size

dev.off()




ape::parafit(stats::cophenetic(human_tree), stats::cophenetic(leprosy_tree), HP = HS_matrix_flipped, nperm = 999, test.links = TRUE, seed = NULL, correction = "lingoes", silent = FALSE)
### Global test:  ParaFitGlobal = 0.0005631084 , p-value = 0.031 ( 999 permutations)/ 0.28 when we used N=10 but then increased permutations to 1000
### and got ParaFitGlobal = 0.0005345115 , p-value = 0.052 ( 999 permutations) another time p value = 0.58 and again 0.063
## after using symmetric =false in max_cong:
## ParaFitGlobal = 0.0005345115 , p-value = 0.056 ( 999 permutations)
## ParaFitGlobal = 0.0005345115 , p-value = 0.06 ( 999 permutations)



PACo_LFc <- max_cong(HS_matrix_flipped, treeH, treeS, n, N,
    method = "paco", symmetric = TRUE,
    ei.correct = "sqrt.D", percentile = 0.01, res.fq = TRUE,
    strat = "parallel", cl = 10
)
head(PACo_LFc)



write.xlsx(PACo_LFc, file = "C:/Users/Desktop/DISS/Human_closest/PACo_LFc2.xlsx")


################################ FOR PLAGUE################

# path to phylogenetic tree files
plague_tree_file <- "C:/Users/Desktop/DISS/yersinia_pestis/ancient_plague_Valtuena.txt"
human_tree_file <- "C:/Users/Desktop/DISS/Human_closest/human_tree_plague_ancientFEB.txt"

# Read the trees
plague_tree <- read.tree(plague_tree_file)
human_tree <- read.tree(human_tree_file)
plague_tree$tip.label
human_tree$tip.label

#### TESTING TANGLEGRAMS

# Create the host-symbiont association matrix (HS)

associations2 <- read.csv("C:/Users/Desktop/DISS/Human_closest/plague_closest_human_febfinal2.csv")

# Initialize the HS matrix
HS <- matrix(0,
    nrow = length(plague_tree$tip.label), ncol = length(human_tree$tip.label),
    dimnames = list(plague_tree$tip.label, human_tree$tip.label)
)


# Populate the HS matrix based on the CSV associations
for (i in 1:nrow(associations2)) {
    pathogen_tip <- associations2$sample_plague[i]
    human_tip <- associations2$sample_human_plague[i]

    if (pathogen_tip %in% plague_tree$tip.label && human_tip %in% human_tree$tip.label) {
        HS[pathogen_tip, human_tip] <- 1
    } else {
        if (!pathogen_tip %in% plague_tree$tip.label) {
            print(paste("Pathogen tip not found:", pathogen_tip))
        }
        if (!human_tip %in% human_tree$tip.label) {
            print(paste("Human tip not found:", human_tip))
        }
    }
}

HS_matrix_df <- as.data.frame(HS)

HS_matrix_flipped <- t(HS_matrix_df)


#################################### Compute max_cong using fqtab, human_tree, leprosy_tree, and HS
n <- 4
N <- 10000

fqtab <- max_cong(HS_matrix_flipped, human_tree, plague_tree, n, N,
    method = "paco", symmetric = TRUE, ei.correct = "sqrt.D",
    percentile = 0.01, res.fq = TRUE, strat = "parallel", cl = 4
)


# Define color gradient
colgrad <- c("#54278f", "#fdae61", "#006837")

# Call the tangle_gram function
treeH <- human_tree
treeS <- plague_tree

distmat_plague <- as.matrix(read.table("C:/Users/Desktop/DISS/yersinia_pestis/ancient_plague_Valtuena_distmat.mdist"))
distmat_human <- as.matrix(read.table("C:/Users/Desktop/DISS/Human_closest/human_dist_mat_plaguematchFEB_clean.mdist"))

plot(tangle_gram(treeH, treeS, HS_matrix_flipped, fqtab, colscale = "diverging", colgrad = colgrad, nbreaks = 50, node.tag = TRUE))

# Plot for plague samples with 3 colours
png(
    file = "C:/Users/Desktop/DISS/Human_closest/plague_tanglegram.png",
    width = 750, height = 750
)
tangle_gram(treeH, treeS, HS_matrix_flipped, fqtab, colscale = "diverging", colgrad = colgrad, nbreaks = 50, node.tag = TRUE, cexpt = 1)
dev.off()

# Identify plague samples with at least one connection
connected_plague_samples <- colnames(HS_matrix_flipped)[colSums(HS_matrix_flipped) > 0]

# Prune the plague tree to keep only the connected tips
treeS <- drop.tip(plague_tree, setdiff(plague_tree$tip.label, connected_plague_samples))

# Update the association matrix to match the pruned tree
HS_matrix_flipped <- HS_matrix_flipped[, connected_plague_samples, drop = FALSE]

# Recompute fqtab to reflect updated HS matrix
fqtab <- max_cong(HS_matrix_flipped, treeH, treeS, n, N,
    method = "paco", symmetric = TRUE, ei.correct = "sqrt.D",
    percentile = 0.01, res.fq = TRUE, strat = "parallel", cl = 4
)

# Plot the updated tanglegram
plot(tangle_gram(treeH, treeS, HS_matrix_flipped, fqtab,
    colscale = "diverging",
    colgrad = colgrad, nbreaks = 50, node.tag = TRUE
))


# Save a high-resolution PNG
png(
    file = "C:/Users/Desktop/DISS/Human_closest/plague_tanglegram_highres.png",
    width = 2000, height = 2000, res = 300
) # Increase resolution (300 dpi) and size

tangle_gram(treeH, treeS, HS_matrix_flipped, fqtab,
    colscale = "diverging", colgrad = colgrad,
    nbreaks = 50, node.tag = TRUE, cexpt = 1
) # Increase text size

dev.off()




#### sorting missing/mismatch values before parafit
dim(stats::cophenetic(human_tree)) # Should be (n.hosts x n.hosts)
dim(stats::cophenetic(plague_tree)) # Should be (n.parasites x n.parasites)
dim(HS_matrix_flipped) # Should be (n.hosts x n.parasites)

# Keep only matching hosts and parasites
common_hosts <- intersect(rownames(HS_matrix_flipped), rownames(stats::cophenetic(human_tree)))
common_parasites <- intersect(colnames(HS_matrix_flipped), rownames(stats::cophenetic(plague_tree)))

# Subset matrices to only include matching samples
HS_matrix_filtered <- HS_matrix_flipped[common_hosts, common_parasites, drop = FALSE]
human_tree_filtered <- drop.tip(human_tree, setdiff(human_tree$tip.label, common_hosts))
plague_tree_filtered <- drop.tip(plague_tree, setdiff(plague_tree$tip.label, common_parasites))

# Re-run parafit
ape::parafit(
    stats::cophenetic(human_tree_filtered),
    stats::cophenetic(plague_tree_filtered),
    HP = HS_matrix_filtered,
    nperm = 999,
    test.links = FALSE,
    correction = "lingoes",
    silent = FALSE
)

#### Global test:  ParaFitGlobal = 0.008052437 , p-value = 0.002 ( 999 permutations)



PACo_LFc <- max_cong(HS_matrix_flipped, treeH, treeS, n, N,
    method = "paco", symmetric = TRUE,
    ei.correct = "sqrt.D", percentile = 0.01, res.fq = TRUE,
    strat = "parallel", cl = 10
)
head(PACo_LFc)



write.xlsx(PACo_LFc, file = "C:/Users/Desktop/DISS/Human_closest/PACo_LFc_plague.xlsx")
