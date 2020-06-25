# Need to clean up package install as not all are used!
# Using checkM tree, chemical data and taxa abundance to produce oridnations.

list.of.packages <- c("vegan", "MASS", "ggplot2", "phyloseq", "dplyr", "ape")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")

library("phyloseq")
packageVersion("phyloseq")
## [1] 1.26.0’
library("vegan")
library(MASS)
packageVersion("vegan")
## [1] ‘2.5.3’
library("ggplot2")
packageVersion("ggplot2")
## [1] ‘3.1.0’
library("grid")
## [1] ‘3.5.1’
library(dplyr)
library(plyr)
library(ape)
#library(ggtree)
#library(tidytree)
detach("package:ggtree", unload=TRUE)
detach("package:treeio", unload=TRUE)


# tree
checkm_tree <- read.tree("tree_qa.tre")
checkm_tree$tip.label <- gsub("-contigs", "", checkm_tree$tip.label)
length(grep("Bin", checkm_tree$tip.label))
# trimming the the checkm tree
is.rooted(checkm_tree)
new_checkm_tree <- ape::drop.tip(checkm_tree, tip=checkm_tree$tip.label[grep("IMG", checkm_tree$tip.label)], trim.internal = TRUE)
bad_bins <- c("Bin_39_8", "Bin_10_36", "Bin_17_9", "Bin_22_3", "Bin_32_1", "Bin_10_16", "Bin_7_11", "Bin_37_1")
matches <-  grep(paste(bad_bins,collapse="|"), new_checkm_tree$tip.label, value=TRUE)
matches <- matches[-3]
cat(paste(shQuote(matches, type="cmd"), collapse=", "))
new_checkm_tree <- ape::drop.tip(new_checkm_tree, tip = matches, trim.internal = TRUE)
length(new_checkm_tree$tip.label)

#reading in taxonmy file
new_RP_MAG_taxonomy <- read.table("COMBINED_bac_arc_bin_phylogeny.csv", sep=",", header=F, na.strings=c("","NA"))
colnames(new_RP_MAG_taxonomy) <- c("label", "taxa")
rownames(new_RP_MAG_taxonomy) <- new_RP_MAG_taxonomy$label
new_RP_MAG_taxonomy$Species <- new_RP_MAG_taxonomy$taxa

#reading otu table 
otu_abund <- read.table("mapped_reads_log10.csv", sep=",", header = T)
rownames(otu_abund) <- otu_abund[,1]
otu_abund <- otu_abund[,-1]
colnames(otu_abund)
sapply(otu_abund, class)
colSums(otu_abund)
otu_abund_edit <- otu_abund
otu_abund_edit[] <- lapply(otu_abund, function(x) ifelse(x>0.5, 0, x))

# not norm here means that the reads have been normalized by sequencing depth
# but not normalized on a log scale as above
otu_abund_not_norm <- read.table("mapped_reads_clean.csv", sep=",", header = T)
rownames(otu_abund_not_norm) <- otu_abund_not_norm[,1]
otu_abund_not_norm <- otu_abund_not_norm[,c(-1,-10,-11)]
sapply(otu_abund_not_norm, class)
colSums(otu_abund_not_norm)
colnames(otu_abund_not_norm)

otu_abund_not_norm_edit <- otu_abund_not_norm
otu_abund_not_norm_edit[] <- lapply(otu_abund_not_norm, function(x) ifelse(x<0.31, 0, x))


#reading in chemistry 
sample_chem <- read.table("chem_with_depth_final.csv", sep=",", header = T)
rownames(sample_chem) <- sample_chem[,1]
sample_chem <- sample_chem[,-1]
sample_chem <- sample_chem[,-12]
sample_chem$WaterDepth <- as.numeric(sample_chem$WaterDepth)
sapply(sample_chem, class)
# remove TDP
sample_chem <- sample_chem[,-8]
sample_chem$ID <- rownames(sample_chem)

fake_sample_chem <- sample_chem
H.soil <- c(1:11)
L.soil <- c(1:11)
C.soil <- c(1:11)
fake_sample_chem <- do.call("rbind", list(sample_chem, H.soil, L.soil, C.soil))
rownames(fake_sample_chem) <- c("H1", "H2", "L1", "L2", "C", "H.soil", "L.soil", "C.soil")
fake_sample_chem$SampleType <- c("H", "H", "L", "L", "C", "Soil", "Soil", "Soil")
fake_sample_chem$ID <- rownames(fake_sample_chem)

##### creating phylo objects
phylo_object <- phyloseq(otu_table(otu_abund_not_norm, taxa_are_rows=TRUE), tax_table(as.matrix(new_RP_MAG_taxonomy)), sample_data(sample_chem), phy_tree(new_checkm_tree))
phylo_bars_all <- phyloseq(otu_table(otu_abund_not_norm, taxa_are_rows=TRUE), tax_table(as.matrix(new_RP_MAG_taxonomy)), phy_tree(new_checkm_tree))

# Create ordination
# phyloseq 
abundance_bars <- plot_bar(phylo_object, fill="taxa")
small_tree_palette <- c("Gammaproteobacteria" = "#2d7eff", 
                          "Patescibacteria" = "#66a561", 
                          "Myxococcota" =   "#682f3f",
                          "Acidobacteriota" =   "#a32fdc", 
                          "Nitrospirota" =  "#6189a5",
                          "Planctomycetota" = "#7092ad",
                          "Eisenbacteria" = "#7a5e1d",
                          "Gemmatimonadota" = "#681c54",
                          "Zixibacteria" = "#3254ff",
                          "Bacteroidota" =   "#d04ee6",
                          "KSB1" = "#f76d3b",
                          "Actinobacteriota" = "#399840",
                          "Chloroflexota" = "#8c6df2",
                          "Armatimonadota" = "#104c00",
                          "Cyanobacteriota" = "#f2c232",
                          "Spirochaetota" = "#007265",
                          "Omnitrophota" = "#837f66",
                          "Lindowbacteria" = "#a5619e",
                          "Methylomirabilota" =   "#ccbb00",
                          "UBP1" =  "#00626d",
                          "UBP10" =  "#35a0dd",
                          "Alphaproteobacteria" = "#74912d",
                          "Verrucomicrobiota" = "#5d97a0", 
                          "Crenarchaeota" = "#77496d", 
                          "Nanoarchaeota" = "#5d97a0", 
                        "Samples" = "#000000")
abundance_bars + scale_color_manual(values=small_tree_palette) + scale_fill_manual(values=small_tree_palette) + theme_minimal()


abundance_bars_all <- plot_bar(phylo_bars_all, fill="taxa")
abundance_bars_all + theme_bw() + scale_color_manual(values=small_tree_palette) + 
  scale_fill_manual(values=small_tree_palette) 

# transforming data prior to ordination 
sed_ord = transform_sample_counts(phylo_object, function(x) x / sum(x) )
#sed_ord_filter = filter_taxa(sed_ord, function(x) mean(x) > 1e-5, TRUE)
#sed_ord = transform_sample_counts(phylo_object, function(x) 1E6 * x/sum(x))
phylum.sum = tapply(taxa_sums(sed_ord), tax_table(sed_ord)[, "taxa"], sum, na.rm=TRUE)
#top15phyla = names(sort(phylum.sum, TRUE))[1:15]
#sed_ord = prune_taxa((tax_table(sed_ord)[, "taxa"] %in% top15phyla), sed_ord)  
phylum.sum = tapply(taxa_sums(sed_ord), tax_table(sed_ord)[, "taxa"], sum, na.rm=TRUE)
top20phyla = names(sort(phylum.sum, TRUE))[1:20]
phylo_object2_top20 = prune_taxa((tax_table(sed_ord)[, "taxa"] %in% top20phyla), sed_ord)

#plotting ordination
sed.ord <- ordinate(sed_ord, "NMDS", "bray")
p1 = plot_ordination(sed_ord, sed.ord, type="taxa", color="taxa", title="taxa")
p1 + scale_color_manual(values=small_tree_palette) + scale_fill_manual(values=small_tree_palette)
#plotting top 20 ordination
sed.ord <- ordinate(phylo_object2_top20, "NMDS", "bray")
p1 = plot_ordination(sed_ord, sed.ord, type="taxa", color="taxa", title="taxa")
p1 + scale_color_manual(values=small_tree_palette) + scale_fill_manual(values=small_tree_palette)

# One plot top 20
p4 <- plot_ordination(phylo_object2_top20, sed.ord, type="biplot", color="taxa", label ="ID") 
p4 <- p4 + theme_bw() + scale_color_manual(values=small_tree_palette) + scale_fill_manual(values=small_tree_palette)
#pdf("nmds_sed_sites_taxa_biplot.pdf", height=6, width = 8)
p4
#dev.off()
ef_RP_MAG <- envfit(sed.ord,sample_data(phylo_object2_top20), permu=30000)
# > ef_RP_MAG
# 
# ***VECTORS
# 
# NMDS1    NMDS2     r2  Pr(>r)  
# Oxygen      0.98707  0.16031 0.7144 0.19167  
# pH         -0.37535  0.92688 0.7298 0.28333  
# Redox       0.96818  0.25026 0.9915 0.03333 *
#   OC         -0.08790 -0.99613 0.5766 0.40833  
# CaCO3      -0.50474  0.86327 0.7494 0.33333  
# NH3        -0.83099  0.55629 0.9743 0.02500 *
#   NO2_.NO3    0.92627  0.37686 0.9847 0.02500 *
#   Cl         -0.09255  0.99571 0.4154 0.56667  
# SO4        -0.98617 -0.16575 0.8762 0.10000 .
# WaterDepth -0.08095 -0.99672 0.5222 0.37500  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Permutation: free
# Number of permutations: 119
# 
# ***FACTORS:
#   
#   Centroids:
#   NMDS1   NMDS2
# IDC   0.3911 -0.0560
# IDH1 -0.5049  0.2321
# IDH2 -0.3772 -0.2833
# IDL1  0.3078  0.1012
# IDL2  0.1832  0.0059
# 
# Goodness of fit:
#   r2 Pr(>r)
# ID  1      1
# Permutation: free
# Number of permutations: 119

ef_RP_MAG.df <- as.data.frame(ef_RP_MAG$vectors$arrows*sqrt(ef_RP_MAG$vectors$r))
ef_RP_MAG.df$species <- rownames(ef_RP_MAG.df)

# Add labels, make a data.frame
#arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
# Define the arrow aesthetic mapping
arrow_map = aes(xend = NMDS1 * 0.5, yend = NMDS2 * 0.5, x = 0, y = 0, shape = NULL, color = NULL, 
                label = species)
label_map = aes(x = 0.65 * NMDS1, y = 0.65 * NMDS2, shape = NULL, color = NULL, 
                label = species)
# Make a new graphic
arrowhead = arrow(length = unit(0.02, "npc"))
p4_chem = p4 + geom_segment(arrow_map, size = 0.5, data = ef_RP_MAG.df, color = "gray", 
                       arrow = arrowhead) + geom_text(label_map, size = 2, data = ef_RP_MAG.df)
p4_chem

plot(sed.ord)
plot(ef_RP_MAG)
plot(ef_RP_MAG, p.max = 0.05, col = "red")





