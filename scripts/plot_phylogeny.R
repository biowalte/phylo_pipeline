#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(phytools)
  library(ggplot2)
  library(dplyr)
})

# ------------------------------------------------------------
# Arguments
# ------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript plot_phylogeny.R <tree.nwk | tree.con.tree>")
}

tree_file <- args[1]

# ------------------------------------------------------------
# File extension validation
# ------------------------------------------------------------
if (!grepl("\\.(nwk|tree|con\\.tree|contree|treefile)$", tree_file)) {
  stop("File must be a valid Newick tree (.nwk, .tree, .contree, .treefile)")
}

# ------------------------------------------------------------
# Read tree (Newick in all cases)
# ------------------------------------------------------------
tree <- read.tree(tree_file)

# ------------------------------------------------------------
# Detect node support type
# ------------------------------------------------------------
support_values <- suppressWarnings(as.numeric(tree$node.label))

support_type <- "Support"

if (all(!is.na(support_values))) {
  if (max(support_values, na.rm = TRUE) <= 1) {
    support_type <- "Posterior Probability"
  } else if (max(support_values, na.rm = TRUE) <= 100) {
    support_type <- "Bootstrap"
  }
}

# ------------------------------------------------------------
# Output directory
# ------------------------------------------------------------
dir.create("plots", showWarnings = FALSE)

# ============================================================
# 1) Simple tree
# ============================================================

p_simple <- ggtree(tree) +
  geom_tiplab(size = 3) +
  ggtitle("Simple phylogenetic tree")

ggsave("plots/tree_simple.png", p_simple, width = 8, height = 6)

# ============================================================
# 2) Tree colored by depth
# ============================================================

p_colored <- ggtree(tree, aes(color = y)) +
  geom_tiplab(size = 3) +
  scale_color_viridis_c() +
  ggtitle("Phylogenetic tree colored by depth")

ggsave("plots/tree_colored.png", p_colored, width = 8, height = 6)

# ============================================================
# 3) Tree with node support values (bootstrap or PP)
# ============================================================

p_support <- ggtree(tree) +
  geom_tiplab(size = 3) +
  geom_nodelab(
    aes(label = label),
    hjust = -0.3,
    size = 3,
    color = "firebrick"
  ) +
  ggtitle(paste("Phylogenetic tree with", support_type))

ggsave("plots/tree_support.png", p_support, width = 8, height = 6)

# ============================================================
# 4) Circular tree
# ============================================================

p_circular <- ggtree(tree, layout = "circular") +
  geom_tiplab(size = 2) +
  ggtitle("Circular phylogenetic tree")

ggsave("plots/tree_circular.png", p_circular, width = 8, height = 8)

# ============================================================
# 5) Automatic clade highlighting
# ============================================================

internal_nodes <- (Ntip(tree) + 1):(Ntip(tree) + min(3, tree$Nnode))

p_clades <- ggtree(tree) +
  geom_tiplab(size = 3)

for (node in internal_nodes) {
  p_clades <- p_clades +
    geom_hilight(
      node = node,
      fill = sample(colors(), 1),
      alpha = 0.3
    )
}

p_clades <- p_clades +
  ggtitle("Automatic clade highlighting")

ggsave("plots/tree_clades.png", p_clades, width = 8, height = 6)

# ============================================================
# 6) phytools (base plotting)
# ============================================================

png("plots/tree_phytools.png", width = 1000, height = 800)
plotTree(tree, ftype = "i", fsize = 0.8)
nodelabels(tree$node.label, frame = "none", cex = 0.7)
title(paste("Phylogenetic tree (phytools) -", support_type))
dev.off()

cat(" Tree file read:", tree_file, "\n")
cat(" Detected support type:", support_type, "\n")
cat(" Plots generated in ./plots\n")
