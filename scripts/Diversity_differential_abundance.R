
#### Script for assignment 3 ####

## Installing and loading packages ##

#BiocManager::install("phyloseq")
library(phyloseq)

#BiocManager::install("biomformat")
library(biomformat)

#BiocManager::install("ANCOMBC")
library(ANCOMBC)
#BiocManager::install("microbiome")
library(microbiome)

library(vegan)
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)

## Read the BIOM file created with kraken2 ##
# Moved to data folder after

biom_data <- read_biom("table.biom")

# Convert to phyloseq object

physeq <- import_biom(biom_data)

# Creating OTU table and species rarefaction curve

otu_table <- as.data.frame(t(otu_table(physeq)))

# Just to check 
#rare_curve <- rarecurve(otu_table, step = 1000)

row.names(otu_table) ##?

# Calculating at relative abundance

physeq_rel <- transform_sample_counts(physeq, function(x)
  x / sum(x))
physeq_phy <- tax_glom(physeq_rel, taxrank = "Rank7") # Rank7 is species

# Add in metadata

sample_names <- c("SRR8146935_bracken_species", "SRR8146938_bracken_species", "SRR8146951_bracken_species", "SRR8146954_bracken_species", "SRR8146936_bracken_species", "SRR8146944_bracken_species", "SRR8146952_bracken_species", "SRR8146956_bracken_species")
diet <- as.factor(c("omnivore", "omnivore", "vegan", "vegan", "omnivore", "vegan", "vegan", "omnivore"))

metadata <- as.data.frame(diet)
rownames(metadata) <- sample_names
metadata

sampledata <- sample_data(metadata)

# Add it to phyloseq object

physeq_labelled <- merge_phyloseq(physeq_phy, sampledata)

# Turn into dataframe for plotting

df <- psmelt(physeq_labelled)

# Plotting only the top 20 most abundant species, add in diet information

top_20 <- df %>%
  slice_max(Abundance, n = 20) %>%
  mutate(Rank6 = sub("^g__", "", Rank6), Rank7 = sub("^s__", "", Rank7)) %>%
  mutate(species_names = paste0(Rank6, " ", Rank7), sample_names = gsub("_bracken_species$", "", Sample)) %>%
  group_by(sample_names) %>%
  mutate(species_names = gsub("^(\\w+)(\\s|$)", "\\1 ", species_names))

# See which taxon has the highest relative abundance to write out to table
top_20_persample <- top_20 %>%
  group_by(Sample) %>%
  slice_max(Abundance) %>%
  ungroup() %>%
  select(sample_names, Abundance, species_names, diet)

write.csv(top_20_persample, "../data/top_rel_abundance_persample.csv", row.names = F)

# Plot  reorder(sample_names, diet)
top_20$Sample <- factor(top_20$Sample, levels = unique(top_20$Sample[order(top_20$diet)]))

rel_abundance_plot <- ggplot(top_20, aes(x=factor(sample_names, levels = c("SRR8146935", "SRR8146936", "SRR8146938", "SRR8146956", "SRR8146944", "SRR8146951", "SRR8146952", "SRR8146954")), y = Abundance, fill = species_names)) + geom_bar(stat = "identity", position = "stack") + labs(y = "Relative Abundance", x = "Diet for each sample", fill="Species names")+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 4)) + scale_x_discrete(labels = c("SRR8146935"="omnivore", "SRR8146936"="omnivore", "SRR8146938"="omnivore", "SRR8146956"="omnivore", "SRR8146944"="vegan", "SRR8146951"="vegan", "SRR8146952"="vegan", "SRR8146954"="vegan"))

ggsave("../figures/rel_abundance_plot.png", rel_abundance_plot)

## Alpha diversity: I want to display Chao1, Berger-Parker, and Shannon
# I add the metadata to the raw counts physeq

physeq_raw <- merge_phyloseq(physeq, sampledata) 

#plot_richness(physeq_raw, color = "diet", measures = c("Chao1", "Shannon")) + scale_x_discrete(labels = c("SRR8146935", "SRR8146936", "SRR8146938", "SRR8146956", "SRR8146944", "SRR8146951", "SRR8146952", "SRR8146954")) + labs(x = "Sample Name", colour="Diet")

# Display in dataframe
alpha_div <- estimate_richness(physeq, measures = c("Chao1", "Shannon"))
alpha_div <- as.data.frame(alpha_div)

# Berger-Parker Index
bpi <- dominance(physeq, index = "DBP", rank = 7, relative = T, aggregate = T)
bpi <- as.data.frame(bpi)

# Merge into table to write out

alpha_diversity_indices <- merge(alpha_div, bpi, by = "row.names", all = T)
alpha_diversity_indices <- merge(alpha_diversity_indices, metadata, by.x = "Row.names", by.y = "row.names", all =T)
  
write.csv(alpha_diversity_indices, file = "../data/alpha_diversity_indices.csv", row.names = T)

# Plotting

chao1 <- ggplot(alpha_diversity_indices, aes(x=Row.names, y=Chao1, colour = diet)) + geom_point() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=6), legend.position = "none")+ scale_x_discrete(labels = c("SRR8146935", "SRR8146936", "SRR8146938", "SRR8146956", "SRR8146944", "SRR8146951", "SRR8146952", "SRR8146954")) + labs(x = "")
shannon <- ggplot(alpha_diversity_indices, aes(x=Row.names, y=Shannon, colour = diet)) + geom_point()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6), legend.position = "top")+ scale_x_discrete(labels = c("SRR8146935", "SRR8146936", "SRR8146938", "SRR8146956", "SRR8146944", "SRR8146951", "SRR8146952", "SRR8146954")) + labs(x = "Sample Name")
bp <- ggplot(alpha_diversity_indices, aes(x=Row.names, y=dbp, colour = diet)) + geom_point()+ scale_x_discrete(labels = c("SRR8146935", "SRR8146936", "SRR8146938", "SRR8146956", "SRR8146944", "SRR8146951", "SRR8146952", "SRR8146954")) + labs(x = "", y = "Berger-Parker")+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size= 6), legend.position = "none")


a_diversity_plot <- chao1 + shannon + bp # Side-by-side

ggsave("../figures/alpha_diversity_plot.png", a_diversity_plot)

### Beta diversity ###

# PCoA with Bray-Curtis
ord.pcoa.bray <- ordinate(physeq_labelled, method="PCoA", distance="bray")

# Plotting
braycurtis <- plot_ordination(physeq_labelled, ord.pcoa.bray, color = "diet") +
  geom_point(aes(color = diet), size = 4, inherit.aes = TRUE) + labs(x="PC1=46.3%", y="PC2=22.9%", color="Diet")

ggsave("../figures/braycurtis_plot.png", braycurtis)

# NMDS with Bray-Curtis and plotting
ord.nmds.bray <- ordinate(physeq_labelled, method="NMDS", distance="bray")

plot_ordination(physeq_labelled, ord.nmds.bray, color="diet") + geom_point(aes(color = diet), size = 4, inherit.aes = TRUE) + labs(color="Diet")

# PcOA with Jaccard distance and plotting
ord.pcoa.jaccard <- ordinate(physeq_labelled, method="PCoA", distance="jaccard")

jaccard <- plot_ordination(physeq_labelled, ord.pcoa.jaccard, color="diet") + geom_point(aes(color = diet), size = 4, inherit.aes = TRUE) + labs(x="PC1=35.4%", y="PC2=20.5%", color="Diet")

ggsave("../figures/jaccard_plot.png", jaccard)

## Permanova, not significant

adonis2(phyloseq::distance(physeq_labelled, method = "bray") ~ diet, data = metadata)

## Differential abundance analysis

#sampledata <- sample_data(metadata)

## Add to phyloseq object which has raw counts and isn't agglomerated

physeq_ancom <- merge_phyloseq(physeq, sampledata)

# I use the Benjamini Hochberg correction

ancombc_out <- ancombc2(data = physeq_ancom, tax_level = "Rank7", fix_formula = "diet", rand_formula = NULL, p_adj_method = "BH", pseudo_sens = TRUE, prv_cut = 0, lib_cut = 0, s0_perc = 0.05, group = "diet", struc_zero = TRUE, neg_lb = TRUE)

# Write out all data in the results section
write.csv(ancombc_out$res, file = "../data/DA_results_all.csv", row.names = F)

# Check which taxa are significantly different

ancombc_sig <- subset(ancombc_out$res, q_dietvegan < 0.05)

# There are none! Maybe some approach 0.05?

subset(ancombc_out$res, q_dietvegan < 0.3) # This is the result with the lowest q-value, really nothing significant

# Looking at which taxa have 0 between samples

zero_taxa <- subset(ancombc_out$zero_ind, `structural_zero (diet = vegan)` != `structural_zero (diet = omnivore)`)
head(zero_taxa)

# Even though abundance is not significantly different, I will look at the taxa with the largest log-fold change in the LFCvegan column (the top 20)

top20_highest_lfc <- ancombc_out$res %>% 
  slice_max(abs(lfc_dietvegan), n=20)

# Write out 
top20_highest_lfc
write.csv(top20_highest_lfc, file = "../data/top20_DAresults_LFC.csv", row.names = F)
 
# Adding a column for naming and plotting

tax_matching <- as.data.frame(physeq_ancom@tax_table)

merged_data <- merge(x = tax_matching, y = top20_highest_lfc, by.x="Rank7", by.y="taxon", all.y = T)

top20_for_plotting <- merged_data %>%
  select(-c(Rank1, Rank2, Rank3, Rank4, Rank5)) %>%
  mutate(Rank7 = sub(".*_g__", "", Rank7)) %>%
  mutate(Rank6 = sub("^g__", "", Rank6), Rank7 = sub("^s__", "", Rank7)) %>%
  mutate(Rank6 = replace_na(Rank6, "")) %>%
  mutate(Rank7 = sub("_s__", " ", Rank7)) %>%
  mutate(species_name = paste0(Rank6, " ", Rank7))


# Plotting log-fold change?

lfc <- ggplot(top20_for_plotting, aes(x = lfc_dietvegan, y = species_name)) + geom_point(size = 2, colour="blue") + geom_vline(xintercept = 0, color = "red") + theme(legend.position = "none", axis.text.y = element_text(size=6)) + labs(y="Species name", x="Log-fold change between omnivores and vegans")

# most species decreased in vegans compared to omnivore, but not significantly

ggsave("../figures/logfoldplot.png", lfc)

## End
