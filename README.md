# Microbiome-differential-abundance

## Introduction

The human gut contains rich microbial communities. Changes in the gut microbiome can affect various aspects of health (McBurney et al., 2019), and have been association with a range of diseases including irritable bowel syndrome, type-2 diabetes, and Alzheimer's disease (Ghaisas, Maher, & Kanthasamy, 2016).
In particular, diet is an incredibly influential factor on the gut microbiome (David et al., 2014). De Filippis et al. (2019) investigated and compared the gut microbiome of individuals who follow a vegan diet, a vegetarian diet, or an omnivore diet. Given that the composition of the gut microbiome can greatly affect health, studying it is of utmost importance. In order to explore the gut microbiome and explore its associations with disease, the genomes of the microbes present must be sequenced, which is known as metagenomics. Sequencing all the microbes present in a sample enables identification of the species present and the abundance of each species (Wang et al., 2015). 

The first step in analysis is classification of species in the sample. Some tools that exist for this purpose include Kraken2 (Wood, Lu, & Langmead, 2019), MetaPhlAn4 (Blanco-Míguez et al., 2023), and Kaiju (Menzel, Ng, & Krogh, 2016). Kraken2 improves upon the earlier Kraken1, which has extremely high memory requirements that represented a barrier to usage. Kraken2 uses k-mer matches, which ensure that it can classify species very quickly and accurately (Wood, Lu, & Langmead, 2019). MetaPhlAn 4 is a taxonomic classifier that achieves more accurate classification via adding more data to its databases (Blanco-Míguez et al., 2023). In a review of taxonomic classifiers by Edwin et al. (2024), MetaPhlAn 3 and 4 were not very sensitive and left 94.5% and 90.5% of species unclassified, respectively. Kaiju is another taxonomic classifier which instead of using k-mer based matching, uses the Burrows-Wheeler transform to perform exact string matching, which greatly improves sensitivity and precision of matches (Menzel, Ng, & Krogh, 2016). Edwin et al. (2024) found that Kaiju performed well, but noted that performance of all classifiers (i.e. ability to classify sequences correctly) was affected by the database used for the classification. In their analysis, the best-performing classifier was Kraken2 used in conjunction with Bracken, which estimates species abundance (Lu et al., 2017). 

The next steps of the workflow often involve computing measures of alpha and beta diversity. Alpha diversity describes the differences in species diversity within a sample, while beta diversity measures the difference in diversity between samples. Cassol, Ibañez, and Bustamante (2025) suggest that alpha diversity can be divided into four categories to provide a comprehensive measure of diversity: richness, dominance, phylogenetics, and information. To represent these categories, I will compute the Chao1 (Chao, 1984, Berger-Parker, Faith, and Shannon metrics, as recommended (Chao1 is particularly recommened for small sample sizes by Kers and Saccenti (2022)).

It is important to compute both alpha and beta diversity (Kers & Saccenti, 2022), so I will also compute Bray-Curtis dissimilarity and the Jaccard measure, as these most commonly used, and Bray-Curtis in particular is suitable for small sample sizes (Kers & Saccenti, 2022). 

Lastly, differential abundance analysis will show whether there are statistically significant differences between samples. Some packages that can be used for this include Aldex3, AncomBC, and established differential expression packages like edgeR and DEseq2. Both edgeR and DEseq2 are useful and accurate tools, but since they are designed for analysis of RNA-seq data of single organisms, are prone to false positives (Nearing et al., 2022, Wirbel et al., 2024). Aldex2 (Fernandes et al., 2013) and ANCOM-BC (Lin & Peddada, 2020) are much more conservative and designed for differential abundance analyses specificically. Nearing et al. 2022 found that Aldex2 was able to correctly identify differentially abundant species across five datasets more consistently than other tools. However, given the similarity in satisfactory performance between Aldex2 and ANCOM-BC (Nearing et al., 2022, Wirbel et al., 2024) and my prior familiarity with ANCOM-BC, as well as the interpretability of its output, I will use this for my analyses.
	
I will use metagenomic sequencing data from De Filippis et al. (2019) to compare gut microbiome diversity and abundance between humans who follow either a vegan or omnivore diet. 

## Methods

I downloaded the raw metagenomic sequencing reads from the gut microbiome of four vegans and four omnivores from NCBI SRA (accessions: SRR8146935, SRR8146936, SRR8146938, SRR8146956, SRR8146944, SRR8146951, SRR8146952, SRR8146954) (De Filippis et al., 2019). I used FastQC v0.12.1 to check read quality and cutadapt v5.2 to trim reads with quality score below 20. 

I used the core_nt database from NCBI (download link: https://benlangmead.github.io/aws-indexes/k2) to classify the reads with Kraken2 v2.1.6, specifying a confidence score of 0.15. I chose the most comprehensive available database and a moderate confidence score as recommended by Liu et al. (2024). I did not specify the memory-mapping flag in order to speed up analyses. I estimated species abundance with Bracken v3.0 (Lu et al., 2017). The resulting species abundance report files were converted into BIOM format for reading into R with kraken-biom v1.0.1.

Once the data was loaded into R with the biomformat package v1.38.3, I calculated and visualised relative abundance and species rarefaction curves with the R packages phyloseq v1.54.2 and vegan v2.7-3. I calculated the following measures of alpha diversity:  Chao1, Shannon, and Berger-Parker indices, and the following measures of beta diversity: Bray-Curtis and Jaccard distance, as well as a PERMANOVA with Bray-Curtis distance.

I calculated differential abundance of microbial species between samples with ANCOM-BC v2.12.1 (Lin & Peddada, 2020). 

## Results

<img width="1432" height="1467" alt="rel_abundance_plot" src="https://github.com/user-attachments/assets/2c9e17ef-9e0a-4976-a541-38630591b572" />





<img width="1432" height="1467" alt="alpha_diversity_plot" src="https://github.com/user-attachments/assets/72d76c81-9bd5-43f4-a304-df6db191befd" />

<img width="1432" height="1467" alt="braycurtis_plot" src="https://github.com/user-attachments/assets/c6974367-78ae-4f5d-ad8b-65a3705d7da0" />

<img width="1432" height="1467" alt="jaccard_plot" src="https://github.com/user-attachments/assets/0843c9e9-5745-4b85-b93f-9c3dee485ed3" />

## Discussion
