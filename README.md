# Microbiome-differential-abundance

## Introduction

The human gut contains rich microbial communities (McBurney et al., 2019). Changes in the gut microbiome can affect various aspects of health and have been association with a range of diseases including irritable bowel syndrome, Type-2 diabetes, and Alzheimer's disease (McBurney et al., 2019; Ghaisas, Maher, & Kanthasamy, 2016). In particular, diet is an incredibly influential factor on the gut microbiome (David et al., 2014). De Filippis et al. (2019) investigated and compared the gut microbiome of individuals who follow a vegan diet, a vegetarian diet, or an omnivore diet, and found that different diets resulted in different strains of the bacterium _Prevotella copri_ being present in the gut, showing that even subtle differences in diet have large effects.

Given that the composition of the gut microbiome can greatly affect health, studying it is of utmost importance. In order to explore the gut microbiome and explore its associations with disease, the genomes of the microbes present must be sequenced, which is known as metagenomics (Wang et al., 2015). Sequencing all the microbes present in a sample enables identification of the species present and the abundance of each species (Wang et al., 2015). 

The first step in analysis is classification of species in the sample. Some tools that exist for this purpose include Kraken2 (Wood, Lu, & Langmead, 2019), MetaPhlAn4 (Blanco-Míguez et al., 2023), and Kaiju (Menzel, Ng, & Krogh, 2016). Kraken2 improves upon the earlier Kraken1, which had extremely high memory requirements that represented a barrier to usage. Kraken2 uses k-mer matches, which ensure that it can classify species very quickly and accurately (Wood, Lu, & Langmead, 2019). MetaPhlAn 4 is a taxonomic classifier that achieves more accurate classification via adding more data to its databases (Blanco-Míguez et al., 2023). In a review of taxonomic classifiers by Edwin et al. (2024), MetaPhlAn 3 and 4 were not very sensitive and left 94.5% and 90.5% of species unclassified, respectively. Kaiju is another taxonomic classifier which instead of using k-mer based matching, uses the Burrows-Wheeler transform to perform exact string matching, which greatly improves sensitivity and precision of matches (Menzel, Ng, & Krogh, 2016). Edwin et al. (2024) found that Kaiju performed well, but noted that performance of all classifiers (i.e. ability to classify sequences correctly) was affected by the database used for the classification. In their analysis, the best-performing classifier was Kraken2 used in conjunction with Bracken, which estimates species abundance (Lu et al., 2017). For this reason, I chose to use Kraken2 and Bracken in this analysis. 
	
The next steps of the workflow often involve computing measures of alpha and beta diversity. Alpha diversity describes the differences in species diversity within a sample, while beta diversity measures the difference in diversity between samples (Kers & Saccenti, 2022). Cassol, Ibañez, and Bustamante (2025) suggest that alpha diversity can be divided into four categories to provide a comprehensive measure of diversity: richness, dominance, phylogenetics, and information. To represent three of these categories, I will compute the Chao1 (Chao, 1984), Berger-Parker (Caruso et al., 2007), and Shannon metrics (Lemos et al., 2011), as recommended. Chao1 is particularly recommended for small sample sizes by Kers and Saccenti (2022). It is important to compute both alpha and beta diversity (Kers & Saccenti, 2022), so I will also compute Bray-Curtis dissimilarity and the Jaccard measure, as these are the most commonly used measures of beta diversity, and Bray-Curtis distance in particular is suitable for small sample sizes (Kers & Saccenti, 2022). 

Lastly, differential abundance analysis will show whether there are statistically significant differences in the abundance of taxa between samples. Some packages that can be used for this include Aldex3, ANCOM-BC, and established differential expression packages like edgeR (Robinson, McCarthy, & Smyth, 2010) and DEseq2 (Love, Huber, & Anders, 2014). Both edgeR and DEseq2 are useful and accurate tools, but since they are designed for analysis of RNA-seq data of single organisms, are prone to false positives (Nearing et al., 2022; Wirbel et al., 2024). Aldex2 (and the recently released Aldex3; Fernandes et al., 2013) and ANCOM-BC (Lin & Peddada, 2020) are much more conservative and designed for differential abundance analyses specificically. Nearing et al. 2022 found that Aldex2 was able to correctly identify differentially abundant species across five datasets more consistently than other tools. However, given the similarity in satisfactory performance between Aldex2 and ANCOM-BC (Nearing et al., 2022, Wirbel et al., 2024) and my prior familiarity with ANCOM-BC, as well as the interpretability of its output, I will use this for my analyses.
	
I will use metagenomic sequencing data from De Filippis et al. (2019) to compare gut microbiome diversity and abundance between eight humans who follow either a vegan or omnivore diet. 


## Methods

I downloaded the raw metagenomic sequencing reads from the gut microbiome of four vegans and four omnivores from NCBI SRA (accessions: SRR8146935, SRR8146936, SRR8146938, SRR8146956, SRR8146944, SRR8146951, SRR8146952, SRR8146954; De Filippis et al., 2019). I used FastQC v0.12.1 to check read quality and cutadapt v5.2 to trim reads with quality score below 20 (Andrews, 2010; Martin, 2011). 

I used the core_nt database from NCBI (download link: https://benlangmead.github.io/aws-indexes/k2) to classify the reads with Kraken2 v2.1.6, specifying a confidence score of 0.15 (Wood, Lu, & Langmead, 2019). I chose the most comprehensive available database and a moderate confidence score as recommended by Liu et al. (2024). I did not specify the memory-mapping flag in order to speed up analyses. I estimated species abundance with Bracken v3.0 (Lu et al., 2017). The resulting species abundance report files were converted into BIOM format for reading into R with kraken-biom v1.0.1 (Dabdoub, 2016).

Once the data was loaded into R with the biomformat package v1.38.3 (McMurdie & Paulson, 2026), I calculated and visualised relative abundance and species rarefaction curves with the R packages phyloseq v1.54.2 and vegan v2.7-3 (McMurdie & Holmes, 2013; Oksanen et al., 2026). I calculated the Chao1 (Chao, 1984), Shannon (Lemos et al., 2011), and Berger-Parker indices (Caruso et al., 2007), and tested for significant differences between diets with a t-test. I also calculated the Bray-Curtis (Bray & Curtis, 1957) and Jaccard distances (Jaccard, 1912), as well as a PERMANOVA for both measures.

I calculated differential abundance of microbial species between samples with ANCOM-BC v2.12.1 (Lin & Peddada, 2020). 


## Results

The reads were of high quality according to FastQC, with a dropoff in read quality in the second read of each sample, as is expected (Tan et al., 2019). I trimmed sequences with Phred score < 20, and the resulting sequences were all of high quality. 
Relative abundance of the top twenty most abundant taxa between samples can be seen in Fig. 1. Two of the four omnivores have _Faecalibacterium prausnitzii_ as the most abundant species and three of four vegans (and one omnivore) have _S. copri_ as the most abundant species.
	
<img width="1432" height="1467" alt="rel_abundance_plot" src="https://github.com/user-attachments/assets/2c9e17ef-9e0a-4976-a541-38630591b572" />

_Figure 1. Relative abundance of microbial species in each sample, with the diet of each individual indicated on the x-axis. Only the top 14 most abundant species across all samples are shown, for visual clarity._

Per sample alpha diversity is shown in Fig. 2. The Chao1 index, Shannon index, and Berger-Parker index were calculated for each sample. There were no significant differences between diets for any diversity measure (Chao1: t(5.9348)= 0.078267, p > 0.05; Shannon: t(4.2584)=0.81739, p > 0.05; Berger-Parker: t(4.1397)=-1.4049, p > 0.05). 

<img width="1432" height="1467" alt="alpha_diversity_plot" src="https://github.com/user-attachments/assets/72d76c81-9bd5-43f4-a304-df6db191befd" />

_Figure 2. The Chao1, Shannon, and Berger-Parker diversity indices calculated for each sample, where each point is coloured by diet._

Measures of beta diversity are shown in Fig. 3 and 4. Principal coordinate analysis plots of Bray-Curtis and Jaccard distance show that there is some clustering of samples by diet, but the PERMANOVA for each measure was not significant (p < 0.05). 

<img width="1432" height="1467" alt="braycurtis_plot" src="https://github.com/user-attachments/assets/c6974367-78ae-4f5d-ad8b-65a3705d7da0" />

_Figure 3. Principal coordinates plot of Bray-Curtis distance calculated across the eight samples, where each point represents a sample and points are coloured by diet. The top two principal components are shown._

<img width="1432" height="1467" alt="jaccard_plot" src="https://github.com/user-attachments/assets/0843c9e9-5745-4b85-b93f-9c3dee485ed3" />

_Figure 3. Principal coordinates plot of Jaccard distance calculated across the eight samples, where each point represents a sample and points are coloured by diet. The top two principal components are shown._

There were no taxa that were significantly differentially abundant between samples. Log-fold change of the top twenty taxa with the largest absolute log-fold change is shown in Fig. 5.
	
<img width="1432" height="1467" alt="logfoldplot" src="https://github.com/user-attachments/assets/0fd743f6-a2f8-436c-b2b8-e3845123e3cc" />

_Figure 5. Log-fold change of the top twenty taxa with the largest absolute log-fold change, shown in individuals that followed a vegan diet relative to individuals that followed an omnivore diet._

## Discussion

The bacterial species with the highest relative abundance in most vegans and one of the omnivores was _Prevotella copri_ (Fig. 1). This finding is supported by the literature, as _P. copri_ is one of the most abundant species in the human gut microbiome, but abundance can also vary between individuals (Falony et al., 2016). _P. copri_ in the gut has been associated with diets higher in fibre (De Filippo et al., 2010). The Prevotella enterotype (at the genus level) has been associated with consumption of a diet high in carbohydrates, and was the most prevalent enterotype in vegetarians (Wu et al., 2011). The association of relative abundance of _P. copri_ with plant-rich diets is likely because _P. copri_ breaks down polysaccharides from plant sources (Fehlner-Peach et al., 2019). Interestingly, relative abundance of _P. copri_ has been associated with both positive and negative disease outcomes (Abdelsalam, Hegazy, & Aziz, 2023). For example, there was a lower relative abundance of _P. copri_ found in the gut microbiome of patients with Parkinson’s disease compared to healthy controls (Bedarf et al., 2017). On the other hand, some studies have found that patients with rheumatoid arthritis had a higher relative abundance of _P. copri_ (Scher et al., 2013). Overall, investigating the relative abundance of species offers a fascinating look into the gut microbiome. 

Measures of alpha diversity were not significantly different in omnivores compared to vegans. Kers and Saccenti (2022) conclude that beta diversity measures are more sensitive to differences between groups than alpha diversity metrics. A principal coordinates analysis of Bray-Curtis distance (Fig. 3) and Jaccard distance (Fig. 4) shows slight clustering of samples by diet. However, a PERMANOVA for both distance measures was not significant. This may indicate that a difference in gut microbiome composition between diets does exist, as De Filippis et al. (2019) found, but the small sample size limited detection of a significant difference.

The small sample size may also account for the lack of significant differences in abundance of species between omnivores and vegans. No taxa were significantly differentially abundant between diets. However, I will discuss some taxa that had the largest log-fold change, although caution should be taken in interpreting these results as they were not statistically significant. I found that _Bacteroides eggerthii_ exhibited a negative log-fold change in vegans compared to omnivores, indicating that the abundance of this species was lower in vegans compared to omnivores. David et al. (2014) found that abundance of the Bacteroides genus was higher in study subjects that adhered to an animal-based diet, and Wu et al. (2011) found that a protein and animal-fat rich diet was associated with the Bacteroides enterotype. Given the literature support for this association, further investigation with a larger sample size could be informative. Overall, the results reflect the strong impact of diet on the gut microbiome and highlight the need for increased sample size to potentially uncover significant associations of taxa to diet. 

## References

Abdelsalam, N. A., Hegazy, S. M., & Aziz, R. K. (2023). The curious case of Prevotella copri. Gut Microbes, 15(2), 2249152. https://doi.org/10.1080/19490976.2023.2249152

Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data [Online]. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/.

Bedarf, J. R., Hildebrand, F., Coelho, L. P., Sunagawa, S., Bahram, M., Goeser, F., Bork, P., & Wüllner, U. (2017). Functional implications of microbial and viral gut metagenome changes in early stage L-DOPA-naïve Parkinson’s disease patients. Genome Medicine, 9, 39. https://doi.org/10.1186/s13073-017-0428-y

Blanco-Míguez, A., Beghini, F., Cumbo, F., McIver, L. J., Thompson, K. N., Zolfo, M., Manghi, P., Dubois, L., Huang, K. D., Thomas, A. M., Nickols, W. A., Piccinno, G., Piperni, E., Punčochář, M., Valles-Colomer, M., Tett, A., Giordano, F., Davies, R., Wolf, J., … Segata, N. (2023). Extending and improving metagenomic taxonomic profiling with uncharacterized species using MetaPhlAn 4. Nature Biotechnology, 41(11), 1633–1644. https://doi.org/10.1038/s41587-023-01688-w

Bray, J. R., & Curtis, J. T. (1957). An ordination of the upland forest communities of southern wisconsin. Ecological Monographs, 27(4), 325–349. https://doi.org/10.2307/1942268

Caruso, T., Pigino, G., Bernini, F., Bargagli, R., & Migliorini, M. (2007). The Berger–Parker index as an effective tool for monitoring the biodiversity of disturbed soils: A case study on Mediterranean oribatid (Acari: oribatida) assemblages. Biodiversity and Conservation, 16(12), 3277–3285. https://doi.org/10.1007/s10531-006-9137-3

Cassol, I., Ibañez, M., & Bustamante, J. P. (2025). Key features and guidelines for the application of microbial alpha diversity metrics. Scientific Reports, 15(1), 622. https://doi.org/10.1038/s41598-024-77864-y

Chao, A. (1984). Nonparametric estimation of the number of classes in a population. Scand. J. Stat. 11 265–270.

David, L. A., Maurice, C. F., Carmody, R. N., Gootenberg, D. B., Button, J. E., Wolfe, B. E., Ling, A. V., Devlin, A. S., Varma, Y., Fischbach, M. A., Biddinger, S. B., Dutton, R. J., & Turnbaugh, P. J. (2014). Diet rapidly and reproducibly alters the human gut microbiome. Nature, 505(7484), 559–563. https://doi.org/10.1038/nature12820

Dabdoub, S.M. (2016). kraken-biom: Enabling interoperative format conversion for Kraken results (Version 1.2) [Software]. Available at https://github.com/smdabdoub/kraken-biom.

De Filippis, F., Pasolli, E., Tett, A., Tarallo, S., Naccarati, A., De Angelis, M., Neviani, E., Cocolin, L., Gobbetti, M., Segata, N., & Ercolini, D. (2019). Distinct genetic and functional traits of human intestinal prevotella copri strains are associated with different habitual diets. Cell Host & Microbe, 25(3), 444-453.e3. https://doi.org/10.1016/j.chom.2019.01.004

De Filippo, C., Cavalieri, D., Di Paola, M., Ramazzotti, M., Poullet, J. B., Massart, S., Collini, S., Pieraccini, G., & Lionetti, P. (2010). Impact of diet in shaping gut microbiota revealed by a comparative study in children from Europe and rural Africa. Proceedings of the National Academy of Sciences of the United States of America, 107(33), 14691–14696. https://doi.org/10.1073/pnas.1005963107

Dos Santos, S. J., & Gloor, G. B. (2026). Incorporating scale uncertainty into differential expression analyses using aldex2. Current Protocols, 6(2), e70307. https://doi.org/10.1002/cpz1.70307

Fehlner-Peach, H., Magnabosco, C., Raghavan, V., Scher, J. U., Tett, A., Cox, L. M., Gottsegen, C., Watters, A., Wiltshire-Gordon, J. D., Segata, N., Bonneau, R., & Littman, D. R. (2019). Distinct polysaccharide utilization profiles of human intestinal Prevotella copri isolates. Cell Host & Microbe, 26(5), 680-690.e5. https://doi.org/10.1016/j.chom.2019.10.013

Fernandes, A. D., Macklaim, J. M., Linn, T. G., Reid, G., & Gloor, G. B. (2013). Anova-like differential expression (Aldex) analysis for mixed population rna-seq. PLOS ONE, 8(7), e67019. https://doi.org/10.1371/journal.pone.0067019

Ghaisas, S., Maher, J., & Kanthasamy, A. (2016). Gut microbiome in health and disease: Linking the microbiome–gut–brain axis and environmental factors in the pathogenesis of systemic and neurodegenerative diseases. Pharmacology & Therapeutics, 158, 52–62. https://doi.org/10.1016/j.pharmthera.2015.11.012

Jaccard, P. (1912). The distribution of the flora in the alpine zone. New Phytol, 11, 37–50. 10.1111/j.1469-8137.1912.tb05611.x

Kers, J. G., & Saccenti, E. (2022). The power of microbiome studies: Some considerations on which alpha and beta metrics to use and how to report results. Frontiers in Microbiology, 12, 796025. https://doi.org/10.3389/fmicb.2021.796025

Lemos, L. N., Fulthorpe, R. R., Triplett, E. W., & Roesch, L. F. W. (2011). Rethinking microbial diversity analysis in the high-throughput sequencing era. Journal of Microbiological Methods, 86(1), 42–51. https://doi.org/10.1016/j.mimet.2011.03.014

Lin, H., & Peddada, S. D. (2020). Analysis of compositions of microbiomes with bias correction. Nature Communications, 11(1), 3514. https://doi.org/10.1038/s41467-020-17041-7

Liu, Y., Ghaffari, M. H., Ma, T., & Tu, Y. (2024). Impact of database choice and confidence score on the performance of taxonomic classification using Kraken2. aBIOTECH, 5(4), 465–475. https://doi.org/10.1007/s42994-024-00178-0

Lu, J., Breitwieser, F. P., Thielen, P., & Salzberg, S. L. (2017). Bracken: Estimating species abundance in metagenomics data. PeerJ. Computer Science, 3, e104. https://doi.org/10.7717/peerj-cs.104

Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.Journal, 17(1), 10. https://doi.org/10.14806/ej.17.1.200

McBurney, M. I., Davis, C., Fraser, C. M., Schneeman, B. O., Huttenhower, C., Verbeke, K., Walter, J., & Latulippe, M. E. (2019). Establishing what constitutes a healthy human gut microbiome: State of the science, regulatory considerations, and future directions. The Journal of Nutrition, 149(11), 1882–1895. https://doi.org/10.1093/jn/nxz154

McMurdie, P., & Paulson, J. (2026). biomformat: An interface package for the BIOM file format (R package version 1.38.3). https://doi.org/10.18129/B9.bioc.biomformat

McMurdie, P., & Holmes, S. (2013). phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. PLOS ONE, 8(4), e61217. https://doi.org/10.1371/journal.pone.0061217

Menzel, P., Ng, K. L., & Krogh, A. (2016). Fast and sensitive taxonomic classification for metagenomics with Kaiju. Nature Communications, 7(1), 11257. https://doi.org/10.1038/ncomms11257

Nearing, J. T., Douglas, G. M., Hayes, M. G., MacDonald, J., Desai, D. K., Allward, N., Jones, C. M. A., Wright, R. J., Dhanani, A. S., Comeau, A. M., & Langille, M. G. I. (2022). Microbiome differential abundance methods produce different results across 38 datasets. Nature Communications, 13, 342. https://doi.org/10.1038/s41467-022-28034-z

Oksanen, J., Simpson, G., Blanchet, F. G., Kindt, R., Legendre, P., Minchin, P. R., O’Hara, R. B., Solymos, P., Stevens, M. H. H., Szoecs, E., Wagner, H., Barbour, M., Bedward, M., Bolker, B., Borcard, D., Borman, T., Carvalho, G., Chirico, M., De Cáceres, M., Durand, S., Evangelista, H., FitzJohn, R., Friendly, M., Furneaux, B., Hannigan, G., Hill, M. O., Lahti, L., Martino, C., McGlinn, D., Ouellette, M., Ribeiro Cunha, E., Smith, T., Stier, A., ter Braak, C. J. F., & Weedon, J. (2026). vegan: Community ecology package (R package version 2.7-3). https://doi.org/10.32614/CRAN.package.vegan

Robinson, M. D., McCarthy, D. J., & Smyth, G. K. (2010). Edger: A bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics, 26(1), 139–140. https://doi.org/10.1093/bioinformatics/btp616

Scher, J. U., Sczesnak, A., Longman, R. S., Segata, N., Ubeda, C., Bielski, C., Rostron, T., Cerundolo, V., Pamer, E. G., Abramson, S. B., Huttenhower, C., & Littman, D. R. (2013). Expansion of intestinal Prevotella copri correlates with enhanced susceptibility to arthritis. eLife, 2, e01202. https://doi.org/10.7554/eLife.01202

Tan, G., Opitz, L., Schlapbach, R., & Rehrauer, H. (2019). Long fragments achieve lower base quality in Illumina paired-end sequencing. Scientific Reports, 9, 2856. https://doi.org/10.1038/s41598-019-39076-7

Truong, D. T., Franzosa, E. A., Tickle, T. L., Scholz, M., Weingart, G., Pasolli, E., Tett, A., Huttenhower, C., & Segata, N. (2015). MetaPhlAn2 for enhanced metagenomic taxonomic profiling. Nature Methods, 12(10), 902–903. https://doi.org/10.1038/nmeth.3589

Wang, W.-L., Xu, S.-Y., Ren, Z.-G., Tao, L., Jiang, J.-W., & Zheng, S.-S. (2015). Application of metagenomics in the human gut microbiome. World Journal of Gastroenterology : WJG, 21(3), 803–814. https://doi.org/10.3748/wjg.v21.i3.803

Wirbel, J., Essex, M., Forslund, S. K., & Zeller, G. (2024). A realistic benchmark for differential abundance testing and confounder adjustment in human microbiome studies. Genome Biology, 25, 247. https://doi.org/10.1186/s13059-024-03390-9

Wood, D. E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. Genome Biology, 20(1), 257. https://doi.org/10.1186/s13059-019-1891-0

Wu, G. D., Chen, J., Hoffmann, C., Bittinger, K., Chen, Y.-Y., Keilbaugh, S. A., Bewtra, M., Knights, D., Walters, W. A., Knight, R., Sinha, R., Gilroy, E., Gupta, K., Baldassano, R., Nessel, L., Li, H., Bushman, F. D., & Lewis, J. D. (2011). Linking long-term dietary patterns with gut microbial enterotypes. Science (New York, N.y.), 334(6052), 105–108. https://doi.org/10.1126/science.1208344


