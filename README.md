# README

This repository contains the datasets and scripts used to create the figures in:

**Hajibabaei et al., 2019. COI metabarcoding primer choice affects richness and recovery of indicator taxa in freshwater systems.  PLoS ONE, accepted.**

## Data analysis outline

1. Raw reads are available from the NCBI SRA # PRJNA545426 and were processed using the SCVUC v2.1 pipeline available from https://github.com/EcoBiomics-Zoobiome/SCVUC_COI_metabarcode_pipeline . 

2. Datasets: The taxonomic assignments are available in taxonomy.matrix.gz .  Site locations are in Sites.csv . FASTA files for the denoised ESVs are in denoised_ESVs_supmat.tar.gz .

3. Figures were pepared in R as follows:
  * Fig 2 was generated with Fig2_multimarker_richness_accumulation.R from taxonomy.matrix.gz . 
  * Fig 3 was generated with Fig3_site_waterquality_indicators.R from taxonomy.matrix.gz . 
  * Fig 4 was generated with Fig4_site_indicator_species_heatmap.R from taxonomy.matrix.gz.  Final editing was done in Inkscape.
  * Fig 5 was generated with Fig5_NMDS.R from taxonomy.matrix.gz . 
  
4. Supplementary figures were prepared in R as follows:
  * Fig S1 was generated with FigS1_map.R from Sites.csv . 
  * Fig S2 was generated with FigS2_TaxSummary.R from taxonomy.matrix.gz . 
  * Fig S3 was generated with FigS3_rarefaction.R from taxonomy.matrix.gz . 
  * Fig S4 was generated with FigS4_MedianRichness.R from taxonomy.matrix.gz . 
  * Fig S5 was generated with FigS5_TotalRichness.R from taxonomy.matrix.gz . 

## Acknowledgements

I would like to acknowledge funding from the Canadian government from the Genomics Research and Development Initiative (GRDI) Ecobiomics project.

Last updated: Sept. 9, 2019.
