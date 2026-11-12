# Wye_pond_metagenomes

Scripts and data used to characterize microbial diversity in metagenomes of agricultureal research pond at Wye Research and Education Center.

Reference: Blaustein RA, Smith JE, Toro M, Pachepsky Y, Stocker MD. 2025. Water metagenomes reflect physicochemical water quality throughout a model agricultural pond. _Frontiers in Microbiology_ 16:1535096. https://doi.org/10.3389/fmicb.2025.1535096

# Contents

**code**
* **command_line**: list of bioinformatics scripts used in the analysis
* **R**: R script to reproduce statistical analysis and generate data visualizations reported in the study

**R_input_data**
* **annotation**: functional content in metagenomes (i.e., metacyc pathways, antimicrobial resistance genes, virulence factors)
* **metadata**: pond water and sediment sample metadata
* **taxa**: kraken2/bracken table of microbial taxa counts
* **tree**: metagenome-assembled genome (MAG) alignment

# Modules Required
_Preprocessing_
* Kneaddata v.0.6.1

_Microbial taxonomic and functional classification_
* Kraken2 v2.1.2
* Braken
* HUMAnN2 v.0.11.1
* ShortBRED v0.9.5

_MAG recovery and characterization_
* SPAdes v.3.15.0
* MetaWRAP v.1.2.2
* GUNC v1.0.5
* dRep v2.6.2
* GTDB-Tk v1.3.0
* FastTree v.2.1.10






_AMRFinderPlus v3.12.8_
