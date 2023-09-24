# Blaustein_RAG_2023

Scripts and data used to characterize microbial diversity in metagenomes of patients with RAG-deficiency.

**Reference:** Blaustein et al. 2023. Expanded microbiome niches of RAG-deficient patients. Cell Reports Medicine. (_accepted_)

# Contents

**code**
* **metagenome_workflow**: List of bioinformatics scripts used in the analysis, i.e., read-based metagenome analysis, metagenome assembly, metagenome-assembled genome (MAG) recovery, viral genome recovery, read-mapping, and supporting files
* **RAG_metagenomes.R**: R scripts to reproduce downstream analysis and generate data visualizations reported in the study

**data** (i.e., working directory for **RAG_metagenomes.R**)

* **ARGs**

    * Antimicrobial resistance gene (ARG) family relative abundances or Reads Per Kilobase per Million mapped reads (RPKMs) in all samples (**shortbred_ARGs_rag_hv.csv**; **shortbred_ARG_hv_ABX_control.txt**)
  
    * Sets of ARGs found to be enriched in metagenomes of adult and pediatric patients with RAG-deficiency as compared with those of healthy volunteers and children (**ARGs_stool_adult_manual_edit.csv**; **ARGs_skin_adult_manual_edit.csv**; **ARGs_skin_ped_manual_edit.csv**)

* **HPVs**
  
    * Human papillomavirus (HPV) RPKMs in metagenomes (**hpviewer_counts.csv**)
  
    * Newick tree of distinct HPV contigs recovered from RAG-deficient patients with context to the Papillomavirus Episteme (PaVE) (**hpv_aln.newick**; **hpv_tree_labels.csv**)

* **MAGs_read_mapping**
  
    * Metagenome-assembled genomes (MAGs) recovered from patients with RAG-deficiency (**rag_mags_taxa.csv**) along with newick tree based on protein sequence alignment (**gtbtk_bac_aln.newick**)

    * Uniquely mapped reads of microbial genomes within metagenomes, i.e., used to predict presence and abundance of microbial taxa (**read_map_unique.txt.zip**; **read_map_taxa_list.csv**)

    * Stats for fractions properly-paired uniquely mapped reads and total Kraken2-classified reads (**read_map_stats.csv**)

* **metadata**

    * Metadata associated with samples (**rag_hv_metadata.csv**; **hv_ABX_control_metadata.txt**)

# Modules Required

_MAG recovery and characterization_
  * SPAdes v3.15.0
  * MetaWRAP v1.2.2
  * GUNC v1.0.5
  * dRep v2.6.2
  * GTDB-Tk v1.3.0
  * FastTree v2.1.10
  * EukCC v2.1.1
  * MASH v2.3
  * MUMmer v4.0.0beta2
 
_Viral genome recovery and characterization_
  * VirFinder v1.1
  * Virsorter2 v2.1
  * CheckV v0.7.0
  * CD-HIT v4.6.8
  * DemoVir
  * prodigal v2.6.3
  * hmmer v3.3.2
  * MUSCLE v5.0.1428

_Read mapping of microbial genomes_
  * BWA-MEM v0.7.17
  * SAMtools v1.15

_Reference-based microbial taxonomic and functional classification_
  * Kraken2 v2.1.2
  * HPViewer
  * ShortBRED v0.9.5

_Metatranscriptomic analysis_
  * ShortBRED v0.9.5
  * Bowtie2 v2.4.5
  * BLAST v2.8.0
  * RagTag v2.1.0
  * Norovirus typing tool v2.0
  * ggtree v3.2.1
  * Clustal Omega
  * STAR v2.7.8a
  * HOMER v4.11.1
  
