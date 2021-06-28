# combat_sda_processing
Python code for the preprocessing of the COvid-19 Multi-omics Blood ATlas (COMBAT) data for SDA tensor and matrix decomposition and the processing of the SDA output data as decribed below. More information can be found in our preprint: https://www.medrxiv.org/content/10.1101/2021.05.11.21256877v1


## Description of data preprocessing: 
A tensor and matrix decomposition method, SDA, as defined in (Hore et al., 2016), was used to integrate 152 samples from gene expression of the bulk-RNA-seq, gene expression from the pseudo-bulk-RNA-seq and cell composition data from CITEseq, cell composition data from the CyTOF all cells and the CyTOF depleted cells, and protein abundances from Luminex and timsTOF.  

The Bulk Total RNAseq (9 missing samples) and the pseudobulk from 10x Genomics Chromium scRNAseq (22 missing samples), was combined into a three dimensional tensor consisting of 152 samples by 9 tissue types by 14,989 genes (which passed QC in both datasets). This expression was normalised by sample in each tissue type by log2-transformation of counts-per-million + 1.  

The number of cells per cell types as defined by 10x Genomics Chromium CITEseq (one two-dimensional matrix of 152 samples by 64 cell types, with 22 samples missing) and CyTOF (a two-dimensional matrix of 152 samples by 10 cell types for all cells, with 21 samples missing, and a two-dimensional matrices of 152 samples by 51 cell types for all cells, with 20 samples missing. We filtered out any samples with fewer than 500 cells in any matrix. The data in each matrix was normalised by a log2 transformation of counts-per-million + 1.  

The proteomics data from luminex (in a two-dimensional matrix of 152 samples by 51 proteins, with 20 samples missing) and timsTOF (in a two-dimensional matrix of 152 samples by 105 proteins, with 17 samples missing) was used, with the data normalised as described in the luminex and timsTOF sections.  

The output data is in the format ready for input to SDA - no column or row labels (they are provided separately here as ```md_*.txt```) and space delimited. ```rna.txt``` is correctly formatted so that SDA can read it as a three-dimensional tensor. The other data files are two-dimensional matrices.

## Description of data analysis:
As described in (Hore et al., 2016), to find robust components we ran the tensor and matrix decomposition ten times for 1000 components. Once again, similar to the Hore et al., 2016 method the absolute correlation (r) was calculated for the sample scores for each pair of components, clustered using Hierarchical clustering on 1-r (dissimilarity measure) and formed flat clusters in which the components in each flat cluster have no greater a cophenetic distance than 0.4. We chose the flat clusters that had components from at least 5 of the 10 runs. The final sample, tissue and gene or protein or cell score was the mean of all the components within the chosen clusters. This resulted in 381 clusters. 

The resulting matrices with the loading scores are provided as ```final_*.csv``` files with the corresponding posterior inclusion probabilities as ```pip_final_*.csv``` files. 
