# Files description

The two R scripts were used to predict mitochondria localization of candidate genes in Zhang et al. manuscript.

GSEA_Matrix.R creates a matrix by collapsing the Normalized Enrichment Scores of GSEA analysis. It requires 2 inputs :
    A list of gene IDs of interest (for example Mitochondrial genes)
    An RNA-seq normalized count matrix in which the analysis will be done.

Mito_prediction_kmeans.R gives score to a candidate list either 1 or 0 for predicted mitochondria location. This program requires the NES matrices of Mitocarta, negative control (for example random transcription factors) and candidates.
