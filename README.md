# SpotEnrich
SpotEnrich is a Python toolkit for **spatial transcriptomics analysis**, specifically designed to interpret **cell-type distribution patterns** produced by tools such as Cell2location. 
It classifies the spots into three categories: 
1. Enrich spots (defined as the specific celltype spots dominated by a single cell type)
2. Mix spot (containing multiple cell types)
3. Unknown (unclassifiable spots)

using a series of adaptive thresholds for categorization.
