# SpotEnrich
SpotEnrich is a Python toolkit for **spatial transcriptomics analysis**, specifically designed to interpret **cell-type distribution patterns** produced by spatial transcriptomics cell-type deconvolution tools such as [Cell2location](https://github.com/BayraktarLab/cell2location).  And then we can use these spots to run Cell-Cell Communication analysis or other analyses.
It classifies the spots into three categories: 
1. **Enrich spots** (defined as the specific celltype spots dominated by a single cell type)
2. **Mix spot** (containing multiple cell types)
3. **Unknown** (unclassifiable spots)

using a series of adaptive thresholds for categorization.

Input data:  
A cell-type distribution patterns result(pd.DataFrame): A DataFrame, is the result produced by Spatial transcriptomics cell-type deconvolution tools, where rows represent spots and columns represent cell types, and the value is the number of each cell type in each spot. 

# Tutorial
## Parameter
**cell_counts** (pd.DataFrame): A DataFrame, where rows represent spots and columns represent cell types, and the value is the number of each cell type in each spot. This file must be provided.    
**top_percentage** (float): The top number of corresponding cell types among all Spots is regarded as spot enrichment, with a default value of ``0.05``. This means that only the top 5 % of ranked spots can be considered as ``Enrich spots``.     
**specific_top_percentages** (dict):  a dictionary used to set the number of corresponding cell types of a specific cell type in all spots to be considered enriched spot. Default is ``None``.  
**ratio_threshold** (float): The ratio threshold of the first and second most cells in the spot, used to determine whether it is a Mix spot or a cell type enrich spot. The default is ``1.1``.
## Return
**result**(pd.Series): A Series indexed by spots, with values representing the classification results of each spot (cell type enrichment, Mix, or Unknown).

## Run the Function
In this example, we assume that Microglia are too sparsely distributed, so we need to increase their ‘specific_top_percentages’ value.   
``
result=SpotEnrich(CellCounts, 
top_percentage=0.15,  
specific_top_percentages={'Microglial' : 0.2},  
ratio_threshold=1.1)
``

Now we use the dictionary `{'Microglia': 0.2}` to set the top-percentage cutoff for Microglia to `0.2`.
