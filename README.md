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
**cell_counts** (pd.DataFrame): A dataframe, where rows represent spots and columns represent cell types, and the value is the number of each cell type in each spot. This file must be provided.    
**top_percentage** (float): The top number of corresponding cell types among all Spots is regarded as spot enrichment, with a default value of ``0.05``. This means that only the top 5 % of ranked spots can be considered as ``Enrich spots``.     
**specific_top_percentages** (dict):  a dictionary used to set the number of corresponding cell types of a specific cell type in all spots to be considered enriched spot. Default is ``None``.  
**ratio_threshold** (float): The ratio threshold of the first and second most cells in the spot, used to determine whether it is a Mix spot or a cell type enrich spot. The default is ``1.1``.
## Return
**result**(pd.DataFrame): A dataframe indexed by spots, with values representing the classification results of each spot (cell type enrichment, Mix, or Unknown).

## Run SpotEnrich
In this example, we assume that we have a cell-type distribution patterns dataframe `CellCounts` and microglial are too sparsely distributed, so we need to increase their `specific_top_percentages` value.   
In `SpotEnrich v1.1.0`, we add new parameter include `specific_ratio_thresholds` and `sample_column`.    
We first use the dictionary `{'Microglia': 0.2}` to set the top-percentage cutoff for Microglia to `0.2`.    
Secondly, we know that bergmann glial cells (BG) and purkinje cells both reside in the purkinje layer, so we set the BG cell-ratio threshold to `1.05` with `{'BG': 1.05}`.        
Finally, `library_id` is used to identify enriched spots per sample.     
Make sure the `SpotEnrich.py` is in your python current working directory, then we can import it.     
```
from SpotEnrich import SpotEnrich
specific_top_percentages={'Microglial':0.2}
specific_ratio_thresholds={'BG':1.05}
result=SpotEnrich(CellCounts, 
                  top_percentage=0.15,  
                  specific_top_percentages={'Microglial' : 0.2},
                  specific_ratio_thresholds=specific_ratio_thresholds,
                  ratio_threshold=1.1,
                  sample_column='library_id'
)
result
```
|      | barcode  | Celltype | library_id |
|  ----  | ----  |  ----  | ----  |
| 13685861458750  | 50465861458750  | Unknown  | control |
| 13685861458751  | 50465861458751  | Mix | treatment |
| 13685861458752  | 50465861458752  | BG | treatment |
| 13685861458752  | 50465861458753  | Microglial | control |

```
# save the result
result=pd.read_csv('EnrichSpot.csv', index_col =0)
```

