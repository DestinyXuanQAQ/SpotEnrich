'''
# Find the celltype enrich spot
# It can find the celltype enrich spot after cell2location, and these spots can use to analysis cell-cell communication.
# Date:    2025-7-14
'''

import pandas as pd 
import scanpy as sc
import os
import numpy

def cell_type_enrichment_spot_analysis(cell_counts, top_percentage=0.05, specific_top_percentages=None, ratio_threshold=1.1): 
    """ 
    Get specific cell type enrich spots based on the results of Cell2location.
 
    Parameter: 
    cell_counts (pd.DataFrame): A DataFrame, where rows represent spots and columns represent cell types, and the value is the number of each cell type in each spot. 
    top_percentage (float): The top number of corresponding cell types among all Spots is regarded as spot enrichment, with a default value of 0.05. 
    specific_top_percentages (dict):  a dictionary used to set the number of corresponding cell types of a specific cell type in all spots to be considered enriched spot. The default is None.
    ratio_threshold (float): The ratio threshold of the first and second most cells in the spot, used to determine whether it is a Mix spot or a cell type enrich spot. The default is 1.1.
 
    Return: 
    pd.Series: A Series indexed by spots, with values representing the classification results of each spot (cell type enrichment, Mix, or Unknown). 
    """ 
    results = [] 
    for index, row in cell_counts.iterrows():  
        # Find out the cell type with the largest number of cells and its quantity in each spot 
        sorted_row = row.sort_values(ascending=False)  
        first_cell_type = sorted_row.index[0]  
        first_count = sorted_row.iloc[0]  
        second_count = sorted_row.iloc[1]  
 
        #  
        ratio = first_count / second_count if second_count > 0 else np.inf  
 
        # Calculate the ratio of the number of the first to the second most numerous cells
        if specific_top_percentages and first_cell_type in specific_top_percentages: 
            current_top_percentage = specific_top_percentages[first_cell_type] 
        else: 
            current_top_percentage = top_percentage 
 
        # Determine whether it is an enriched spot
        threshold = cell_counts[first_cell_type].quantile(1 - current_top_percentage) 
        if ratio > ratio_threshold and first_count >= threshold: 
            results.append(first_cell_type)  
        elif 1 <= ratio <= ratio_threshold: 
            results.append('Mix')  
        else: 
            results.append('Unknown')  
 
    return pd.Series(results, index=cell_counts.index)  
