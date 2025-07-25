'''
# SpotEnrich
# It can classify the celltype enrich spot after cell2location, and these spots can use to analysis cell-cell communication.
# Update:  2025-7-20 
Updated content: 
    1. Support independent ratio thresholds for specific cell types (e.g., neurons need to be >1.3)
    2. Retain the global threshold and the cell type-specific enrichment threshold
    3. Support group analysis by specified sample columns
    4. Add sample information to increase the number of enrich spots for each slice
'''

import pandas as pd 
import numpy as np 
import scanpy as sc
import os

def _process_spot(row, thresholds, specific_ratio_thresholds, ratio_threshold): 
    sorted_row = row.sort_values(ascending=False)  
    top2 = sorted_row.iloc[:2]  
    cell_type1, count1 = top2.index[0],  top2.iloc[0]  
    count2 = top2.iloc[1]  if len(top2) > 1 else 0 
 
    # Determine the ratio threshold for the current cell type
    if specific_ratio_thresholds and cell_type1 in specific_ratio_thresholds: 
        current_ratio_threshold = specific_ratio_thresholds[cell_type1] 
    else: 
        current_ratio_threshold = ratio_threshold 
 
    ratio = count1 / count2 if count2 > 0 else np.inf  
    thresh = thresholds[cell_type1] 
 
    # apply the thresholds to classify the spot 
    if ratio > current_ratio_threshold and count1 >= thresh: 
        return cell_type1 
    elif ratio_threshold >= ratio >= 1: 
        return 'Mix' 
    return 'Unknown' 
 
def SpotEnrich( 
    cell_counts, 
    sample_column=None,  # add: set the sample id information(a colname)
    top_percentage=0.05, 
    specific_top_percentages=None, 
    ratio_threshold=1.1, 
    specific_ratio_thresholds=None, 
): 
    """   
    参数: 
    cell_counts (pd.DataFrame): A DataFrame, where rows represent spots and columns represent cell types, and the value is the number of each cell type in each spot. 
    sample_column(str) : A colname, in order to split the data base on the sample id.  
    top_percentage (float): Percentage threshold for classify enriched spots, default is 0.05.
    specific_top_percentages (dict):  a dictionary used to set the number of corresponding cell types of a specific cell type in all spots to be considered enriched spot. Default is None.
    ratio_threshold (float): The ratio threshold of the first and second most cells in the spot, used to determine whether it is a Mix spot or a cell type enrich spot. Default is 1.1. 
    specific_ratio_thresholds (dict): 特定细胞类型比值阈值 (如 {'Neuron': 1.3}) 
    specific_ratio_thresholds (dict): A dictionary used to set the different ratio threshold for specific cell types. For example: {'Neuron': 1.3}. Default is None.
    """ 
    if sample_column is not None: 
        # split the data based on the sample id
        grouped = cell_counts.groupby(sample_column)  
        results = [] 
        for _, group in grouped: 
            # remove the sample column
            cell_count_group = group.drop(columns=[sample_column])  
            # calculate the quantile thresholds
            thresholds = {} 
            for ct in cell_count_group.columns:  
                if specific_top_percentages and ct in specific_top_percentages: 
                    thresh = cell_count_group[ct].quantile(1 - specific_top_percentages[ct]) 
                else: 
                    thresh = cell_count_group[ct].quantile(1 - top_percentage) 
                thresholds[ct] = thresh 
 
            group_results = [_process_spot(row, thresholds, specific_ratio_thresholds, ratio_threshold) for _, row in cell_count_group.iterrows()]  
            results.extend(group_results)  
    else: 
        # if sample id is None, just process the whole data without splitting
        # remove the sample column
        if sample_column in cell_counts.columns:  
            cell_counts = cell_counts.drop(columns=[sample_column])  
        # calculate the quantile thresholds
        thresholds = {} 
        for ct in cell_counts.columns:  
            if specific_top_percentages and ct in specific_top_percentages: 
                thresh = cell_counts[ct].quantile(1 - specific_top_percentages[ct]) 
            else: 
                thresh = cell_counts[ct].quantile(1 - top_percentage) 
            thresholds[ct] = thresh 
 
        results = [_process_spot(row, thresholds, specific_ratio_thresholds, ratio_threshold) for _, row in cell_counts.iterrows()]  

    # tansfrom to dataframe and set the column name
    result_df = pd.DataFrame({ 
        'barcode': cell_counts.index,  
        'Celltype': results 
    }) 
    result_df.index  = cell_counts.index  
    result_df=result_df.loc[cell_counts.index,:]
    result_df[sample_column]=cell_counts[sample_column]
    return result_df 