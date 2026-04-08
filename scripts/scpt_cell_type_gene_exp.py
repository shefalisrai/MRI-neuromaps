#Script courtesy of Dr. Golia Shfiei, PennLINCS, University of Pennsylvania
#This script uses the abagen toolbox to get gene expression data for the Schaefer 100 parcellation. 
# The output is a csv file with 100 rows (one for each parcel) and 20,737 columns (one for each gene). 
# The values in the csv file represent the expression level of each gene in each parcel. 


#Modified by SR, April 2026

import numpy as np
import pandas as pd
import abagen

parcellationDir = '/Users/shefalirai/Desktop/UCSD_Research/TemperamentData/Temperament_analysis/neuromaps/parcellationDir/'
outDir = '/Users/shefalirai/Desktop/UCSD_Research/TemperamentData/Temperament_analysis/neuromaps/results/'

schaefer_mni = (parcellationDir + 'Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.nii.gz')

expression = abagen.get_expression_data(schaefer_mni, missing='interpolate',
                                         lr_mirror='bidirectional')

# Save output
expression.to_csv(outDir + 'schaefer100_gene_expression.csv')
print(expression.shape)
print(expression.head())