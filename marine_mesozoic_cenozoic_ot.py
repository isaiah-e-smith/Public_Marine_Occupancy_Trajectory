#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sun Jun 12 20:47:25 2022

@author: isaiahsmith

"""

# Isaiah E. Smith
# Last updated: 13 June 2022
# pre-processing of data for input into models etc.
# If you have questions about this project, or would like to see other code relating to this project, please reach out to me.


#%%

import pandas as pd
import math
import numpy as np
from tqdm import tqdm


#%%

# Read in data 

ori_dat = pd.read_csv('py_prep_marine_mesozoic_cenozoic.csv')

print(ori_dat.head())

#%%

print(f"Number of columns: {len(ori_dat.columns)}\n")

print(ori_dat.columns)

    
#%%

print(f"Number of rows: {len(ori_dat)}\n")


#%%

# Drop rows that are missing paleocoordinate data, age data, and genus/species data, as these will be useless for this analysis 
# Also reset index 

ori_dat = ori_dat.dropna(subset=["Estimated.Paleo.Longitude", 
                                 "Estimated.Paleo.Latitude",
                                 "Age..Ma..Gradstein.et.al..2012",
                                 "Genus",
                                 "Species"])

ori_dat = ori_dat.reset_index()

print(f"Number of rows, after removing rows with missing data: {len(ori_dat)}\n")


#%%

# Now take the species column and the genus column and make a combined column

gen_sp = [species + '_' + ori_dat['Genus'][index] for index, species in enumerate(ori_dat['Species'])]

ori_dat["gen_sp"] = gen_sp 


#%%

print(ori_dat.head())

print(f"Number of columns: {len(ori_dat.columns)}\n")

print(ori_dat.columns)


#%%

print(f"Number of rows: {len(ori_dat)}")


#%%

# Assign each occurence to a time bin 
# First, select bin size 
# Then, calculate bins. 

binsize = 1000000 # change depending on preference

min_val = math.floor(min(ori_dat["Age..Ma..Gradstein.et.al..2012"]))
max_val = math.ceil(max(ori_dat["Age..Ma..Gradstein.et.al..2012"]))

# create bin edges
bin_cuts = [i for i in np.arange(min_val, max_val, (binsize/1000000))] # divide by 1000000 since ages are listed in MYA
    
print(f"Number of bins: {len(bin_cuts)}\n")

# bin age data 
bin_dat = np.digitize(ori_dat["Age..Ma..Gradstein.et.al..2012"], bin_cuts, right=True)

print(f"Number of bin assignments: {len(bin_dat)}\n")

#%%

ori_dat['bin'] = bin_dat
    
#%%

print(ori_dat.head())

print(f"Number of columns: {len(ori_dat.columns)}\n")

print(ori_dat.columns)


#%%

# create empty dataframe to fill:

un_factor = pd.DataFrame(columns = ["bin", "taxon", "delta_1", "raw_prop_sampled", "extinction", "number_cells_sampled"])

print(un_factor) 

#%%

# create list of unique taxa

un_tax = ori_dat.gen_sp.unique()

# how many unique taxa
print(f'Number of unique taxa: {len(un_tax)}\n')

#%%
count = 0

for tx in tqdm(un_tax): 
          
      # create list of bins in which that taxon occurs, without repeating bins

      binlist = ori_dat['bin'][ori_dat['gen_sp'] == tx].unique()
      
      # sort the list from small to large
      binlist.sort()

      # for every occupied bin for each taxon
     
      for index, _bin in enumerate(binlist):
        
        if index != len(binlist)-1:
            prev_bin = binlist[index+1]
        else:
            prev_bin = float('nan')

        
        # all geographic cells occupied (by any taxon) during this given bin
        a = ori_dat['cells'][ori_dat['bin'] == _bin].unique()

        # all geographic cells occupied (by any taxon) during last occupied bin
        b = ori_dat['cells'][ori_dat['bin'] == prev_bin].unique()
    
        
        # now we can work with unique taxon/bin pairings
        
        # what unique cells are occupied by this taxon at this time bin
        un = ori_dat['cells'][(ori_dat['gen_sp'] == tx) & (ori_dat['bin'] == _bin)].unique()
      

        # what unique cells are occupied by this taxon at the last occupied time bin (i.e. skip empty bins and "bridge gap")
        un_before = ori_dat['cells'][(ori_dat['gen_sp'] == tx) & (ori_dat['bin'] == prev_bin)].unique()
                
        # find raw proportion for this taxon in this bin
        rp = len(un)/len(a)

        # find raw proportion for this taxon in previously occupied bin
        if ((len(un_before)==0) | (len(b)==0)):
          rpb = float('nan')
        else:
          rpb = len(un_before)/len(b)
        
        
        prop = (rp)/(rpb)
        
        if np.isinf(prop)==True or np.isnan(prop) == True or prop == 0:
          ans = float('nan')
        else:
          ans = prop
        
        
        
        # this is the delta value, which can be saved in the row for each taxon/bin combo
        lans = np.log(ans)

        # save to data frame

        un_factor.at[count, 'bin'] = _bin
        
        un_factor.at[count, 'taxon'] = tx
        
        un_factor.at[count, 'delta_1'] = lans
        
        un_factor.at[count, 'raw_prop_sampled'] = rp
        
        un_factor.at[count, 'number_cells_sampled'] = len(un)
        
        count += 1
        
#%%

print(un_factor.head())


#%%

# Fill Extinction column (1 = goes extinct, or 0 = does not go extinct)
    
# start by filling the Extinction column with all 0's (default: not going extinct in given bin)

un_factor['extinction'] = 0

print(len(un_factor))
    

#%%
   
# Now, for last occurrences (extinctions), assign Extinction value of 1
unique_taxon = un_factor['taxon'].unique()

for taxon in tqdm(unique_taxon):
      
    #for each unique species, create a list of all of the occurences
    spc = un_factor[un_factor['taxon'] == taxon]
    
    # last temporal bin (closest to present)
    spbin = list(spc['bin'])[0]
    
    # taxon name
    sptax = list(spc['taxon'])[0]
      
    # Redefines the Extinction variable as 1 (going extinct) for the first numerical bin occurrence (last temporal occurrence) of the taxon
    un_factor.loc[(un_factor['bin'] == spbin) & (un_factor['taxon'] == sptax), 'extinction'] = 1
        
    
#%%

print(un_factor.tail())


#%%

print(list(un_factor['extinction']))
