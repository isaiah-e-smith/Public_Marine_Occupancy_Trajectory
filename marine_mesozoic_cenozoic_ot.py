#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sun Jun 12 20:47:25 2022

@author: isaiahsmith

"""

# Isaiah E. Smith
# Last updated: 12 June 2022
# pre-processing of data for input into models etc.
# If you have questions about this project, or would like to see other code relating to this project, please reach out to me.


#%%

import pandas as pd
import math
import numpy as np


#%%

# Read in data 

ori_dat = pd.read_csv('py_prep_marine_mesozoic_cenozoic.csv')

print(ori_dat.head())

#%%

print("Number of columns:",len(ori_dat.columns), "\n")

print(ori_dat.columns)

    
#%%

print("Number of rows:",len(ori_dat), "\n")


#%%

# Drop rows that are missing paleocoordinate data, age data, and genus/species data, as these will be useless for this analysis 
# Also reset index 

ori_dat=ori_dat.dropna(subset=["Estimated.Paleo.Longitude"])
ori_dat=ori_dat.dropna(subset=["Estimated.Paleo.Latitude"])
ori_dat=ori_dat.dropna(subset=["Age..Ma..Gradstein.et.al..2012"])
ori_dat=ori_dat.dropna(subset=["Genus"])
ori_dat=ori_dat.dropna(subset=["Species"])
ori_dat = ori_dat.reset_index()

print("Number of rows, after removing rows with missing data:",len(ori_dat), "\n")


#%%

# Now take the species column and the genus column and make a combined column

gen_sp = []

for i in range(0, len(ori_dat)):
    tmp_sp = ori_dat['Species'][i]
    tmp_gen = ori_dat['Genus'][i]
    tmp_list = [tmp_gen, tmp_sp]
    tmp_gen_sp = ("_".join(tmp_list))
    gen_sp.append(tmp_gen_sp)

ori_dat["gen_sp"] = gen_sp 


#%%

print(ori_dat.head())

print("Number of columns:",len(ori_dat.columns), "\n")

print(ori_dat.columns)


#%%

print("Number of rows:",len(ori_dat), "\n")


#%%

# Assign each occurence to a time bin 
# First, select bin size 
# Then, calculate bins. 

binsize = 1000000

bin_cuts = []

# create bin edges
for i in np.arange(math.floor(min(ori_dat["Age..Ma..Gradstein.et.al..2012"])), math.ceil(max(ori_dat["Age..Ma..Gradstein.et.al..2012"])), (binsize/1000000)):
    bin_cuts.append(i)
    
print("Number of bins:",len(bin_cuts), "\n")

# bin age data 
bin_dat = np.digitize(ori_dat["Age..Ma..Gradstein.et.al..2012"],bin_cuts,right=True)

print("Number of bin assingments:",len(bin_dat), "\n")

#%%

bin = []

for i in range(0, len(bin_dat)):
    tmp_bin = bin_dat[i]
    bin.append(tmp_bin)
    
ori_dat['bin'] = bin
    


#%%

print(ori_dat.head())

print("Number of columns:",len(ori_dat.columns), "\n")

print(ori_dat.columns)


#%%

# create empty dataframe to fill:

un_factor = pd.DataFrame(columns = ["bin", "taxon", "delta_1", "raw_prop_sampled", "extinction", "number_cells_sampled"])

print(un_factor, sep='\n')

#%%

# create list of unique taxa

un_tax = ori_dat.gen_sp.unique()

# how many unique taxa
print("Number of unique taxa:",str(un_tax.size), "\n")


#%%

count = 0

for i in range(0,len(un_tax)-1):
    
      # specific taxon
      tx = un_tax[i]
      
      # creat list of bins in which that taxon occurs, without repeating bins
      binlist = ori_dat['bin'][ori_dat['gen_sp'].str.contains(tx)].unique()
      
      # sort the list from small to large
      binlist.sort()
        
      # for every occupied bin for each taxon
      for j in range(0,(len(binlist)-1)):
        
        # specific bin from bin list
        bin = binlist[j]
        
        if j != len(binlist)-1:
            prev_bin = binlist[j+1]
        else:
            prev_bin = float("nan")
        
        # all geographic cells occupied (by any species) during this given bin
        a = ori_dat['cells'][(ori_dat['bin']==bin)].unique()
        
        # all geographic cells occupied (by any species) during last occupied bin
        b = ori_dat['cells'][ori_dat['bin']==prev_bin].unique()
    
        
        # now we can work with unique taxon/bin pairings
        
        # what unique cells are occupied by this taxon at this time bin
        un = ori_dat['cells'][(ori_dat['gen_sp'].str.contains(tx)) & (ori_dat['bin']==bin)].unique()
        
        # what unique cells are occupied by this taxon at the last occupied time bin (i.e. skip empty bins and "bridge gap")
        un_before = ori_dat['cells'][(ori_dat['gen_sp'].str.contains(tx)) & (ori_dat['bin']==prev_bin)].unique()
        
        # find raw proportion for this taxon in this bin
        rp = len(un)/len(a)
        
        
        
        # find raw proportion for this taxon in previously occupied bin
        if ((len(un_before)==0) | (len(b)==0)):
          rpb = float("nan")
        else:
          rpb = len(un_before)/len(b)
        
        
        prop = (rp)/(rpb)
        
        if (np.isinf(prop)==True):
          ans = float("nan")
        elif (np.isnan(prop) == True):
          ans = float("nan")
        elif (prop==0):
          ans = float("nan") 
        else:
          ans = prop
        
        
        
        # this is the delta value, which can be saved in the row for each taxon/bin combo
        lans = np.log(ans)
        
        # save to data frame
        un_factor.at[count, 'bin'] = bin
        
        un_factor.at[count, 'taxon'] = tx
        
        un_factor.at[count, 'delta_1'] = lans
        
        un_factor.at[count, 'raw_prop_sampled'] = rp
        
        un_factor.at[count, 'number_cells_sampled'] = len(un)
        
        
        count = count + 1
        
        
#%%

print(un_factor.head())


#%%

# Fill Extinction column (1 = goes extinct, or 0 = does not go extinct)
    
# start by filling the Extinction column with all 0's (default: not going extinct in given bin)

un_factor['extinction'] = 0
    

#%%
   
# Now, for last occurrences (extinctions), assign Extinction value of 1
    
for i in range(0,len(un_factor['taxon'].unique())-1):
      
    #for each unique species, create a list of all of the occurences
    spc = un_factor[(un_factor['taxon'] == un_factor['taxon'].unique()[i])]
    
    # last temporal bin (closest to present)
    spbin = list(spc['bin'])[0]
    
    # taxon name
    sptax = list(spc['taxon'])[0]
      
    # Redefines the Extinction variable as 1 (going extinct) for the first numerical bin occurrence (last temporal occurrence) of the taxon
    un_factor.loc[(un_factor['bin'] == spbin) & (un_factor['taxon'] == sptax), 'extinction'] = 1
        
    
#%%

print(un_factor.tail())


#%%

