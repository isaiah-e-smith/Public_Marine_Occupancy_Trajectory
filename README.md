# Public_Marine_Occupancy_Trajectory

## Original code: April 2022
## Minor updates: October 2023

### Background
Selected code from the marine organism occupancy trajectory analysis project. Although geographic range size is widely acknowledged as a predictor of extinction risk in marine species, this relationship is not well quantified. Additionally, the role of change in geographic range through time is even less well-quantified. This project aims to create a better understanding of these spatio-temporal relationships and how they relate to extinction risk in marine organisms. This script simply preprocesses data for future analyses.

### Analysis
This script preprocesses mesozoic marine fossil data for later spatio-temporal analyses. Namely, it assigns each individual fossil occurrence to a time bin (bin size = 1 Ma) and calculates proportional occupancy of geographic cells (based on paleo-geographic coordinates). Proportional geographic cell occupancy is one way that historical geographic species range can be reconstructed from the fossil record. Additionally, change in proportional geographic occupancy is also calculated (delta = 1 time bin). To correct for sampling biases (which are common in the fossil record; Vilhena and Smith, 2013) a "paired-cell approach" (as used in Kiessling and Kocsis, 2016) is used when calculating proportional geographic occupancy and the change in proportional geographic occupancy. 

The final product is a data frame that contains one unique row for each species/time-bin pairing. Each row contains the proportional occupancy, the change in proportional occupancy going from the previous bin to the current bin, and a binary extinction value (0 or 1). This data set can then be used to train a logistic model (or carry out other paleontological analyses). 

### Conclusions
This is only a preprocessing script, so there are no real "conclusions" to be drawn from this script alone. This project (including many unpublished scripts) is currently under peer review/revision for publication. The entirety of the code will be made public upon publication of the associated research paper. In the meantime, if you have questions about this project, or would like to see other code relating to this project, please reach out to me directly.

### Sources
Kiessling, W., & Kocsis, Á. T. (2016). Adding fossil occupancy trajectories to the assessment of modern extinction risk. Biology letters, 12(10), 20150813.   
Vilhena, D. A., & Smith, A. B. (2013). Spatial bias in the marine fossil record. PLoS One, 8(10), e74470.   
