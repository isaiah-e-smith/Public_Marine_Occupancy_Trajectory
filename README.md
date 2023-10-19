# Public_Marine_Occupancy_Trajectory

Selected code from the marine organism occupancy trajectory analysis project. 

This script preprocesses mesozoic marine fossil data for later spatio-temporal analyses. Namely, it assigns each individual fossil occurrence to a time bin and calculates proportional occupancy of geographic cells (based on paleo-geographic coordinates). Proportional geographic cell occupancy is one way that historical geographic species range can be reconstructed from the fossil record.

The final product is a data frame that contains one unique row for each species/time-bin pairing. Each row contains the proportional occupancy, the change in proportional occupancy going from the previous bin to the current bin, and a binary extinction value (0 or 1). This data set can then be used to train a logistic model (or carry out other paleontological analyses). 

If you have questions about this project, or would like to see other code relating to this project, please reach out to me.
