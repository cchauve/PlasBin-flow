# PlasBin-flow: experiments

## Overview

This directory contains the scripts to evaluate the performance of plasmid binning methods. 

The `exp` directory contains the evaluations for the two experiments described in the paper *PlasBin-flow: A flow-based MILP algorithm for plasmid contigs binning*. The `exp/main_results` directory contains the results of various plasmid binning methods which have been included in the our analysis in the main paper. The `exp/weighted_objective` directory contains the results of different weighting schemes used for the PlasBin-flow objective function. 

The `evaluations` directory contains the text files with details on statistics for all the samples in our test set. The `basepair_level` and `contig_level` directories contain the weighted and unweighted statistics as described in the paper for three different minimum contig length thresholds (`0, 100, 1000`) each.

The output files for both the experiments can be accessed at https://doi.org/10.5281/zenodo.7807303.

## Obtaining statistics

The scripts to obtain and visualize the statistics for both the experiments have been included in the `scripts` subdirectory for the respective experiment. 