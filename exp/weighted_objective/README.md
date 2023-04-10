# PlasBin-flow: changing weights in the objective function

## Overview

This directory contains the results of different weighting schemes used for the PlasBin-flow objective function. 

Evaluation types (et): Evaluations have been carried out at `basepair` as well as `contig` level. 

Length thresholds (th): Finally, we evaluated the results for various minimum contig length thresholds (`0, 100, 1000`). Contigs below these length threshold have been left out of the respective evaluations.

## Obtaining statistics

The scripts to obtain and visualize the statistics have been included in the `scripts` subdirectory. 

We first save the statistics for various combinations of (et, th) in the `eval` subfolder in a CSV file for each combination:
```
python scripts/eval_to_file.py
```

We then visualize the performance of various binning methods:
```
python scripts/create_fig.py -et EVAL_TYPE -th Y 
```
where `EVAL_TYPE` refers to the desired type of evaluation mentioned above and `Y` is the desired length threshold. 