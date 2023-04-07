# PlasBin-flow: comparison with other plasmid binning methods

## Overview

This directory contains the results of various plasmid binning methods which have been included in the our analysis in the main paper.

Sample set (sset): For this experiments, we have divided the samples in our test dataset into 2 sets. The bacterial species in the first set of samples is supported by gplas (using mlplasmids) whereas that in the second dataset are not.

Evaluation types (et): Evaluations have been carried out at `basepair` as well as `contig` level. 

Length thrsholds (th): Finally, we evaluated the results for various minimum contig length thresholds (`0, 100, 1000`). Contigs below these length threshold have been left out of the respective evaluations.

## Obtaining statistics

The scripts to obtain and visualize the statistics have been included in the `scripts` subdirectory. 

We first save the statistics for various combinations of (sset, et, th) in the `eval` subfolder in a CSV file for each combination:
```
python scripts/eval_to_file.py
```

We then visualize the performance of various binning methods:
```
python scripts/create_fig.py -sset X -et EVAL_TYPE -th Y 
```
where `X` is `1` or `2` based on set of samples, `EVAL_TYPE` refers to the desired type of evaluation mentioned above and `Y` is the desired length threshold. 