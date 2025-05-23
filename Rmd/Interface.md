
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Tutorial for simuation with `sCCIgen` interactive interface.

In R (or R Studio), open the Rshiny interface by running:

``` r
# Set up to your working directory first 
setwd("~path")

# Open the interface
library(sCCIgen)
run_interactive_sCCIgen()
```

It will show you three tasks to choose from:

#### Task 1: I want to download a pre-simulated SRT dataset.

`sCCIgen` provides three pre-simulated datasets. You can quickly check
the data format and explore your analyses before spending efforts on
learning the package and simulating data. These datasets are simulated
based on normal breast snRNAseq, mouse brain SeqFISH+, and human ovarian
cancer MERFISH.

Here, select one example data and it will be downloaded to your working
directory.

#### Task 2: I want to create a parameter file.

If you select Task 2, `sCCIgen` will first walk you through (1) the
selection of reference dataset, and (2) the additional parameter
determination questions to help you generate a parameter file.

Playing with the decoy data may help get familiar with the simulator,
while spending little time waiting for the simulation results. For
simulation based on real datasets, users can also use our built-in real
data by selection (5-7) or upload input data by selecting “user input”.

#### Task 3: I have a parameter file and want to run a simulation.

Finally, users have the parameter file and can conduct the simulations.
Just tell `sCCIgen` the path to the parameter file, and it will provide
the simulated data and documentation!
