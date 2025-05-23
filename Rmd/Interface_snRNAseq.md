
<!-- README.md is generated from README.Rmd. Please edit that file -->

Run `sCCIgen` with an Interactive Interface on Docker

### 1. Installation

1.  An installation of `Docker Desktop` is highly recommended,
    especially for beginners. Docker can be downloaded at
    <https://www.docker.com>.

2.  `sCCIgen` uses command line prompts, such as through `Terminal` in
    Mac. You might want to download your favorite software to run
    command line. My favorite is Visual Studio Code (`VSC`). If `VSC` is
    used, Docker needs to be installed in its Extension. I will use
    `VSC` for this tutorial.

<!-- 3. Clone `sCCIgen` repository to your local machine, such as on `Terminal` type
```
git clone https://github.com/songxiaoyu/sCCIgen
```
-->

### 2. Run Docker

1.  Open Docker Desktop.

2.  On `VSC` terminal, pull image from the Docker Hub as follows:

<!-- -->

    docker pull songxiaoyu152/st_simulator_test 

3.  Run the Docker container directly with a working directory (WORKDIR)
    bound to your local machine. Note WORKDIR will be the location of
    all of your input data and your outputs.

<!-- -->

    docker run --mount type=bind,source="${WORKDIR}",target=/working_directory -it songxiaoyu152/st_simulator_test

My working directory is
“/Users/songxiaoyu152/Dropbox/SpatialTranscriptomics/Paper_Simulator/UseDocker”,
so this is my code:

    docker run --mount type=bind,source=/Users/songxiaoyu152/Dropbox/SpatialTranscriptomics/Paper_Simulator/UseDocker,target=/working_directory -it songxiaoyu152/st_simulator_test

### 3. Open the interacive interface and choose the tasks of interest.

In R (or R Studio), run

``` r
library(sCCIgen)
run_interactive_sCCIgen()
```

It will show you three tasks to choose from:

#### Task 1: I want to download a pre-simulated spatial transcriptomics dataset.

`sCCIgen` provides three pre-simulated datasets. You can quickly check
the data format and explore your analyses before spending efforts on
learning the package and simulating data. These datasets are simulated
based on normal breast snRNAseq, mouse brain SeqFISH+, and human ovarian
cancer MERFISH.

Here, select one example data and it will be downloaded to your working
directory.

#### Task 2: I want to generate a parameter file using command line prompts.

If you select Task 2, `sCCIgen` will first walk you through a number of
parameter selection questions to help you generate a parameter file. The
first question is to select the data for simulation. Playing with the
decoy data may help you to get familiar with the simulator, while
spending little time on waiting for the simulation results. You can also
use our built-in real data by selection (4-6) or your input data by
selecting “user input”.

#### Task 3: I have a parameter file and want to run a simulation.

Finally, you have the parameter file and can conduct the simulations.
Just tell `sCCIgen` your parameter file, and it will provide you the
data and document you need!
