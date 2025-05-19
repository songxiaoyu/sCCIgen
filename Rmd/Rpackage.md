
## Tutorial for simuation with `sCCIgen` R package (without the interactive interface).

### 1. Download R package

While R docker with a user interface is attractive, experienced users
might like the flexibility of R pacakge for running specialized tasks.

You can install the latest version directly from GitHub with
[devtools](https://github.com/hadley/devtools):

``` r
install.packages("devtools")
devtools::install_github("songxiaoyu/sCCIgen")
library(sCCIgen)
# If the built-in sCCIgen data is needed for your simulation, load this package as well.
devtools::install_github("songxiaoyu/sCCIgen_data")
library(sCCIgen_data)
```

### 2. Develop a parameter file

The sets of parameters needed for simulation depend on multiple factors
such as the input data (e.g. single-cell RNAseq vs SRT), the purpose of
the simulation (e.g. adding de novo patterns or not), and the output
(e.g. number of simulated data sets, single-cell vs multi-cell
resoultion).

It is highly recommended to use <u> the interactive interface </u> for
at least once to get the parameter file structure based on your
simulation situations.

Alternatively, users can use the sample files
[here](sample_parameter_file) to fill in the parameters of interest for
their simulation.

Last but not least, users can use the help manual for each function to
construct/revise the parameter file from the scratch.

### 3. Perform the entire simulation and save the results.

Assuming you already have a parameter file, you can run the entire
simulation using codes like this:

``` r
input="PathToParameterFile"
ParaSimulation(input=input)
```

### 4. Split the simulation to multiple steps for expedition.

Alternative to 3, users can split the simulation to multiple steps. The
primary gain is to expedite the simulations, especially if it runs for
more than one times.

Specifically, the `sCCIgen` simulation includes two steps: (1)
estimating model parameters from the input data, and (2) simulating
datasets using the estimated parameters. The first step is usually time
consuming, especially when the input data is large. By splitting the
simulation into multiple steps, users can estimate it only once and save
the results for future use. This is especially helpful if users may need
to run simulations based on the same input data for multiple times .

Therefore, users can split the simulation into <u> model fitting </u>
and <u> simulation loading fitted model </u> steps. Here are the
step-by-step tutorial for the fast simulation:

#### Step 1: Estimate Gaussian Copula (if gene-gene correlation is considered in simulation).

This is part of the model fitting to capture gene-gene correlation.
Users can save this output for for simulating correlated genes.

``` r
# Assuming InputData has been downloaded from sCCIgen_data
load("InputData/cell_feature_data/snRNAseq_breast_cellfeature_033023.RData")
load("InputData/expression_data/snRNAseq_breast_expr_033023.RData")
# Estiamte Gaussian Copula and save results
CopulaEst=Est_GeneCopula(expr=expr, 
                         anno=cell_feature[,1], 
                         min_nonzero_num =3,
                         zp_cutoff=0.8, ncores=10)
save(CopulaEst, file="InputData/copula_data/snRNAseq_breast_Copula_033023.RData")
```

#### Step 2: Fit the marginal models of the gene expressions.

``` r

# Digest parameters
para=ParaDigest(input)

# If gene-gene correlation is not needed in simulation, specify:
# CopulaEst=NA
# Fitted expression model for each gene using the estimated CopulaEst and save results
ModelFitFile=ParaFitExpr(para=para, 
                         expr=expr, 
                         feature=cell_feature,
                         CopulaEst=CopulaEst, ncores=10, save=T)
```

#### Step 3: Simulation loading the fitted models.

``` r
# Run simulation without the need for refitting the expression models.
ParaSimulation(input=input, ModelFitFile=ModelFitFile)
```

### 5. Run nested functions for simulation byproducts.

#### Task 1: Plot the spatial regions simulated by `sCCIgen`.

``` r
win=RandomRegionWindow(nRegion=2, seed=123)
plot(win$window[[1]], col="pink")
plot(win$window[[2]], col="blue", add=T)
plot(win$window[[3]], col="orange", add=T)
```
