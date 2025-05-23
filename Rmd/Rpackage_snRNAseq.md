
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

### 2. Load and clean sample data

``` r
# Download sample data from https://github.com/songxiaoyu/sCCIgen_data/tree/main/input_data. 

load("Github/sCCIgen_data/input_data/snRNAseq_breast_2025_expr.Rdata")
anno=colnames(expr)

dim(expr)
expr[1:3,1:3]
```

### 3. Analysis of the existing data to provide insights into the parameters of the simulation.

Users can split the simulation into into <u> model fitting </u> and <u>
simulation loading fitted model </u> steps to expedite the simulations,
especially if it runs for more than one times.

Specifically, the `sCCIgen` simulation includes two steps: (1)
estimating model parameters from the input data, and (2) simulating
datasets using the estimated and user-provided parameters. The first
step can be time consuming, especially when the input data is large. By
splitting the simulation into multiple steps, users can estimate model
parameters only once and save the results for future use. This is
especially helpful if users may need to run simulations based on the
same input data for multiple times.

#### Task 1: Estimate model parameters from the snRNAseq for simulation.

This is part is to fit the expression data. When sim_method==“copula”,
it will fit both the gene marignal distribution and gene-gene
correlation. When sim_method==“ind”, it will only fit the gene marginal
distribution.

Note: If the number of genes or cells are large, model fitting may take
some time. It is suggested to selecte a reasonably large sample size per
cell type before the model fitting.

``` r

# model fitting 
ModelEst=Est_ModelPara(expr=expr, anno=anno, sim_method='copula', ncores=10)
saveRDS(ModelEst, file="Github/sCCIgen_data/real_data_est/snRNAseq_est/snRNAseq_breast_2025_fit_w_cor.RDS")
```

Note: Additional tasks are available for simulations with spatial input
(not not snRNAseq based) [here](Rpackage_SRT.md).

### 4. Develop a parameter file

Users need to develop a parameter file. The sample parameter file for
snRNAseq based simulation is
[here](https://github.com/songxiaoyu/sCCIgen_data/tree/main/sample_parameter_file/snRNAseq)
for downloading and filling in to perform simulations.

### 5. Run simulation.

``` r
# load parameter file
input="Github/sCCIgen_data/sample_parameter_file/snRNAseq/scRNAseq_default.tsv"

# The default parameter file does not provide estimated model parameters. 
# Run simulation, with model parameters added in with ModelFitFile.
model_param_path="Github/sCCIgen_data/real_data_est/snRNAseq_est/snRNAseq_breast_2025_fit_w_cor.RDS"
ParaSimulation(input=input, ModelFitFile=model_param_path)

# Run simulation including the estimation of the model parameters.
ParaSimulation(input=input)

```

### 6. Run nested functions to obtain simulation byproducts.

#### Task 1: Plot the spatial regions simulated by sCCIgen.

If users are interested to obtain the simulated regions, a nested
function RandomRegionWindow can be used as folllows:

``` r
# The parameter file specifies that nRegion=2 and seed=1234
win=RandomRegionWindow(nRegion=2, seed=1234)
plot(win$window[[1]], col="pink")
plot(win$window[[2]], col="blue", add=T)
```
