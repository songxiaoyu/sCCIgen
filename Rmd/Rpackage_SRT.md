
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

load("Github/sCCIgen_data/input_data/SeqFishPlusCortex_2025_expr.Rdata")
load("Github/sCCIgen_data/input_data/SeqFishPlusCortex_2025_spatial.Rdata")
anno=colnames(expr)

dim(expr)
expr[1:3,1:3]
dim(spatial)
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
ModelEst=Est_ModelPara(expr=expr, anno=anno, sim_method='ind', ncores=10)
saveRDS(ModelEst, file="Github/sCCIgen_data/SeqFISH_plus_est/SeqFishPlusCortex_2025t_2025_fit_wo_cor.RDS")
```

#### Task 2: Estimate CCIs in the input data following the existing Giotto pipeline.

``` r

# Preprocess data with Giotto pipeline
dat=preprocessGiotto(expr_data=expr, spatial_data=spatial, run_hvg=T, 
                      run_kNN_network=T, run_Delaunay_network=T) 

# Estimate cell-cell attraction and inhibition patterns, and save in pre-defined folder
cellProximityTable(gobject=dat, abs_enrichm=0.3, p_adj = 0.05, 
                  save_folder="Github/sCCIgen_data/SeqFISH_plus_est")
                  
# Estimate gene expressions of cells impacted by their neighbors, and save in pre-defined folder                 
ExprDistanceTable(gobject=dat, in_hvg=T, region_specific=F, abs_log2fc_ICG=0.25, p_adj = 0.05,
                  save_folder="Github/sCCIgen_data/SeqFISH_plus_est")                 

# Estimate gene expressions of cells impacted by gene expressions of neighboring cells, narrow to 
#  known pairs, such as ligand and receptor pairs, and save in pre-defined folder  

LRTable(gobject=dat, database="mouse", region_specific=F, p_adj=0.05, abs_log2fc_LR=0.25,
        save_folder="Github/sCCIgen_data/SeqFISH_plus_est")
```

#### Task 3: Estimate spatial region specific genes.

``` r

SpatialTable(gobject, top_num=2, fdr_cut=0.05,
        save_folder="Github/sCCIgen_data/SeqFISH_plus_est")
```

### 4. Develop a parameter file

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

### 5. Perform the entire simulation and save the results.

Assuming you already have a parameter file, you can run the entire
simulation using codes like this:

``` r
input="PathToParameterFile"
ParaSimulation(input=input)
```

### 6. Run nested functions to obtain simulation byproducts.

#### Task 1: Plot the spatial regions simulated by `sCCIgen`.

If users are interested to obtain the simulated regions, a nested
function `RandomRegionWindow` can be used as folllows:

``` r
win=RandomRegionWindow(nRegion=2, seed=123)
plot(win$window[[1]], col="pink")
plot(win$window[[2]], col="blue", add=T)
plot(win$window[[3]], col="orange", add=T)
```
