
## Tutorial for simuation with `sCCIgen` R package (without the interactive interface) based on SRT data.

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

Users can split the `sCCIgen` simulation into (1) <u> model fitting </u>
and (2) <u> simulation using fitted model and user-provided parameters
</u> steps to expedite the simulations.

It is especially helpful if the number of genes and/or cells are very
large and users want to run simulation for more than once. By splitting
the simulation into these two steps, users can estimate model parameters
only once and save the results for multiple use.

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
saveRDS(ModelEst, file="Github/sCCIgen_data/real_data_est/SeqFishPlus/SeqFishPlusCortex_2025_fit_wo_cor.RDS")
```

#### Task 2: Estimate CCIs in the input data following the existing Giotto pipeline.

``` r

# Preprocess data with Giotto pipeline
db=preprocessGiotto(expr_data=expr, spatial_data=spatial, run_hvg=T, 
                      run_kNN_network=T, run_Delaunay_network=T) 

# Estimate cell-cell attraction and inhibition patterns, and save in pre-defined folder
cellProximityTable(gobject=db, abs_enrichm=0.3, p_adj = 0.05, 
                  save_folder="Github/sCCIgen_data/real_data_est/SeqFishPlus")
                  
# Estimate gene expressions of cells impacted by their neighbors, and save in pre-defined folder                 
ExprDistanceTable(gobject=db, in_hvg=T, region_specific=F, abs_log2fc_ICG=0.25, p_adj = 0.05,
                  save_folder="Github/sCCIgen_data/real_data_est/SeqFishPlus")                 

# Estimate gene expressions of cells impacted by gene expressions of neighboring cells, narrow to 
#  known pairs, such as ligand and receptor pairs, and save in pre-defined folder  

ExprExprTable(gobject=db, database="mouse", region_specific=F, p_adj=0.05, abs_log2fc_LR=0.25,
        save_folder="Github/sCCIgen_data/real_data_est/SeqFishPlus")
```

#### Task 3: Estimate spatial region specific genes.

``` r
# Estimate spatial region specific genes.
SpatialTable(gobject, top_num=2, fdr_cut=0.05, save_folder="Github/sCCIgen_data/real_data_est/SeqFishPlus")
```

### 4. Develop a parameter file

Users need to develop a parameter file. The sample parameter file for
snRNAseq based simulation is
[here](https://github.com/songxiaoyu/sCCIgen_data/tree/main/sample_parameter_file/SRT)
for downloading and filling in to perform simulations.

### 5. Perform the entire simulation and save the results.

Assuming you already have a parameter file, you can run the entire
simulation using codes like this:

``` r
model_param_path="Github/sCCIgen_data/real_data_est/SeqFishPlus/SeqFishPlusCortex_2025_fit_wo_cor.RDS"

# Simulate default data - using existing cells but simulate expression with ground truth
input="Github/sCCIgen_data/sample_parameter_file/SRT/SeqFishPlus_default.tsv"
ParaSimulation(input=input, ModelFitFile=model_param_path)

# Simulate new cells and genes where genes have region specific expressions 
input="Github/sCCIgen_data/sample_parameter_file/SRT/SeqFishPlus_RegionDiffGenes.tsv"
ParaSimulation(input=input, ModelFitFile=model_param_path)

# Simulate new cells and genes where expression of genes are associated with distances to other cells. 
input="Github/sCCIgen_data/sample_parameter_file/SRT/SeqFishPlus_ICGs.tsv"
ParaSimulation(input=input, ModelFitFile=model_param_path)


# Simulate new cells and genes where expression of genes are associated with distances to other cells. 
input="Github/sCCIgen_data/sample_parameter_file/SRT/SeqFishPlus_LR.tsv"
ParaSimulation(input=input, ModelFitFile=model_param_path)

```
