
## Tutorial for simuation with `sCCIgen` R package (without the interactive interface) based on SRT data.

### 1. Load R package

``` r
library(sCCIgen)
```

### 2. Load and clean sample data

``` r
# Download sample data from https://github.com/songxiaoyu/sCCIgen_data/tree/main/input_data. 

load("Github/sCCIgen_data/input_data/SeqFishPlusCortex_2025_expr.Rdata")
load("Github/sCCIgen_data/input_data/SeqFishPlusCortex_2025_spatial.Rdata")

dim(expr)
expr[1:3,1:3]
dim(spatial)

anno=colnames(expr)
region=spatial[,4]
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
it will fit both the gene marginal distribution and gene-gene
correlation. When sim_method==“ind”, it will only fit the gene marginal
distribution.

Note: If the number of genes or cells are large, model fitting may take
some time. It is suggested to select a reasonable sample size
(e.g. \<10000 per cell type) before the model fitting, as more cells may
not be needed improve the estimation. Similarly, if some genes are
extremely zero-inflated, narrowing the simulation to reasonbly variable
genes is an option.

``` r

# model fitting 
ModelEst=Est_ModelPara(expr=expr, anno=anno, sim_method='ind', region=region, ncores=14)
saveRDS(ModelEst, file="Github/sCCIgen_data/real_data_est/SeqFishPlus/SeqFishPlusCortex_fit_wo_cor_region.RDS")
```

#### Task 2: Estimate CCIs in the input data following the existing Giotto pipeline.

``` r
library(Giotto)
# Preprocess data with Giotto pipeline
db=preprocessGiotto(expr_data=expr, spatial_data=spatial, run_hvg=T, 
                      run_kNN_network=T, run_Delaunay_network=T) 

# Estimate cell-cell attraction and inhibition patterns, and save in pre-defined folder
cellProximityTable(gobject=db, abs_enrichm=0.3, p_adj = 0.05, 
                  output_file="Github/sCCIgen_data/real_data_est/SeqFishPlus/est_CCI_dist_dist_region.csv")
                  
# Estimate gene expressions of cells impacted by their neighbors, and save in pre-defined folder                 
ExprDistanceTable(gobject=db, in_hvg=T, region_specific=T, abs_log2fc_ICG=0.3, p_adj = 0.05,
                  output_file="Github/sCCIgen_data/real_data_est/SeqFishPlus/est_CCI_dist_expr_region.csv")                 

# Estimate gene expressions of cells impacted by gene expressions of neighboring cells, narrow to 
#  known pairs, such as ligand and receptor pairs, and save in pre-defined folder  

ExprExprTable(gobject=db, database="mouse", region_specific=T, p_adj=0.05, abs_log2fc_LR=0.3,
        output_file="Github/sCCIgen_data/real_data_est/SeqFishPlus/est_CCI_expr_expr_region.csv")
```

#### Task 3: Estimate spatial region specific genes.

``` r
# Estimate spatial region specific genes.
SpatialTable(gobject=db, top_num=2, fdr_cut=0.05, output_file="Github/sCCIgen_data/real_data_est/SeqFishPlus/est_CCI_expr_expr_region.csv")
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
model_param_path="Github/sCCIgen_data/real_data_est/SeqFishPlus/SeqFishPlusCortex_fit_wo_cor_region.RDS"

input="Github/sCCIgen_data/sample_parameter_file/SRT/SeqFishPlus_Region_AllPattern.tsv"
ParaSimulation(input=input, ModelFitFile=model_param_path)

```
