
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `sCCIgen`: A high-fidelity spatially resolved transcriptomics data simulator for cell-cell interaction studies.

## 1. Introduction

Spatially resolved transcriptomics (SRT) provides an invaluable avenue
for examining cell-cell interactions within native tissue environments.
The development and evaluation of analytical tools for SRT data
necessitate tools for generating synthetic datasets with known ground
truth of cell-cell interaction induced features.

We introduce `sCCIgen`, a novel real-data-based simulator tailored to
generate high-fidelity SRT data with a focus on cell-cell interactions.
`sCCIgen` preserves transcriptomic and spatial characteristics in SRT
data, while comprehensively models various cell-cell interaction
features, including cell colocalization, spatial dependence among gene
expressions, and gene-gene interactions between nearby cells.

Reference: XS… ““.

## 2. Options

Users have two options to perform SRT simulations using `sCCIgen`.
Depending on preference, Users can choose to use:

- **An interactive interface.** The interface helps users to choose
  tasks, data, and parameters, and lets `sCCIgen` run simulations to
  provide simulated data, together with parameters and documentations,
  saved in the working directory. As simulations usually involve many
  parameters, this option is highly recommended for the first time user
  or researchers with limited programming expertise.
  - [Tutorials for using the interface to simulate SRT based on snRNAseq
    data.](https://github.com/songxiaoyu/sCCIgen/tree/main/Rmd/Interface_snRNAseq.md)
  - [Tutorials for using the interface to simulate SRT based on
    single-cell SRT
    data.](https://github.com/songxiaoyu/sCCIgen/tree/main/Rmd/Interface_SRT.md)
  - [Tutorials for using the interface to simulate SRT based on
    single-cell expression and unpaired spatial
    data.](https://github.com/songxiaoyu/sCCIgen/tree/main/Rmd/Interface_unpaired.md)
- **R package.** A direct use of `sCCIgen` package offers greater
  flexibility, allowing users to adjust more parameters that would
  otherwise be fixed at default values. It also enables the output of
  simulation byproducts, such as spatial windows. This approach is
  recommended for users with goode programming experience who wish to
  perform customized simulations.
  - [Tutorials for using the R package to simulate SRT based on snRNAseq
    data.](https://github.com/songxiaoyu/sCCIgen/tree/main/Rmd/Rpackage_snRNAseq.md)
  - [Tutorials for using the R package to simulate SRT based on
    single-cell SRT
    data.](https://github.com/songxiaoyu/sCCIgen/tree/main/Rmd/Rpackage_SRT.md)
  - [Tutorials for using the R package to simulate SRT based on
    single-cell expression and unpaired spatial
    data.](https://github.com/songxiaoyu/sCCIgen/tree/main/Rmd/Rpackage_unpaired.md)
