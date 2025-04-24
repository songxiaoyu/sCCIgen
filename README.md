
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `sCCIgen` for ‘de novo’ spatial patterns

## 1. Introduction

Existing Spatially Resolved Transcriptomics (SRT) data include limited
spatial patterns. `sCCIgen` generates de novo spatial patterns using
snRNAseq or single-cell SRT as reference.

`sCCIgen` provides both a R package for simulating SRT data and an
interface to guide users through task selection, parameter
specification, simulation, and documentation.

Reference: XS… ““.

## 2. Options

Users have two options to perform SRT simulations using `sCCIgen`.
Depending on preference, you can choose to use:

- **An interactive interface.** The interface helps you to choose tasks,
  data, and parameters, and lets `sCCIgen` run to provide you simulated
  data, with parameters and documentations, saved in your working
  directory. As simulations usually involve many parameters, this option
  is highly recommended for first time user or researchers with limited
  programming expertise.
  - [Tutorials for running the interface.](Interface.md)
- **R package.** A direct use of `sCCIgen` package offers greater
  flexibility, allowing users to adjust more parameters that would
  otherwise be fixed at default values. It also enables the output of
  simulation byproducts, such as spatial windows. This approach is
  recommended for users with goode programming experience who wish to
  perform customized simulations.
  - [Tutorials for running the R package.](Rpackage.md)
