# Supplementary Materials for scMaSigPro: Differential Expression Analysis along Single Cell Trajectories

This repository contains the supplementary material for the [scMaSigPro: Differential Expression Analysis along Single Cell Trajectories], which details our evaluation, benchmarks and analysis. 

## Docker

The docker container to re-run all the analysis can be found on [DockerHub](https://hub.docker.com/r/spriyansh29/sc_masigpro)

## Directory Structure
1. Directory __'benchmarks':__ Scripts for used to benchmark scMaSigPro can be found in this directory.
```
benchmarks
    |
    +--00_Parameter_Estimation: Splatter script to estimate parameters from real dataset.
    |
    +--01_Sparsity: Evaluation of datasets with increasing Sparsity.
    |
    +--02_Skewness: Evaluation of datasets with heterogeneous distribution of cells along the branching paths.
    |
    +--03_Different_length: Evaluation of datasets with different lengths of the branching paths.

comparison
    |
    +--01_ComparisonWithTradeSeq: Rscripts to compare with tradeSeq.
    |
    +--02_SpeedTimeComparisonWithTradeSeq: Time Comparison with tradeSeq.
    
analysis_public_data
    |
    +--01_Raw_to_processed.r: Raw counts to processed counts.
    |
    +--02_Cell_Type_Annotation.R: Cell Type inference with Azimuth.
    |
    +--03_SubSampling: SubSampling based on cell types.
    |
    +--04_TI_Monocle3.R: Trajectory inference with monocle3.
    |
    +--04.1_ScMaSigPro_Input.R: Selection of branhcing paths.
    |
    +--05_scMaSigPro.R: ScMaSigPro Analysis.
    |
    +--06_Exploration.R: GO enrichment.

Rscripts
    |
    +--combine: Rscripts to combine figures and tables for manuscript.
    |
    +--helper_function: Additional function.
```


