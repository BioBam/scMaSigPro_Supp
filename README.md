# Supplementary Materials for scMaSigPro: Differential Expression Analysis along Single Cell Trajectories

This repository contains the supplementary material for the [scMaSigPro: Differential Expression Analysis along Single Cell Trajectories], which details our evaluation, benchmarks and analysis. 

## Docker

The docker container to re-run all the analysis can be found on [DockerHub](https://hub.docker.com/repository/docker/spriyansh29/sc_masigpro/general)

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
    +--03_
```
