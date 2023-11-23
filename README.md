# Supplementary Materials for [scMaSigPro: Detection of differential genes along pseudotime with MaSigPro]

This repository contains the supplementary material for the [scMaSigPro: Detection of differential genes along pseudotime with MaSigPro], which details our bioinformatics analysis. The repository is structured as follows:

## Directory Structure
- `Analysis_Public_Data`: This directory contains scripts and data related to the public datasets we analyzed i.e. [Setty et al, 2019](https://www.nature.com/articles/s41587-019-0068-4). 
- `Article_Image`: Images used in the main article are stored here.
- `benchmarks`: Benchmark tests for our analysis scripts can be found in this directory.
- `Figures`: Generated figures from our analysis that are included in the article.
- `LICENSE`: The license file for the code and data used in this project.
- `R_Scripts`: R scripts used for statistical analysis and data processing.
- `README.md`: The file you're currently reading, which explains the repository structure and how to use the files.
- `References`: Bibliographic data and reference materials used in the article.
- `RMDs`: R Markdown files used to generate the report and figures.
- `scMaSigPro_Supp.Rproj`: R Project file for supplementary data analysis.
- `Tables`: Tables generated from the analysis that are included in the article.

## Docker Container
We also provide a Docker container to facilitate the replication of our analysis environment. The conatiner is avaible here [spriyansh29/sc_masigpro_supp](https://hub.docker.com/repository/docker/spriyansh29/sc_masigpro/general). This ensures that you can run our scripts with all the necessary dependencies already installed. Please make sure that you have docker already installed on your system.

This image is using *[rocker/rstudio:4](https://hub.docker.com/layers/rocker/rstudio/4/images/sha256-f8c7260993558a5683ae6874c602233c6ceff962b351476625e149ae38a3a41e?context=explore)* as the base image. Please install and verify that you have docker up and running, visit [docker installation](https://docs.docker.com/engine/install/) for more details.

### Pull the image
Depending on the CPU architecture pull the one which is compatible with yours using the following commands.

1. "AMD64" (Most cases)
```
docker pull spriyansh29/sc_masigpro_supp:amd64
```

2. "ARM64" (Apple M Chips )
```
spriyansh29/sc_masigpro_supp:arm64
```

### To use the Docker container, follow these steps:

Run the following command in the command line or powershell ()windows
```
docker run -p 8888:8787 -d --rm -e USERID=$(id -u) \
    -v path/to/local/directory:/supp_data \
    spriyansh29/sc_masigpro:amd64
```

Replace `path/to/local/directory` with the path to the directory on your host system.

4. After running the command, RStudio should be accessible via your web browser at `http://localhost:8888`. Both the username and password are `admin`.

5. Once inside the RStudio interface, open the new project located at `scMaSigPro_Supp/scMaSigPro_Supp.Rproj`.

6. Now you should be able to run the scripts.

## Contributing

Please refer to `LICENSE` for the terms of use before you use or contribute to this repository.

## Questions

If you have any questions or need further clarification, please file an issue in this repository, and we will get back to you as soon as possible.

## Citation

If you use the materials from this repository, please cite our article as follows: 