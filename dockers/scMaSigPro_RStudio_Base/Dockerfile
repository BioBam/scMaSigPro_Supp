# Use rocker/rstudio as the base image
FROM rocker/rstudio:4.3

# Build Command
## docker build --build-arg GITHUB_PAT=PAT_HERE -t spriyansh29/sc_masigpro:rstudio_base .

# Maintainer
LABEL maintainer = Priyansh Srivastava <spriyansh29@gmail.com>

# Install Required OS libs
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libpng-dev \
    gdal-bin \
    libgdal-dev \
    libfontconfig1-dev \
    libudunits2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libglpk40 \
    texlive \
    texlive-latex-extra \
    texlive-fonts-recommended \
    libxt6 \
    libcairo2-dev \
    patch \
    libgeos-dev \
    python \
    python3 \
    python3-pip \
    libxt-dev \
    cmake \
    libgsl-dev \
    build-essential \
    libfftw3-dev \
    git \
    pip \
    pip3 \
    libqpdf28 \
    qpdf \
    libmagick++-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install Leiden Lang
RUN pip install leidenalg
RUN pip install phate
RUN pip install pydiffmap

# Installation in all in one
RUN Rscript -e 'install.packages(c("reticulate", "diffusionMap", "igraph", "terra", "intergraph", "leidenbase", "tidyverse", "BiocManager", "metap", "intergraph", "spdep", "eHOF", "devtools", "phateR", "doParallel", "MatrixExtra", "ggpubr", "ggrastr", "Cairo", "cellranger", "mutoss", "rvest", "xml2", "lava", "qqconf", "revealjs", "paletteer", "ggpubr", "markdown", "magick", "moments","entropy","venn","mclust","UpSetR", "coop", "styler", "ggupset", "mclust", "ComplexUpset", "RColorConesa", "wesanderson"), repos = "https://cloud.r-project.org/", clean = TRUE, dependencies = TRUE, Ncpus = 24, quiet = TRUE)'

# BioConductor
RUN Rscript -e "BiocManager::install('qvalue', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('multtest', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('BiocGenerics', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('DelayedArray', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('DelayedMatrixStats', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('limma', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('S4Vectors', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('SingleCellExperiment', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('SummarizedExperiment', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('batchelor', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('Matrix.utils', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('rhdf5', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('HDF5Array', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('slingshot', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('miloR', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('tradeSeq', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('scater', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('glmGamPoi', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('SingleR', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('splatter', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('DESeq2', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('scran', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('limma', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('org.Hs.eg.db', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('GO.db', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('BiocStyle', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('BiocCheck', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('maSigPro', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('iCOBRA', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('clusterProfiler', update = TRUE, ask = FALSE, Ncpus = 24)" && \
    Rscript -e "BiocManager::install('biomaRt', update = TRUE, ask = FALSE, Ncpus = 24)"

# Receive GitHub PAT as a build argument
ARG GITHUB_PAT

# Receive GitHub PAT as a build argument
ARG GITHUB_PAT

# Install specific version of Seurat
RUN Rscript -e "Sys.setenv(GITHUB_PAT = '${GITHUB_PAT}'); \
                remotes::install_github('satijalab/seurat', 'seurat5', quiet = TRUE, Ncpus = 24)"

# Install seurat-data from satijalab on GitHub
RUN Rscript -e "Sys.setenv(GITHUB_PAT = '${GITHUB_PAT}'); \
                remotes::install_github('satijalab/seurat-data', quiet = TRUE, Ncpus = 24)"

# Install azimuth from satijalab on GitHub
RUN Rscript -e "Sys.setenv(GITHUB_PAT = '${GITHUB_PAT}'); \
                remotes::install_github('satijalab/azimuth', quiet = TRUE, Ncpus = 24)"

# Install seurat-disk from mojaveazure on GitHub
RUN Rscript -e "Sys.setenv(GITHUB_PAT = '${GITHUB_PAT}'); \
                remotes::install_github('mojaveazure/seurat-disk')"

# Install monocle3 from cole-trapnell-lab on GitHub
RUN Rscript -e "Sys.setenv(GITHUB_PAT = '${GITHUB_PAT}'); \
                remotes::install_github('cole-trapnell-lab/monocle3')"

# Fix irlba and Matrix
RUN Rscript -e "remove.packages('irlba'); \
                install.packages(c('Matrix', 'irlba'), repos = 'https://cloud.r-project.org/', \
                clean = TRUE, dependencies = TRUE, Ncpus = 24, quiet = TRUE)"

# Create Users
RUN mkdir /supp_data
RUN mkdir /additional_data
RUN mkdir /admin
RUN mkdir /priyansh
RUN mkdir /priyansh/gitRepos
RUN mkdir /guest
RUN mkdir /esr

# Expose port 8787 (RStudio runs on this port)
EXPOSE 8787

# Init
CMD ["/init"]

## How to Run
#docker run -p 8003:8787 --rm -e USERID=$(id -u) -e USER=admin -e PASSWORD=admin spriyansh29/sc_masigpro:rstudio_base