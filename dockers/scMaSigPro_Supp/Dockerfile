# Use rocker/rstudio as the base image
FROM spriyansh29/sc_masigpro:rstudio_base

# Build Command
## $docker build --build-arg GITHUB_TOKEN=PAT_HERE -t spriyansh29/sc_masigpro:supp_amd64_devel .

# Maintainer
LABEL maintainer = Priyansh Srivastava <spriyansh29@gmail.com>

# Change Working directory
WORKDIR "/priyansh/gitRepos"

# Receive GitHub PAT as a build argument
ARG GITHUB_PAT

# Pull the latest version
RUN git clone https://${GITHUB_PAT}@github.com/spriyansh/scMaSigPro_Supp.git

# Change Working directory
WORKDIR /priyansh/gitRepos/scMaSigPro_Supp

# Expose port 8787 (RStudio runs on this port)
EXPOSE 8787

# Init
CMD ["/init"]

## How to Run
#docker run -p 8003:8787 --rm -e USERID=$(id -u) -v /data/priyansh_data/scMaSigPro_supp_data/:/supp_data -e USER=priyansh -e PASSWORD= spriyansh29/sc_masigpro:supp_amd64_devel