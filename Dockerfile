# Use Ubuntu 20.04 LTS as the base image
FROM ubuntu:20.04

# Set working directory
WORKDIR /scratch

# Set the maintainer label
LABEL maintainer="nrminor@wisc.edu"

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=America/New_York
ENV HOME=/opt
ENV ~=/opt

# Set default command to be the bash shell
ENTRYPOINT ["bash"]

# run a few apt installs
RUN apt-get update && \
    apt-get install -y curl wget && \
    rm -rf /var/lib/apt/lists/* && \
    mkdir /dependencies && \
    dpkg -l > /dependencies/apt-get.lock

# Install everything else with Pixi:
# --------------------------------
# 1) copy the required dependency and configuration file into the image
COPY pyproject.toml $HOME/pyproject.toml
COPY pixi.lock $HOME/pixi.lock

# 2) install pixi
RUN cd $HOME && PIXI_ARCH=x86_64 curl -fsSL https://pixi.sh/install.sh | bash

# 3) make sure pixi and pixi installs are on the $PATH
ENV PATH=$PATH:$HOME/.pixi/bin

# 4) install everything else (save for a crates.io rust dependency) with pixi
RUN cd $HOME && pixi install --frozen && pixi clean cache --assume-yes

# 5) install the rust dependency with cargo, which will be slow because it involves compilation
RUN cd $HOME && pixi run install-scidataflow

# 6) Make sure the cargo-installed binaries are on $PATH
ENV PATH=$PATH:$HOME/.cargo/bin

# 7) Set up the global scidataflow configuration
RUN sdf config --name "Docker Bot" --email "nrminor@wisc.edu" --affiliation "University of Wisconsin - Madison"

# 7) modify the shell config so that each container launches within the pixi env
RUN echo "export PATH=$PATH:$HOME/.pixi/envs/default/bin" >> $HOME/.bashrc

# 8) modify some nextflow environment variables
RUN echo "export NXF_CACHE_DIR=/scratch" >> $HOME/.bashrc
RUN echo "export NXF_HOME=/scratch" >> $HOME/.bashrc

