# Use Ubuntu 20.04 LTS as the base image
FROM ubuntu:20.04

# Set working directory
WORKDIR /scratch

# Set the maintainer label
LABEL maintainer="nrminor@wisc.edu"

# Set time zone
ENV TZ America/New_York

# Set environment variables to non-interactive (this prevents some prompts)
ENV DEBIAN_FRONTEND=non-interactive

# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    autoconf \
    automake \
    libtool \
    pkg-config \ 
    zstd \
    pigz \
    wget \
    curl \
    default-jre && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    mkdir /dependencies && \
    dpkg -l > /dependencies/apt-get.lock

# Install VSEARCH
RUN wget https://github.com/torognes/vsearch/archive/v2.25.0.tar.gz && \
    tar xzf v2.25.0.tar.gz && \
    cd vsearch-2.25.0 && \
    ./autogen.sh && \
    ./configure CFLAGS="-O3" CXXFLAGS="-O3" && \
    make -j 8 && \
    make install

# Install SeqKit
RUN wget https://github.com/shenwei356/seqkit/releases/download/v2.1.0/seqkit_linux_amd64.tar.gz && \
    tar -xvzf seqkit_linux_amd64.tar.gz && \
    chmod +x seqkit && \
    mv seqkit /usr/local/bin/ && \
    rm seqkit_linux_amd64.tar.gz

# Install BBMap
RUN wget https://sourceforge.net/projects/bbmap/files/BBMap_38.90.tar.gz && \
    tar -xvzf BBMap_38.90.tar.gz && \
    rm BBMap_38.90.tar.gz

# Set BBMap environment variable
ENV PATH="/opt/bbmap:${PATH}"

# Install csvtk
RUN wget https://github.com/shenwei356/csvtk/releases/download/v0.23.0/csvtk_linux_amd64.tar.gz && \
    tar -xvzf csvtk_linux_amd64.tar.gz && \
    chmod +x csvtk && \
    mv csvtk /usr/local/bin/ && \
    rm csvtk_linux_amd64.tar.gz

# Install Nextflow
RUN curl -s https://get.nextflow.io | bash && \
    chmod +x /scratch/nextflow && \
    mv /scratch/nextflow /usr/local/bin/nextflow

# Set default command
CMD ["bash"]
