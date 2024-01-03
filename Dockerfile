# Use Ubuntu 20.04 LTS as the base image
FROM ubuntu:20.04

# Set working directory
WORKDIR /scratch

# Set the maintainer label
LABEL maintainer="nrminor@wisc.edu"

# Set time zone
ENV TZ America/New_York

# set home
ENV HOME=/opt
ENV ~=/opt

# Set environment variables to non-interactive (this prevents some prompts)
ENV DEBIAN_FRONTEND=non-interactive

# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    autoconf \
    automake \
    gcc \
    make \
    cmake \
    libtool \
    pkg-config \ 
    zstd \
    pigz \
    unzip \
    wget \
    curl \
    git \
    default-jre \
    python3.12 \
    python3-pip \
    nim && \
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

# Install SeqFu
RUN git clone https://github.com/telatin/seqfu2 && \
    cd seqfu2 && \
    nimble build && \
    cd .. && \
    mv seqfu2 /opt
ENV PATH="${PATH}:/opt/seqfu2/bin"

# Install BBMap
RUN wget https://sourceforge.net/projects/bbmap/files/BBMap_38.90.tar.gz && \
    tar -xvzf BBMap_38.90.tar.gz && \
    rm BBMap_38.90.tar.gz && \
    mv bbmap /opt

# Set BBMap environment variable
ENV PATH="/opt/bbmap:${PATH}"

# Install csvtk
RUN wget https://github.com/shenwei356/csvtk/releases/download/v0.23.0/csvtk_linux_amd64.tar.gz && \
    tar -xvzf csvtk_linux_amd64.tar.gz && \
    chmod +x csvtk && \
    mv csvtk /usr/local/bin/ && \
    rm csvtk_linux_amd64.tar.gz

# install Rust
RUN mkdir -m777 /opt/rust /opt/.cargo
ENV RUSTUP_HOME=/opt/rust CARGO_HOME=/opt/.cargo PATH=/opt/.cargo/bin:$PATH
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | bash -s -- -y && \
    bash "/opt/.cargo/env"

# Install Sylph
RUN git clone https://github.com/bluenote-1577/sylph && \
    cd sylph && \
    cargo install --path . --root /opt/.cargo

# Install fastqc-rs
RUN cargo install fastqc-rs --root /opt/.cargo

# Install Sourmash
RUN pip install sourmash==4.8.4

# Install multiqc
RUN pip install multiqc==1.19

# Install Nextflow
RUN curl -s https://get.nextflow.io | bash && \
    chmod 777 nextflow && \
    mv nextflow /usr/local/bin/ && \
    chmod 777 /usr/local/bin/nextflow

# Set default command
CMD ["bash"]
