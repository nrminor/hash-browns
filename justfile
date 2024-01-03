default:
    just --list

# Install MacOS packages available via homebrew
homebrew:
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)" && \
    -brew install \
    wget \
    curl \
    git \
    zstd \
    pigz \
    unzip \
    make \
    java \
    seqkit \
    csvtk \
    vsearch \
    r \
    nim
    -brew install --cask docker
alias brew := homebrew

# Install the Rust toolchain and the crates used by Amplityper
rust:
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y && \
    cargo install sylph
    cargo install fastqc-rs
    cargo install nanoq
alias rs := rust

# Install Python tools (currently assumes pip is already installed)
python:
    python3 -m pip install sourmash==4.8.4
    python -m pip install sourmash==4.8.4
    pip install sourmash==4.8.4
alias py := python

# Packages to build from source localled (bbtools, seqfu, nextflow, and centrifuge)
local-builds:
    @echo "Installing the bbmap suite."
    touch ~/.zprofile
    -mkdir ~/bioinformatics
    wget -q https://sourceforge.net/projects/bbmap/files/latest/download -O ~/bioinformatics/bbmap.tar.gz
    -tar -xzf ~/bioinformatics/bbmap.tar.gz -C ~/bioinformatics
    rm ~/bioinformatics/bbmap.tar.gz
    chmod +x ~/bioinformatics/bbmap/*
    echo "export PATH=$PATH:~/bioinformatics/bbmap" >> ~/.zprofile
    source ~/.zprofile
    @echo "Installing seqfu"
    git clone https://github.com/telatin/seqfu2 ~/bioinformatics/seqfu2
    (cd ~/bioinformatics/seqfu2 && nimble build)
    echo "export PATH=$PATH:~/bioinformatics/seqfu2/bin" >> ~/.zprofile
    source ~/.zprofile
    @echo "Installing Nextflow"
    wget -qO- https://get.nextflow.io | bash
    -mkdir ~/bioinformatics/nextflow
    mv nextflow ~/bioinformatics/nextflow
    echo "export PATH=$PATH:~/bioinformatics/nextflow" >> ~/.zprofile
    source ~/.zprofile
    @echo "Installing Centrifuge"
    git clone https://github.com/DaehwanKimLab/centrifuge
    cd centrifuge
    -make
    -sudo make install prefix=/usr/local
    cd ..
    rm -rf centrifuge
alias lb := local-builds

all-macos:
    just brew
    just rs
    just py
    just lb
alias mac := all-macos
