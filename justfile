default:
    just --list

# Install the Rust toolchain and the crates used by Amplityper
rust:
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y && \
    cargo install sylph
    cargo install fastqc-rs
    cargo install nanoq
alias rs := rust

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
    r
    -brew install --cask docker
alias brew := homebrew

# Packages to build from source localled (bbtools and centrifuge)
local-builds:
    touch ~/.zprofile
    -mkdir ~/bioinformatics
    wget -q https://sourceforge.net/projects/bbmap/files/latest/download -O ~/bioinformatics/bbmap.tar.gz
    -tar -xzf ~/bioinformatics/bbmap.tar.gz -C ~/bioinformatics
    rm ~/bioinformatics/bbmap.tar.gz

    git clone https://github.com/DaehwanKimLab/centrifuge
    cd centrifuge
    -make
    -sudo make install prefix=/usr/local
    cd ..
    rm -rf centrifuge

