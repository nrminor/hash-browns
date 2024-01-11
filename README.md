# Hash Browns
[![Open Source Starter Files](https://github.com/nrminor/hash-browns/actions/workflows/open-source-starter.yaml/badge.svg)](https://github.com/nrminor/hash-browns/actions/workflows/open-source-starter.yaml) [![Docker CI](https://github.com/nrminor/hash-browns/actions/workflows/docker-image.yaml/badge.svg)](https://github.com/nrminor/hash-browns/actions/workflows/docker-image.yaml)

Metagenomic viral read classification well done.

### Overview
Hash-Browns* implements three tools for [MinHashing]((https://en.wikipedia.org/wiki/MinHash)) and frequency-weighting metagenomic classifications: [BBSketch](https://www.biostars.org/p/234837/) (from the `bbmap` suite, sometimes called BBTools), [`sylph`](https://github.com/bluenote-1577/sylph), and [`sourmash`](https://github.com/sourmash-bio/sourmash). Users can choose to use any or none of these tools to classify metagenomic sequencing reads based on the NCBI _nt_ database. At this stage, BBSketch is the most stable and requires the least RAM, though this is subject to change in the near term. 

A few points to take note of:
1. Hash-Browns allows users to download _nt_ and the associated taxonomic information from NCBI and store it locally for use later. This means users can rapidly classify samples once the NCBI datasets are downloaded. But it also means users must have 300-500 GB of disk space available for long-term caching. We use a 2TB external SSD to store the NCBI datasets alongside any sketches thereof, but users may choose whichever solution is best for them.
2. By default, Hash-Browns assumes that all the software it needs is installed and available locally on the machine it's run on. However, we have also prepared a Docker image that can be pulled from Docker Hub and used to run containers without installing software locally. To use containers, simply add `-profile docker` to your Nextflow run command. 
   - Alternatively, if you're using a MacOS machine, a recipe in the provided `justfile` will install all Hash-Browns dependencies locally so you can use it without Docker overhead. To invoke the recipe, make sure [the Rust toolchain](https://www.rust-lang.org/tools/install) is installed, run `cargo install just`, and then run `just mac` in the same directory as the `justfile`. Read more about `just` at [https://just.systems/](https://just.systems/).
3. **Hash-Browns currently only supports Oxford Nanopore long reads**. We recommend users look into [nf-core/taxprofiler](https://nf-co.re/taxprofiler) or [nf-core/mag](https://nf-co.re/mag) for short read applications, though we may add support for short reads in the future.

* - If you must know, HASH-BROWNS is short for **H**igh-throughput **A**mplified **S**equence **H**ash **B**ioinformatics for **R**ealtime **O**xford **N**anopore **W**eighted **N**o-bias **S**ketching
