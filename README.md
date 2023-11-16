# Hash Browns
[![Open Source Starter Files](https://github.com/nrminor/hash-browns/actions/workflows/open-source-starter.yaml/badge.svg)](https://github.com/nrminor/hash-browns/actions/workflows/open-source-starter.yaml) [![Open Source Starter Files](https://github.com/nrminor/hash-browns/actions/workflows/open-source-starter.yaml/badge.svg)](https://github.com/nrminor/hash-browns/actions/workflows/open-source-starter.yaml)
Metagenomic read classification well done.

### Overview
At this stage, Hash Browns is simply a parallelized Nextflow refactor of the bbmap `fetchNt.sh` pipeline, which comes with the bbmap source whenever downloaded. It implements a [MinHash](https://en.wikipedia.org/wiki/MinHash) to classify metagenomic sequencing reads based on NCBI's full nucleotide database (nt).
