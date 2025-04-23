process MEGABLAST {
    tag "${sample_id}"

    cpus 4

    input:
    tuple val(sample_id), path(human_virus_contigs), path(blast_db)

    output:
    tuple val(sample_id), path("${sample_id}.txt")

    script:
    def outfmt = "6 qseqid qlen sseqid stitle length pident evalue bitscore sscinames staxids"
    def max_target_seqs = 5
    def task = "megablast"
    """
    export BLASTDB=resources/
    blastn -task ${task} \
    -db ${blast_db} \
    -query ${human_virus_contigs} \
    -num_threads ${task.cpus} \
    -outfmt "${outfmt}" \
    -max_target_seqs ${max_target_seqs} \
    -out .
    """
}

process ANNOTATE_MEGABLAST_RESULTS {

    tag "${sample_id}"

    cpus 1

    input:
    tuple val(sample_id), path(blast_txt)

    output:
    tuple val(sample_id), path("${sample_id}.txt"), emit: hits
    path "*.sqlite", emit: taxa_sqlite

    script:
    // TODO: This script is hardcoded only to work with snakemake
    """
    annotate_blast_results.py
    """
}

process FILTER_NON_VIRUS_MEGABLAST_NODES {

    tag "${sample_id}"

    cpus 1

    input:
    tuple val(sample_id), path(annotated_blast)

    output:
    tuple val(sample_id), path("${sample_id}.txt")

    script:
    // TODO: This script is hardcoded only to work with snakemake
    """
    filter_non_virus_megablast_nodes.py
    """
}

process REMOVE_MEGABLAST_MAPPED_CONTIGS {

    tag "${sample_id}"

    cpus 1

    input:
    tuple val(sample_id), path(filtered_megablast_nodes), path(human_virus_contigs)

    output:
    tuple val(sample_id), path("${sample_id}.classified.txt"), path("${sample_id}.pruned.fasta")

    script:
    // TODO: This script is hardcoded only to work with snakemake
    """
    remove_megablast_mapped_contigs.py
    """
}

process BLASTN_CLASSIFY {
    tag "${sample_id}"

    cpus 4

    input:
    tuple val(sample_id), path(classified), path(pruned_contigs), path(blast_db)

    output:
    tuple val(sample_id), path("${sample_id}.txt")

    script:
    def outfmt = "6 qseqid qlen sseqid stitle length pident evalue bitscore sscinames staxids"
    def max_target_seqs = 5
    def task = "blastn"
    """
    export BLASTDB=resources/
    blastn -task ${task} \
    -db ${blast_db} \
    -query ${pruned_contigs} \
    -num_threads ${task.cpus} \
    -outfmt "${outfmt}" \
    -max_target_seqs ${max_target_seqs} \
    -out .
    """
}

process ANNOTATE_BLASTN_RESULTS {

    tag "${sample_id}"

    cpus 1

    input:
    tuple val(sample_id), path(blastn_txt), path(taxa_sqlite)

    output:
    tuple val(sample_id), path("${sample_id}.txt")

    script:
    // TODO: This script is hardcoded only to work with snakemake
    """
    annotate_blast_results.py
    """
}

process FILTER_NON_VIRUS_BLASTN_NODES {

    tag "${sample_id}"

    cpus 1

    input:
    tuple val(sample_id), path(annotated_blast)

    output:
    tuple val(sample_id), path("${sample_id}.txt")

    script:
    // TODO: This script is hardcoded only to work with snakemake
    """
    filter_non_virus_megablast_nodes.py
    """
}

process EXTRACT_UNCLASSIFIED_CONTIGS {

    tag "${sample_id}"

    cpus 1

    input:
    tuple val(sample_id), path(contigs), path(megablast_txt), path(blastn_txt)

    output:
    tuple val(sample_id), path("${sample_id}.unclassified.fasta")

    script:
    // TODO: This script is hardcoded only to work with snakemake
    """
    extract_unclassified_contigs.py
    """
}
