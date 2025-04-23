process EXTRACT_HUMAN_VIRUS_READS {

    tag "${sample_id}"

    cpus 2

    input:
    tuple val(sample_id), val(sample_type), path(fastq), path(stat_dbss), path(human_virus_taxlist)

    output:
    tuple val(sample_id), val(sample_type), path("*.f*q.gz")

    script:
    """
    seqkit fq2fa --threads 1 ${fastq} \
    | aligns_to \
    -dbss ${stat_dbss} \
    -num_threads ${task.cpus} \
    -tax_list ${human_virus_taxlist} \
    stdin \
    | cut -f1 \
    | seqkit grep -f - \
    ${fastq} \
    -o ${sample_id}.human_virus.fastq.gz
    """
}

process CLASSIFY_CONTIGS_FIRST_PASS {

    tag "${sample_id}"

    input:
   tuple val(sample_id), path(filtered_contigs), path(stat_index)

    output:
    tuple val(sample_id), path(filtered_contigs), path("${sample_id}.firstpass.txt")

    script:
    """
    aligns_to -dbs ${stat_index} \
    -num_threads ${task.cpus} \
    ${filtered_contigs} \
    > "${sample_id}.firstpass.txt"
    """
}

process GENERATE_CONTIGS_TAXA_LIST {

    tag "${sample_id}"

    input:
    tuple val(sample_id), path(filtered_contigs), path(first_pass_file)

    output:
    tuple val(sample_id), path(filtered_contigs), path("${sample_id}.report")

    script:
    """ 
    generate_tax_list.py ${first_pass_file} ${sample_id}.report
    """
}

process CLASSIFY_CONTIGS_SECOND_PASS {

    tag "${sample_id}"

    cpus 1

    input:
    tuple val(sample_id), path(filtered_contigs), path(tax_list), path(stat_dbss)

    output:
    tuple val(sample_id), path("${sample_id}.secondpass.txt")

    script:
    """
    aligns_to \
    -tax_list ${tax_list} \
    -dbss {stat_dbss} \
    -num_threads ${task.cpus} \
    ${filtered_contigs} \
    > "${sample_id}.secondpass.txt"
    """
}


process GENERATE_STAT_CONTIG_REPORT {

    tag "${sample_id}"

    input:
    tuple val(sample_id), path(secondpass_file)

    output:
    tuple val(sample_id), path("${sample_id}.report")

    script:
    """
    hits_to_report.py \
    --cutoff-percent ${params.cutoff_percent} \
    ${secondpass_file} \
    ${sample_id}.report
    """
}

process IDENTIFY_HUMAN_VIRUS_FAMILY_CONTIGS {

    tag "${sample_id}"

    input:
    tuple val(sample_id), path(secondpass_file)

    output:
    tuple val(sample_id), path("${sample_id}.report")

    script:
    // TODO: This script is hardcoded only to work with snakemake
    """
    extract_taxa_spots.py
    """
}
