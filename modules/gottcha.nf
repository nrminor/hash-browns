process GOTTCHA2_PROFILE_NANOPORE {

    tag "${sample_id}"
    publishDir params.gottcha_sam, mode: 'copy', overwrite: false, pattern: "*.sam"
    publishDir params.gottcha_stats, mode: 'copy', overwrite: false, pattern: "*.tsv"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }

    cpus params.gottcha2_cpus

    input:
    tuple path(ref_mmi), path(stats), path(tsv), val(sample_id), path(fastq)

    output:
    tuple path(ref_mmi), path(stats), path(tsv), val(sample_id), path("${sample_id}*.sam"), emit: aligned
    tuple path(ref_mmi), path(stats), path(tsv), val(sample_id), path("${sample_id}*.full.tsv"), emit: full_tsv
    path "*.tsv", emit: all_stats

    script:
    def ref_prefix = file(ref_mmi).getBaseName().toString().replace(".mmi", "")
    """
    gottcha2.py  \
    --database ${ref_prefix} \
    --prefix ${sample_id} \
     -t ${task.cpus} \
    -i ${fastq} \
     --nanopore
    """
}

process GOTTCHA2_PROFILE_ILLUMINA {

    tag "${sample_id}"
    publishDir params.gottcha_sam, mode: 'copy', overwrite: false, pattern: "*.sam"
    publishDir params.gottcha_stats, mode: 'copy', overwrite: false, pattern: "*.tsv"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }

    cpus params.gottcha2_cpus

    input:
    tuple path(ref_mmi), path(stats), path(tsv), val(sample_id), path(fastq1), path(fastq2)

    output:
    tuple path(ref_mmi), path(stats), path(tsv), val(sample_id), path("${sample_id}*.sam"), emit: aligned
    tuple path(ref_mmi), path(stats), path(tsv), val(sample_id), path("${sample_id}*.full.tsv"), emit: full_tsv
    path "*.tsv", emit: all_stats

    script:
    def ref_prefix = file(ref_mmi).getBaseName().toString().replace(".mmi", "")
    """
    gottcha2.py -d ${ref_prefix} -t ${task.cpus} -i ${fastq1} ${fastq2}
    """
}

process GENERATE_FASTA {

    tag "${sample_id}"
    publishDir params.gottcha_fasta, mode: 'copy', overwrite: false, pattern: "*.fasta"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }

    cpus params.gottcha2_cpus

    input:
    tuple path(ref_mmi), path(stats), path(tsv), val(sample_id), path("*extract*"), emit: extracted_reads
    path "*.*", emit: all_files

    output:
    path "*.fasta"

    script:
    def ref_prefix = file(ref_mmi).getBaseName().toString().replace(".mmi", "")
    """
    gottcha2.py \
    -s ${sample_id}.gottcha_species.sam \
    -ef \
    --nanopore \
    --database ${ref_prefix}
    """
}
