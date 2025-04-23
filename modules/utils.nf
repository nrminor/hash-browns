process PREPARE_NVD_INPUT {
    input:
    tuple val(sample_id), path(fastq)

    output:
    tuple val(sample_id), val(sample_type), path(fastq)

    script:
    def sample_type = sample_id.contains("SRR") ? 'ont' : 'illumina'
}
