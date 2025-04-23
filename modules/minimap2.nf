process MAP_READS_TO_CONTIGS {

    tag "${sample_id}"

    cpus 4

    input:
    tuple val(sample_id), path(contigs), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai")

    script:
    def preset = sample_id.contains('ont') || sample_id.contains('sra')
        ? 'map-ont'
        : 'sr'
    """
    minimap2 -ax ${preset} -t ${task.cpus} ${contigs} ${reads} | \
    samtools view -b -F 4 | \
    samtools sort -@ ${task.cpus} -o ${sample_id}.bam && \
    samtools index ${sample_id.bam}
    """
}

process COUNT_MAPPED_READS {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(bam), path(bam_index)

    output:
    tuple val(sample_id), path("${sample_id}.filtered.bam"), path("${sample_id}.filtered.bam.bai")
    path "${sample_id}_mapped_counts.txt"

    script:
    """
    samtools view -F 2304 -b ${bam} > ${sample_id}.filtered.bam
    samtools index ${sample_id}.filtered.bam
    samtools idxstats ${sample_id}.filtered.bam | awk '{{print $1, $3}}' > "${sample_id}_mapped_counts.txt"
    """
}

