process GOTTCHA2_PROFILE_NANOPORE {

    tag "${sample_id}"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }

    cpus 12
    maxForks 2

    input:
    tuple val(sample_id), path(fastq), path(ref_fna), path(ref_mmi), path(ref_stats), path(ref_tsv)

    output:
    tuple val(sample_id), path("${sample_id}*.sam"), path(ref_fna), path(ref_mmi), path(ref_stats), path(ref_tsv), emit: aligned
    tuple val(sample_id), path("${sample_id}*.full.tsv"), path(ref_fna), path(ref_mmi), path(ref_stats), path(ref_tsv), emit: full_tsv
    path "*.tsv", emit: all_stats

    when:
    params.all || params.gottcha2 || (params.tools && params.tools.contains("gottcha2") || params.tools.contains("gottcha"))

    script:
    def ref_prefix = file(ref_mmi).getBaseName().toString().replace(".mmi", "")
    """
    gottcha2.py  \
    --database ./${ref_prefix} \
    --noCutoff --dbLevel strain --threads ${task.cpus} \
    --nanopore \
    --input ${fastq}
    """
}

process GOTTCHA2_PROFILE_ILLUMINA {

    tag "${sample_id}"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }

    cpus 12
    maxForks 2

    input:
    tuple val(sample_id), path(fastq), path(ref_fna), path(ref_mmi), path(ref_stats), path(ref_tsv)

    output:
    tuple val(sample_id), path("${sample_id}*.sam"), path(ref_fna), path(ref_mmi), path(ref_stats), path(ref_tsv), emit: aligned
    tuple val(sample_id), path("${sample_id}*.full.tsv"), path(ref_fna), path(ref_mmi), path(ref_stats), path(ref_tsv), emit: full_tsv
    path "*.tsv", emit: all_stats

    when:
    params.all || params.gottcha2 || (params.tools && params.tools.contains("gottcha2") || params.tools.contains("gottcha"))

    script:
    def ref_prefix = file(ref_mmi).getBaseName().toString().replace(".mmi", "")
    """
    gottcha2.py \
    --database ./${ref_prefix} \
    --noCutoff --dbLevel strain --threads ${task.cpus} \
    --input ${fastq}
    """
}

process GENERATE_FASTA {

    tag "${sample_id}"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }

    input:
    tuple val(sample_id), path(sam), path(ref_fna), path(ref_mmi), path(ref_stats), path(ref_tsv)

    output:
    path "*"

    when:
    params.all || params.gottcha2 || (params.tools && params.tools.contains("gottcha2") || params.tools.contains("gottcha"))

    script:
    def ref_prefix = file(ref_mmi).getBaseName().toString().replace(".mmi", "")
    """
    gottcha2.py \
    --noCutoff --dbLevel strain --threads ${task.cpus} -ef \
    --database ./${ref_prefix} \
    --sam ${sam}
    """
}
