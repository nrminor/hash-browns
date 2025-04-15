process GOTTCHA2_PROFILE_NANOPORE {

    tag "${sample_id}"
    // publishDir params.gottcha_sam, mode: 'copy', overwrite: false, pattern: "*.sam"
    // publishDir params.gottcha_stats, mode: 'copy', overwrite: false, pattern: "*.tsv"

    maxForks { params.max_tasks ? params.max_tasks : params.available_cpus }
    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }

    cpus 12

    input:
    tuple val(sample_id), path(fastq), path(ref_db)

    output:
    tuple val(sample_id), path("${sample_id}*.sam"), path(ref_db), emit: aligned
    tuple val(sample_id), path("${sample_id}*.full.tsv"), path(ref_db), emit: full_tsv
    path "*.tsv", emit: all_stats

    when:
    (params.tools && params.tools.contains("gottcha2") || params.tools.contains("gottcha")) || params.all || params.gottcha2

    script:
    def ref_prefix = file(ref_db[0]).getBaseName().toString().replace(".mmi", "")
    """
    gottcha2.py  \
    --database ${ref_prefix} \
    --prefix ${sample_id} \
    --noCutoff ---dbLevel strain --threads ${task.cpus} \
    --nanopore \
    --input ${fastq}
    """
}

process GOTTCHA2_PROFILE_ILLUMINA {

    tag "${sample_id}"
    // publishDir params.gottcha_sam, mode: 'copy', overwrite: false, pattern: "*.sam"
    // publishDir params.gottcha_stats, mode: 'copy', overwrite: false, pattern: "*.tsv"

    maxForks { params.max_tasks ? params.max_tasks : params.available_cpus }
    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }

    cpus 12

    input:
    tuple val(sample_id), path(fastq), path(ref_db)

    output:
    tuple val(sample_id), path("${sample_id}*.sam"), path(ref_db), emit: aligned
    tuple val(sample_id), path("${sample_id}*.full.tsv"), path(ref_db), emit: full_tsv
    path "*.tsv", emit: all_stats

    when:
    (params.tools && params.tools.contains("gottcha2") || params.tools.contains("gottcha")) || params.all || params.gottcha2

    script:
    def ref_prefix = file(ref_db[0]).getBaseName().toString().replace(".mmi", "")
    """
    gottcha2.py \
    --database ${ref_prefix} \
    --prefix ${sample_id} \
    --noCutoff ---dbLevel strain --threads ${task.cpus} \
    --input ${fastq}
    """
}

process GENERATE_FASTA {

    tag "${sample_id}"
    // publishDir params.gottcha_fasta, mode: 'copy', overwrite: false, pattern: "*.fasta"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }

    input:
    tuple val(sample_id), path(sam), path(ref_db)

    output:
    path "*"

    when:
    (params.tools && params.tools.contains("gottcha2") || params.tools.contains("gottcha")) || params.all || params.gottcha

    script:
    def ref_prefix = file(ref_db[0]).getBaseName().toString().replace(".mmi", "")
    """
    gottcha2.py \
    --noCutoff ---dbLevel strain --threads ${task.cpus} -ef \
    --database ${ref_prefix} \
    --prefix ${sample_id} \
    --sam ${sam}
    """
}
