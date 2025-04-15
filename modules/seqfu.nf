process VALIDATE_SEQS {

	/*
    */

	tag "${sample_id}"
	label "general"
	publishDir params.read_checks, pattern: "*.tsv", mode: 'copy', overwrite: true

	errorStrategy 'ignore'

	cpus 1

	input:
	tuple val(sample_id), val(platform), path(reads)

	output:
	tuple val(sample_id), val(platform), path(reads), path("${sample_id}_seqfu_report.tsv")

	script:
	"""
	seqfu check \
	--deep --thousands \
	`realpath ${reads}` > ${sample_id}_seqfu_report.tsv
	"""
}
