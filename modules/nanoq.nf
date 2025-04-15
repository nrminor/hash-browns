process READ_QC {

	/*
	*/

	tag "${sample_id}"
	label "general"
	// publishDir params.filtered, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 1

	input:
	tuple val(sample_id), val(platform), path(reads), path(report)

	output:
	tuple val(sample_id), val(platform), path("${sample_id}_nanoq.fastq.gz")

	script:
	"""
	nanoq -i `realpath ${reads}` \
	--min-len 200 --min-qual 10 \
	-r ${sample_id}_nanoq_report.txt \
	> ${sample_id}_nanoq.fastq.gz
	"""
}
