process FASTQC_REPORT {

    /*
    */

	tag "${sample_id}"
    label "general"
	publishDir params.fastqc_results, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 1

	input:
	tuple val(sample_id), path(reads)

	output:
	path "${sample_id}_qc.html", emit: html
	path "${sample_id}/", emit: multiqc_data

	script:
	"""
	fqc -q ${reads} -s . > ${sample_id}_qc.html
	mkdir ${sample_id}
	mv fastqc_data.txt ${sample_id}/fastqc_data.txt
	"""

}

process MULTIQC_REPORT {

    /*
    */
	
    label "multiqc"
	publishDir params.preprocessing, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 1

	input:
	path fastqc_files

	output:
	path("*.html")

	script:
	"""
	multiqc ${fastqc_files}
	"""

}
