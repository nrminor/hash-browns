process SYLPH_SKETCH_DB {

	/* */

	storeDir params.db_cache

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	cpus params.available_cpus

	input:
	path input_db

	output:
	path "*.syldb"

	when:
	(params.tools && params.tools.contains("sylph")) || params.all || params.sylph

	script:
	"""
	sylph sketch -t ${task.cpus} -k 31 -i -c 200 -g ${input_db} -o ${input_db}.syldb
	"""
}

process SYLPH_SKETCH_SAMPLE {

	/*
	-c is the subsampling rate, also referred to as the compression parameter.
	A higher -c is faster but less sensitive at low coverage 
	*/

	tag "${sample_id}"

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	cpus 4

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}*.sylsp")

	when:
	(params.tools && params.tools.contains("sylph")) || params.all || params.sylph

	script:
	"""
	sylph sketch -t ${task.cpus} -k ${params.k} -c 100 -r ${reads} -o ${sample_id}
	"""
}

process CLASSIFY_WITH_SYLPH {

	/* */

	tag "${sample_id}"

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	input:
	tuple val(sample_id), path(sample_sketches), path(syldb)

	output:
	tuple val(sample_id), path("${sample_id}*.tsv")

	when:
	(params.tools && params.tools.contains("sylph")) || params.all || params.sylph

	script:
	"""
	sylph profile \
	-t ${task.cpus} --minimum-ani 90 --estimate-unknown -M 3 --read-seq-id 0.80 \
	${sample_sketches} ${syldb} > ${sample_id}_sylph_results.tsv
	"""
}

process SYLPH_TAX_DOWNLOAD {

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	output:
	val "ready"

	when:
	(params.tools && params.tools.contains("sylph")) || params.all || params.sylph

	script:
	"""
	if [ ! -d ${params.sylph_tax_dir} ]; then
		mkdir -p ${params.sylph_tax_dir}
	fi
	sylph-tax download --download-to ${params.sylph_tax_dir}
	"""
}

process SYLPH_TAXPROF {

	/* */

	tag "${sample_id}"

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	input:
	tuple val(sample_id), path(tsv_dir), val(ready)

	script:
	"""
	sylph-tax taxprof sylph_results/*.tsv -t ${params.sylph_taxonomy} -o ${sample_id}
	"""
}

process EXTRACT_HUMAN_VIRUSES {

	/* */

	tag "${sample_id}"

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	input:
	tuple val(sample_id), path(tsv)

	output:
	tuple val(sample_id), "*.tsv"

	when:
	(params.tools && params.tools.contains("sylph")) || params.all || params.sylph

	script:
	"""
	csvtk sort -t -k "5:nr" -l ${tsv} \
	| csvtk grep -t --ignore-case -f "Contig_name" -r -p virus \
	| csvtk grep -t --ignore-case -f "Contig_name" -r -p human \
	-o ${sample_id}_sylph_human_virus_only_results.tsv
	"""
}
