#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //

workflow {

	// assert params.taxpath : "Please set path on your system for storing taxonomy files."
	// assert params.nt_storedir : "Please set path on your system for storing the NCBI nt database."

	// prints to the screen and to the log
	// https://patorjk.com/software/taag/#p=testall&t=HASH-BROWNS
	log.info(
		"""
			=======================================================================================================================
			=  ====  =====  ======      ===  ====  ============      ====       =====    ====  ====  ====  ==  =======  ===      ==
			=  ====  ====    ====  ====  ==  ====  ============  ===  ===  ====  ===  ==  ===  ====  ====  ==   ======  ==  ====  =
			=  ====  ===  ==  ===  ====  ==  ====  ============  ====  ==  ====  ==  ====  ==  ====  ====  ==    =====  ==  ====  =
			=  ====  ==  ====  ===  =======  ====  ============  ===  ===  ===   ==  ====  ==  ====  ====  ==  ==  ===  ===  ======
			=        ==  ====  =====  =====        ==        ==      ====      ====  ====  ==   ==    ==  ===  ===  ==  =====  ====
			=  ====  ==        =======  ===  ====  ============  ===  ===  ====  ==  ====  ===  ==    ==  ===  ====  =  =======  ==
			=  ====  ==  ====  ==  ====  ==  ====  ============  ====  ==  ====  ==  ====  ===  ==    ==  ===  =====    ==  ====  =
			=  ====  ==  ====  ==  ====  ==  ====  ============  ===  ===  ====  ===  ==  =====    ==    ====  ======   ==  ====  =
			=  ====  ==  ====  ===      ===  ====  ============      ====  ====  ====    =======  ====  =====  =======  ===      ==
			=======================================================================================================================

			HASH-BROWNS: Metagenomic read classification well-done.
			(version 0.1.0)
			===================================
			fastq directory : ${params.fastq_dir}
			results dir     : ${params.results}

			Storage directories:
			-----------------------------------
			

			Chosen tools:
			-----------------------------------
			Sylph           : ${params.sylph}
			Sourmash        : ${params.sourmash}
			GOTTCHA2        : Coming soon!
			Strobealign     : Coming soon!
			KMCP            : Coming soon!
			Mmseqs2         : Coming soon!
			CLARK           : Coming soon!
			MetaPHlAn4      : Coming soon!
			Centrifuge      : Coming soon!
			Megan           : Coming soon!
			mOTUs           : Coming soon!

			Run settings:
			-----------------------------------
			[realtime_dir   : ${params.realtime_dir}]
			[cleanup        : ${params.cleanup}]
			[download only? : ${params.download_only}]
			[available cpus : ${params.available_cpus}]
			[run date       : ${params.date}]

			""".stripIndent()
	)

	// input channels
	if (params.realtime_dir) {
		ch_fastqs = Channel.watchPath("${params.realtime_dir}/**/*.fastq*", 'create,modify')
			.map { fastq -> tuple(file(fastq), file(fastq).countFastq()) }
			.filter { it[1] >= 20 }
			.map { fastq -> tuple(file(fastq).getSimpleName(), file(fastq)) }
	}
	else {
		ch_fastqs = Channel.fromPath("${params.fastq_dir}/*.fastq*")
			.map { fastq -> tuple(file(fastq), file(fastq).countFastq()) }
			.filter { it[1] >= 20 }
			.map { fastq, count -> tuple(file(fastq).getSimpleName(), file(fastq)) }
	}

	ch_urls = Channel.fromList(params.accession2taxid_urls)
		.flatten()

	ch_data_manifest = Channel.fromPath(params.data_manifest)

	// Workflow steps
	VALIDATE_SEQS(
		ch_fastqs
	)

	READ_QC(
		VALIDATE_SEQS.out
	)

	FASTQC_REPORT(
		READ_QC.out
	)

	MULTIQC_REPORT(
		FASTQC_REPORT.out.multiqc_data.collect()
	)
}
// --------------------------------------------------------------- //



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //
// Additional parameters that are derived from parameters set in nextflow.config

// preprocessing results
params.preprocessing            = params.results + "/preprocessing"
params.read_checks              = params.preprocessing + "/1_read_checks"
params.filtered                 = params.preprocessing + "/2_filtered_reads"
params.fastqc_results           = params.preprocessing + "/3_FastQC_reports"

// bbsketch results
params.bbsketch_results         = params.results + "/bbsketch"
params.sorted_nt                = params.bbsketch_results + "/sorted_nt"
params.bbsketches               = params.bbsketch_results + "/sketches"
params.bbsketch_classifications = params.bbsketch_results + "/classifications"

// sylph results
params.sylph_results            = params.results + "/sylph"
params.sylph_sketches           = params.sylph_results + "/sketches"
params.sylph_classifications    = params.sylph_results + "/classifications"

// sourmash results
params.sourmash_results         = params.results + "/sourmash"
params.sourmash_sketches        = params.sourmash_results + "/sketches"
params.sourmash_classifications = params.sourmash_results + "/classifications"

// CPUs to use when sharing
params.shared_cpus              = Math.floor(params.available_cpus / 2)

// --------------------------------------------------------------- //




// PROCESS SPECIFICATION 
// --------------------------------------------------------------- //

process VALIDATE_SEQS {

	/*
    */

	tag "${sample_id}"
	label "general"
	publishDir params.read_checks, pattern: "*.tsv", mode: 'copy', overwrite: true

	errorStrategy 'ignore'

	cpus 1

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path(reads), path("${sample_id}_seqfu_report.tsv")

	script:
	"""
	seqfu check \
	--deep --thousands \
	`realpath ${reads}` > ${sample_id}_seqfu_report.tsv
	"""
}

process READ_QC {

	/*
	*/

	tag "${sample_id}"
	label "general"
	publishDir params.filtered, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 1

	input:
	tuple val(sample_id), path(reads), path(report)

	output:
	tuple val(sample_id), path("${sample_id}_nanoq.fastq.gz")

	script:
	"""
	nanoq -i `realpath ${reads}` \
	--min-len 200 --min-qual 10 \
	-r ${sample_id}_nanoq_report.txt \
	> ${sample_id}_nanoq.fastq.gz
	"""
}

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
	path "*.html"

	script:
	"""
	multiqc ${fastqc_files}
	"""
}

process FETCH_FAST_MODE_DB {

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	input:
	path data_manifest

	output:
	path "resources/human_virus_db.fa.gz"

	when:
	params.fast_mode == true

	script:
	"""
	sdf pull --urls --overwrite && \
	sdf status > ${params.date}_status_check.txt
	"""
}

process SKETCH_DB_WITH_SYLPH {

	/* */

	storeDir params.nt_storedir

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	cpus params.available_cpus

	input:
	path nt_db

	output:
	path "nt_c200_k31.syldb"

	when:
	params.download_only == false && params.fast_mode == false && params.sylph == true

	script:
	"""
	sylph sketch -t ${task.cpus} -k 31 -i -c 200 -g ${nt_db} -o nt_c200_k31
	"""
}

process SKETCH_FAST_DB_WITH_SYLPH {

	/* */

	storeDir params.nt_storedir

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	cpus params.available_cpus

	input:
	path fast_db

	output:
	path "human_virus_c200_k31.syldb"

	when:
	params.download_only == false && params.fast_mode == true && params.sylph == true

	script:
	"""
	sylph sketch -t ${task.cpus} -k 31 -i -c 200 -g ${fast_db} -o human_virus_c200_k31
	"""
}

process SKETCH_SAMPLE_WITH_SYLPH {

	/* */

	tag "${sample_id}"
	// publishDir params.sylph_sketches, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	cpus 4

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}*.sylsp")

	when:
	params.download_only == false && params.sylph == true

	script:
	"""
	sylph sketch -t ${task.cpus} -k 31 -c 100 -r ${reads} -o ${sample_id}
	"""
}

process CLASSIFY_WITH_SYLPH {

	/* */

	tag "${sample_id}"
	publishDir params.sylph_classifications, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	cpus params.shared_cpus

	input:
	each path(nt_syldb)
	tuple val(sample_id), path(sample_sketches)

	output:
	path "${sample_id}*.tsv"

	script:
	"""
	sylph profile \
	-t ${task.cpus} --minimum-ani 90 --estimate-unknown -M 3 --read-seq-id 0.80 \
	${sample_sketches} ${nt_syldb} > ${sample_id}_sylph_results.tsv

	csvtk sort -t -k "5:nr" -l ${sample_id}_sylph_results.tsv \
	| csvtk grep -t --ignore-case -f "Contig_name" -r -p virus \
	| csvtk grep -t --ignore-case -f "Contig_name" -r -p human \
	-o ${sample_id}_sylph_human_virus_only_results.tsv
	"""
}

process SKETCH_DB_WITH_SOURMASH {

	/* */

	storeDir params.nt_storedir

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries

	memory 100.GB

	input:
	path nt_db

	output:
	path "nt_k31.sig.gz"

	when:
	params.download_only == false && params.fast_mode == false && params.sourmash == true

	script:
	"""
	sourmash sketch dna -p k=31,k=51,scaled=1000,abund --singleton -f -o nt_k31.sig.gz ${nt_db}
	"""
}

process SKETCH_FAST_DB_WITH_SOURMASH {

	/* */

	storeDir params.nt_storedir

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries

	memory 100.GB

	input:
	path fast_db

	output:
	path "human_virus_k31.sig.gz"

	when:
	params.download_only == false && params.fast_mode == true && params.sourmash == true

	script:
	"""
	sourmash sketch dna \
	-p k=31,k=51,scaled=1000,abund --singleton -f \
	-o human_virus_k31.sig.gz ${fast_db}
	"""
}

process SKETCH_SAMPLE_WITH_SOURMASH {

	/* */

	tag "${sample_id}"
	publishDir params.sourmash_sketches, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}_reads.sig")

	when:
	params.download_only == false && params.sourmash == true

	script:
	"""
	sourmash sketch dna -p scaled=1000,k=31,k=51 ${reads} -o ${sample_id}_reads.sig
	"""
}

process SOURMASH_GATHER {

	/* */

	tag "${sample_id}"
	publishDir params.sourmash_classifications, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	cpus params.available_cpus

	input:
	tuple val(sample_id), path(sample_sigs)
	each path(nt_sigs)

	output:
	path "*"

	when:
	params.download_only == false && params.sourmash == true

	script:
	"""
	sourmash gather \
	--prefetch --estimate-ani-ci --create-empty-results -k 31 \
	${sample_sigs} ${nt_sigs} -o ${sample_id}_sourmash_results.csv
	"""
}
