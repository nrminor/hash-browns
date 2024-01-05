#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //

// prints to the screen and to the log
log.info	"""
			HASH-BROWNS (version 0.1.0)
			===================================
			fastq_dir       : ${params.fastq_dir}
			results_dir     : ${params.results}
			query_fasta     : ${params.query_fasta}

			Storage directories:
			-----------------------------------
			Taxonomy dir    : ${params.taxpath}
			NCBI NT dir     : ${params.nt_storedir}

			Chosen tools:
			-----------------------------------
			BBSketch        : ${params.bbsketch}
			Sylph           : ${params.sylph}
			Sourmash        : ${params.sourmash}
			Centrifuge      : Coming soon!

			Run settings:
			-----------------------------------
			[fast_mode      : ${params.fast_mode}]
			[realtime_dir   : ${params.realtime_dir}]
			[cleanup        : ${params.cleanup}]
			[download only? : ${params.download_only}]
			[available cpus : ${params.available_cpus}]
			[run date       : ${params.date}]
			"""
			.stripIndent()

workflow {

	assert params.taxpath : "Please set path on your system for storing taxonomy files."
	assert params.nt_storedir : "Please set path on your system for storing the NCBI nt database."
	
	// input channels
	if ( params.realtime_dir ) {
		ch_fastqs = Channel
			.watchPath( "${params.realtime_dir}/**/*.fastq*", 'create,modify' )
			.map { fastq -> tuple( file(fastq).getSimpleName(), file(fastq) ) }
	} else {
		ch_fastqs = Channel
			.fromPath( "${params.fastq_dir}/*.fastq*" )
			.map { fastq -> tuple( file(fastq).getSimpleName(), file(fastq) ) }
	}
    
	ch_urls = Channel
		.fromList( params.accession2taxid_urls )
		.flatten()

		ch_data_manifest = Channel
			.fromPath( params.data_manifest )
	
	
	// Workflow steps
	VALIDATE_SEQS (
		ch_fastqs
	)
	
	READ_QC (
		VALIDATE_SEQS.out
	)

	FASTQC_REPORT (
		READ_QC.out
	)

	MULTIQC_REPORT (
		FASTQC_REPORT.out.multiqc_data.collect()
	)

	FETCH_FAST_MODE_DB ( )

	FETCH_ACCESSION2TAXID (
		ch_urls
	)

	FETCH_TAXONOMY ( )

	UNZIP_TAXONOMY (
		FETCH_TAXONOMY.out
	)

	CONSTRUCT_TAX_TREE (
		UNZIP_TAXONOMY.out
	)

	CONSTRUCT_GITABLE (
		FETCH_ACCESSION2TAXID.out.collect()
	)

	ANALYZE_ACCESSIONS (
		FETCH_ACCESSION2TAXID.out.collect()
	)

    FETCH_NT ( )

	// GI2TAXID (
	// 	FETCH_NT.out
	// )

    // SORT_BY_NAME (
    //     GI2TAXID.out,
	// 	CONSTRUCT_TAX_TREE.out
    // )

    // SKETCH_BLACKLIST (
    //     SORT_BY_NAME.out
    // )

    SKETCH_DB_WITH_BBSKETCH (
        FETCH_NT.out
			.mix(
				FETCH_FAST_MODE_DB.out
			)
    )

    CLASSIFY_WITH_BBSKETCH (
        READ_QC.out,
        SKETCH_DB_WITH_BBSKETCH.out.collect()
    )

	SKETCH_DB_WITH_SYLPH (
        FETCH_NT.out
			.mix(
				FETCH_FAST_MODE_DB.out
			)
	)

	SKETCH_SAMPLE_WITH_SYLPH (
		READ_QC.out
	)

	CLASSIFY_WITH_SYLPH (
		SKETCH_DB_WITH_SYLPH.out,
		SKETCH_SAMPLE_WITH_SYLPH.out
	)

	SKETCH_DB_WITH_SOURMASH (
		FETCH_NT.out
			.mix(
				FETCH_FAST_MODE_DB.out
			)
	)

	SKETCH_SAMPLE_WITH_SOURMASH (
		READ_QC.out
	)

	SOURMASH_GATHER (
		SKETCH_SAMPLE_WITH_SOURMASH.out,
		SKETCH_DB_WITH_SOURMASH.out
	)
	
}
// --------------------------------------------------------------- //



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //
// Additional parameters that are derived from parameters set in nextflow.config
params.preprocessing = params.results + "/preprocessing"
params.read_checks = params.preprocessing + "/1_read_checks"
params.filtered = params.preprocessing + "/2_filtered_reads"
params.fastqc_results = params.preprocessing + "/3_FastQC_reports"

// bbsketch results
params.bbsketch_results = params.results + "/bbsketch"
params.sorted_nt = params.bbsketch_results + "/sorted_nt"
params.bbsketches = params.bbsketch_results + "/sketches"
params.bbsketch_classifications = params.bbsketch_results + "/classifications"

// sylph results
params.sylph_results = params.results + "/sylph"
params.sylph_sketches = params.sylph_results + "/sketches"
params.sylph_classifications = params.sylph_results + "/classifications"

// sourmash results
params.sourmash_results = params.results + "/sourmash"
params.sourmash_sketches = params.sourmash_results + "/sketches"
params.sourmash_classifications = params.sourmash_results + "/classifications"

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

// process FILTLONG {}

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
	path("*.html")

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
	path "human_virus_db.fa.gz"

	when:
	fast_mode == true

	script:
	"""
	sdf pull --urls --overwrite && \
	sdf status > ${params.date}_status_check.txt
	"""

}

process FETCH_ACCESSION2TAXID {

	storeDir params.taxpath

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	cpus 2

	input:
	val url

	output:
	path "shrunk.${file_name}"

	when:
	params.fast_mode == false

	script:
	file_name = url.toString().split("/")[-1]
	"""
	wget -q -O - ${url} \
	| shrinkaccession.sh in=stdin.txt.gz out=shrunk.${file_name} zl=9 t=${task.cpus}
	"""

}

process FETCH_TAXONOMY {

	storeDir params.taxpath

	output:
	path "taxdmp.zip"

	script:
	"""
	wget -nv ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip
	"""
}

process UNZIP_TAXONOMY {

	storeDir params.taxpath

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	input:
	path taxdmp_zip

	output:
	path "*.dmp"

	script:
	"""
	unzip -o ${taxdmp_zip}
	"""

}

process CONSTRUCT_TAX_TREE {

	storeDir params.taxpath

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	memory 16.GB

	input:
	path taxdmp_files

	output:
	path "tree.taxtree.gz"

	script:
	"""
	taxtree.sh names.dmp nodes.dmp merged.dmp tree.taxtree.gz -Xmx16g
	"""

}

process CONSTRUCT_GITABLE {

	storeDir params.taxpath

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	input:
	path shrunk_taxid2accessions

	output:
	path "gitable.int1d.gz"
	
	script:
	"""
	gitable.sh \
	shrunk.dead_nucl.accession2taxid.gz,shrunk.dead_prot.accession2taxid.gz,shrunk.dead_wgs.accession2taxid.gz,shrunk.nucl_gb.accession2taxid.gz,shrunk.nucl_wgs.accession2taxid.gz,shrunk.pdb.accession2taxid.gz,shrunk.prot.accession2taxid.gz \
	gitable.int1d.gz
	"""

}

process ANALYZE_ACCESSIONS {

	storeDir params.taxpath

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	input:
	path shrunk_taxid2accessions

	output:
	path "patterns.txt"

	script:
	"""
	analyzeaccession.sh shrunk.*.accession2taxid.gz out=patterns.txt
	"""

}

process FETCH_NT {
	
	/* */

	storeDir params.nt_storedir

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1
	
	output:
    path "nt.fa.gz"

	when:
	params.fast_mode == false
	
	script:
	"""
	wget -q -O - ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nt.gz \
	| gi2taxid.sh -Xmx1g \
    in=stdin.fa.gz out=nt.fa.gz \
    pigz=32 unpigz=t bgzip=t preferbgzip=t zl=8 server=f ow shrinknames maxbadheaders=5000 \
    badheaders=badHeaders.txt taxpath=${params.taxpath}
	"""
}

process GI2TAXID {
	
	/* */

	storeDir params.nt_storedir

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	input:
	path nt
	
	output:
    path "nt.fa.gz"
	
	script:
	"""
	gi2taxid.sh -Xmx1g \
    in=`realpath ${nt}` out=nt.fa.gz \
    pigz=32 unpigz=t bgzip=t preferbgzip=t zl=8 server=f ow shrinknames maxbadheaders=5000 \
    badheaders=badHeaders.txt taxpath=${params.taxpath}
	"""
}

// process SORT_BY_NAME {
	
// 	/* */

// 	storeDir params.nt_storedir

// 	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
// 	maxRetries 1

// 	memory 64.GB
	
// 	input:
// 	path nt_fasta
// 	path taxtree
	
// 	output:
// 	path "nt_sorted.fa.gz"
	
// 	script:
// 	"""
// 	sortbyname.sh -Xmx64g \
//     in=`realpath ${nt_fasta}` out=nt_sorted.fa.gz \
//     ow taxa taxpath=${params.taxpath} tree="tree.taxtree.gz" fastawrap=1023 zl=9 fixjunk \
// 	pigz=32 minlen=60 bgzip unbgzip
// 	"""
// }

// process SKETCH_BLACKLIST {
	
// 	/* */
	
// 	tag "${tag}"
// 	publishDir params.results, mode: 'copy'
	
// 	input:
// 	path nt_sorted
	
// 	output:
// 	path "blacklist_nt_genus_100.sketch"
	
// 	script:
// 	"""
// 	sketchblacklist.sh -Xmx31g \
//     in=`realpath ${nt_sorted}` out=blacklist_nt_genus_100.sketch \
//     prepasses=1 tree="tree.taxtree.gz" taxa taxlevel=genus ow mincount=120 k=32,24 depth taxpath=${params.taxpath}
// 	"""
// }

process SKETCH_DB_WITH_BBSKETCH {
	
	/* */

	storeDir params.nt_storedir

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	memory 32.GB
	
	input:
	path nt_fasta
	
	output:
	path "taxa*.sketch"

	when:
	params.download_only == false && params.bbsketch == true
	
	script:
	"""
	bbsketch.sh -Xmx32g \
    in=`realpath ${nt_fasta}` \
    out=taxa.sketch \
    k=32,24 autosize=t depth=t minsize=300 \
    prefilter=t tossjunk=t ow unpigz files=1 \
    mode=taxa taxpath=${params.taxpath} \
    tree=${params.taxpath}/tree.taxtree.gz
	"""
}

process CLASSIFY_WITH_BBSKETCH {
	
	/* */
	
	tag "${sample_id}"
	publishDir params.bbsketch_classifications, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	cpus 3
	memory 32.GB
	
	input:
    tuple val(sample_id), path(fastq)
	each path(nt_sketches)
	
	output:
	path "*"
	
	script:
	"""
	seqkit seq -v ${fastq} \
	| comparesketch.sh -Xmx32g \
	in=stdin.fastq out=${sample_id}_profiled.tsv \
	tree=${params.taxpath}/tree.taxtree.gz taxa.sketch \
	k=32,24 mode=sequence level=1 format=3 records=1 ow sortbyani=t \
	printtaxa=t printdepth=t sortbydepth=t printunique=t printunique2=t && \
	cat ${sample_id}_profiled.tsv | awk 'NR==1 || /virus/' > ${sample_id}.virus_only.bbmap_profiled.tsv
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
	params.download_only == false && params.sylph == true

	script:
	"""
	sylph sketch -t ${task.cpus} -k 31 -i -c 200 -g ${nt_db} -o nt_c200_k31
	"""

}

process SKETCH_SAMPLE_WITH_SYLPH {
	
	/* */
	
	tag "${sample_id}"
	// publishDir params.sylph_sketches, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	cpus 3

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

	cpus params.available_cpus

	input:
	each path(nt_syldb)
	tuple val(sample_id), path(sample_sketches)

	output:
	path "*"

	script:
	"""
	sylph profile \
	-t ${task.cpus} --minimum-ani 90 --estimate-unknown -M 3 --read-seq-id 0.80 \
	${sample_sketches} ${nt_syldb} > ${sample_id}_sylph_results.tsv && \
	csvtk sort -t -k "5:nr" -l ${sample_id}_sylph_results.tsv \
	| csvtk grep -t -f "Contig_name" -r -p "virus" -o ${sample_id}_sylph_virus_only_results.tsv
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
	params.download_only == false && params.sourmash == true

	script:
	"""
	sourmash sketch dna -p k=31,scaled=1000,abund -f -o nt_k31.sig.gz ${nt_db}
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
	tuple val(sample_id), path("${sample_id}_reads.sig"), path("${sample_id}_reads.sbt.zip")

	when:
	params.download_only == false && params.sourmash == true

	script:
	"""
	sourmash sketch dna -p scaled=1000,k=31 ${reads} -o ${sample_id}_reads.sig && \
	sourmash index ${sample_id}_reads ${sample_id}_reads.sig
	"""

}

process SOURMASH_GATHER {
	
	/* */
	
	tag "${sample_id}"
	publishDir params.sourmash_sketches, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	input:
	tuple val(sample_id), path(sample_sigs), path(index)
	each path(nt_sigs)

	output:

	when:
	params.download_only == false && params.sourmash == true

	script:
	"""
	sourmash gather -p abund ${nt_sigs} ${sample_sigs} -o ${sample_id}_sourmash_results.csv
	"""

}

// --------------------------------------------------------------- //
