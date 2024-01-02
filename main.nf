#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
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
			.fromPath( "${params.fastq_dir}/**/*.fastq*" )
			.map { fastq -> tuple( file(fastq).getSimpleName(), file(fastq) ) }
	}
    
	ch_urls = Channel
		.fromList( params.accession2taxid_urls )
		.flatten()
	
	
	// Workflow steps
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

	GI2TAXID (
		FETCH_NT.out
	)

    SORT_BY_NAME (
        GI2TAXID.out,
		CONSTRUCT_TAX_TREE.out
    )

    // SKETCH_BLACKLIST (
    //     SORT_BY_NAME.out
    // )

    SKETCH_NT_WITH_BBSKETCH (
        GI2TAXID.out
    )

    CLASSIFY_WITH_BBSKETCH (
        ch_fastqs,
        SKETCH_WITH_BBSKETCH.out.collect()
    )

	SKETCH_NT_WITH_SYLPH (
        GI2TAXID.out
	)

	SKETCH_SAMPLE_WITH_SYLPH (
		ch_fastqs
	)

	CLASSIFY_WITH_SYLPH (
		SKETCH_NT_WITH_SYLPH.out,
		SKETCH_WITH_SYLPH.out
	)

	SKETCH_NT_WITH_SOURMASH (
		GI2TAXID.out
	)

	SKETCH_SAMPLE_WITH_SOURMASH (
		ch_fastqs
	)

	SOURMASH_GATHER (
		SKETCH_WITH_SOURMASH.out,
		SKETCH_NT_WITH_SOURMASH.out
	)
	
}
// --------------------------------------------------------------- //



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //
// Additional parameters that are derived from parameters set in nextflow.config

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

process FETCH_ACCESSION2TAXID {

	storeDir params.taxpath

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	cpus 2

	input:
	val url

	output:
	path "shrunk.${file_name}"

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
    path "nt.gz"
	
	script:
	"""
	wget -q ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nt.gz
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

process SKETCH_NT_WITH_BBSKETCH {
	
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
	params.download_only == false && bbsketch
	
	script:
	"""
	bbsketch.sh \
    in=${nt_fasta} \
    out=taxa#.sketch \
    k=32,24 autosize=t depth=t minsize=300 \
    server=f prefilter=t tossjunk=t ow unpigz \
    mode=taxa taxpath=${params.taxpath} \
    tree=${params.taxpath}/tree.taxtree.gz files=31 ow unpigz \
    minsize=300 prefilter autosize k=32,24 depth
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
	comparesketch.sh -Xmx${task.memory}g \
	in=${fastq} out=${sample_id}_profiled.tsv \
	tree=${params.taxdir}/tree.taxtree.gz taxa*.sketch \
	k=32,24 mode=sequence level=1 format=3 records=1 printtaxa=t ow \
	exclude=1923094,Potexvirus,Virgaviridae,Bromoviridae,191289,Tymoviridae,Carlavirus && \
	cat ${sample_id}_profiled.tsv | awk 'NR==1 || /virus/' > 02_clean.virus_only.bbmap_profiled.tsv
	"""
}

process SKETCH_NT_WITH_SYLPH {
	
	/* */

	storeDir params.nt_storedir

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	cpus 8

	input:
	path nt_db

	output:
	tuple val(sample_id), path("nt_k31.syldb")

	when:
	params.download_only == false && sylph

	script:
	"""
	sylph sketch \
	--threads ${task.cpus} -k 31 \
	${nt_db} -o nt_k31
	"""

}

process SKETCH_SAMPLE_WITH_SYLPH {
	
	/* */
	
	tag "${sample_id}"
	publishDir params.sylph_sketches, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}.sylsp")

	when:
	params.download_only == false && sylph

	script:
	"""
	sylph sketch \
	--threads ${task.cpus} -k 31 --individual-records \
	${reads} -o ${sample_id}
	"""

}

process CLASSIFY_WITH_SYLPH {
	
	/* */
	
	tag "${sample_id}"
	publishDir params.sylph_classifications, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	input:
	tuple val(sample_id), path(sample_sketches)
	each path(nt_syldb)

	output:
	path "${sample_id}_sylph_results.tsv"

	script:
	"""
	sylph profile \
	-t 12 --minimum-ani 75 --estimate-unknown \
	${sample_sketches} ${nt_syldb} \
	| csvtk sort -t -k "13:nr" -l > ${sample_id}_sylph_results.tsv
	"""
	
}

process SKETCH_NT_WITH_SOURMASH {
	
	/* */

	storeDir params.nt_storedir

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	input:
	path nt_db

	output:
	path "nt_k31.sig.gz"

	when:
	params.download_only == false && sourmash

	script:
	"""
	sourmash sketch dna -p k=31,scaled=1000,abund -f -o nt_k31.sig.gz ${nt_db}
	sourmash sketch dna -p scaled=1000,k=31 ${reads} -o ${sample_id}_reads.sig && \
	sourmash index ${sample_id}_reads ${sample_id}_reads.sig
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
	params.download_only == false && sourmash

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
	params.download_only == false && sourmash

	script:
	"""
	sourmash gather -p abund ${nt_sigs} ${sample_sigs} -o ${sample_id}_sourmash_results.csv
	"""

}

// --------------------------------------------------------------- //
