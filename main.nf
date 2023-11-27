#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {

	assert params.taxpath : "Please set path on your system for storing taxonomy files."
	assert params.nt_storedir : "Please set path on your system for storing the NCBI nt database."
	
	// input channels
    ch_fastqs = Channel
        .fromPath( "${params.fastq_dir}/**/*.fastq" )
        .map { fastq -> tuple( file(fastq).getSimpleName(), file(fastq) ) }
	
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

    SORT_BY_NAME (
        FETCH_NT.out,
		CONSTRUCT_TAX_TREE.out
    )

    // SKETCH_BLACKLIST (
    //     SORT_BY_NAME.out
    // )

    SKETCH_WITH_BBSKETCH (
        SORT_BY_NAME.out
    )

    CLASSIFY_WITH_BBSKETCH (
        ch_fastqs,
        SKETCH_WITH_BBSKETCH.out.collect()
    )

	// SKETCH_WITH_SYLPH ()

	// CLASSIFY_WITH_SYLPH ()

	// SKETCH_WITH_SOURMASH ()

	// CLASSIFY_WITH_SOURMASH ()
	
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
	
	output:
    path "nt.fa.gz"
	
	script:
	"""
	wget -q -O - ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nt.gz \
    | gi2taxid.sh -Xmx1g \
    in=stdin.fa.gz out=nt.fa.gz \
    pigz=32 unpigz=t bgzip=t preferbgzip=t zl=8 server ow shrinknames maxbadheaders=5000 \
    badheaders=badHeaders.txt taxpath=${params.taxpath}
	"""
}

process SORT_BY_NAME {
	
	/* */

	storeDir params.nt_storedir

	memory 32.GB
	
	input:
	path nt_fasta
	path taxtree
	
	output:
	path "nt_sorted.fa.gz"
	
	script:
	"""
	sortbyname.sh -Xmx32g \
    in=`realpath ${nt_fasta}` out=nt_sorted.fa.gz \
    ow taxa tree=`realpath ${taxtree}` fastawrap=1023 zl=9 pigz=32 minlen=60 bgzip unbgzip
	"""
}

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
//     prepasses=1 tree="${params.taxpath}/tree.taxtree.gz" taxa taxlevel=genus ow mincount=120 k=32,24 depth taxpath=${params.taxpath}
// 	"""
// }

process SKETCH_WITH_BBSKETCH {
	
	/* */
	
	publishDir params.bbsketches, mode: 'copy', overwrite: true

	memory 32.GB
	
	input:
	path nt_sorted
	
	output:
	path "taxa*.sketch"

	when:
	params.download_only == false
	
	script:
	"""
	bbsketch.sh -Xmx32g 
    in=`realpath ${nt_sorted}` out=taxa#.sketch \
    mode=taxa tree="${params.taxpath}/tree.taxtree.gz" files=31 ow unpigz \
	minsize=300 prefilter autosize k=32,24 depth taxpath=${params.taxpath} # \
    # blacklist=blacklist_nt_genus_100.sketch 
	"""
}

process CLASSIFY_WITH_BBSKETCH {
	
	/* */
	
	publishDir params.bbsketch_classifications, mode: 'copy', overwrite: true

	cpus 3
	
	input:
    tuple val(sample_id), path(fastq)
	each path(nt_sorted)
	
	output:
	path "*"

	when:
	params.download_only == false
	
	script:
	"""
	comparesketch.sh \
    in=`realpath ${fastq}` k=32,24 tree="${params.taxpath}/tree.taxtree.gz" taxa*.sketch \
    exclude=1923094,Potexvirus,Virgaviridae,Bromoviridae,191289,Tymoviridae,Carlavirus # \
    # blacklist=blacklist_nt_genus_100.sketch
	"""
}

// process SKETCH_WITH_SYLPH {
	
// 	/* */
	
// 	publishDir params.sylph_sketches, mode: 'copy', overwrite: true

// }

// process CLASSIFY_WITH_SYLPH {
	
// 	/* */
	
// 	publishDir params.sylph_classifications, mode: 'copy', overwrite: true
	
// }

// process SKETCH_WITH_SOURMASH {
	
// 	/* */
	
// 	publishDir params.sourmash_sketches, mode: 'copy', overwrite: true

// }

// process CLASSIFY_WITH_SOURMASH {
	
// 	/* */
	
// 	publishDir params.sourmash_classifications, mode: 'copy', overwrite: true

// }

// --------------------------------------------------------------- //