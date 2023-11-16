#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {
	
	
	// input channels
    ch_fastqs = Channel
        .fromPath( "${params.fasq_dir}/**/*.fastq" )
        .map { fastq -> tuple( file(fastq).getSimpleName(), file(fastq) ) }
	
	
	// Workflow steps 
    FETCH_NT ()

    SORT_BY_NAME (
        FETCH_NT.out
    )

    // SKETCH_BLACKLIST (
    //     SORT_BY_NAME.out
    // )

    GENERATE_SKETCHES (
        SORT_BY_NAME.out
    )

    CLASSIFY_WITH_SKETCHES (
        ch_fastqs,
        GENERATE_SKETCHES.out.collect()
    )
	
	
}
// --------------------------------------------------------------- //



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //
// Additional parameters that are derived from parameters set in nextflow.config

// --------------------------------------------------------------- //




// PROCESS SPECIFICATION 
// --------------------------------------------------------------- //

process FETCH_NT {
	
	/* */
	
	output:
    path "${params.date}_nt.fa.gz"
	
	when:
    params.download_nt == true
	
	script:
	"""
	wget -q -O - ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nt.gz \
    | gi2taxid.sh -Xmx1g \
    in=stdin.fa.gz out=${params.date}_nt.fa.gz \
    pigz=32 unpigz bgzip zl=8 server ow shrinknames maxbadheaders=5000 \
    badheaders=badHeaders.txt taxpath=${params.taxpath}
	"""
}

process SORT_BY_NAME {
	
	/* */
	
	input:
	path nt_fasta
	
	output:
	path "${params.date}_nt_sorted.fa.gz"
	
	script:
	"""
	sortbyname.sh -Xmx96g \
    in=`realpath ${nt_fasta}` out=${params.date}_nt_sorted.fa.gz \
    ow taxa tree=auto fastawrap=1023 zl=9 pigz=32 minlen=60 bgzip unbgzip
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
//     prepasses=1 tree=auto taxa taxlevel=genus ow mincount=120 k=32,24 depth taxpath=${params.taxpath}
// 	"""
// }

process GENERATE_SKETCHES {
	
	/* */
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	input:
	path nt_sorted
	
	output:
	path "taxa*.sketch"
	
	script:
	"""
	bbsketch.sh -Xmx31g 
    in=`realpath ${nt_sorted}` out=taxa#.sketch \
    mode=taxa tree=auto files=31 ow unpigz minsize=300 prefilter autosize k=32,24 depth taxpath=${params.taxpath} # \
    # blacklist=blacklist_nt_genus_100.sketch 
	"""
}

process CLASSIFY_WITH_SKETCHES {
	
	/* */
	
	tag "${sample_id}"
	publishDir params.results, mode: 'copy'
	
	input:
    tuple val(sample_id), path(fastq)
	each path(nt_sorted)
	
	output:
	path "taxa*.sketch"
	
	script:
	"""
	comparesketch.sh \
    in=`realpath ${fastq}` k=32,24 tree=auto taxa*.sketch \
    exclude=1923094,Potexvirus,Virgaviridae,Bromoviridae,191289,Tymoviridae,Carlavirus # \
    # blacklist=blacklist_nt_genus_100.sketch
	"""
}

// --------------------------------------------------------------- //