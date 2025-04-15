process MERGE_PAIRS {

	/* */

	tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 4

	input:
	tuple val(sample_id), path(reads1), path(reads2)

	output:
	tuple val(sample_id), val("illumina"), path("${sample_id}.merged.fastq.gz")

	script:
	"""
	bbmerge.sh \
	in=`realpath ${reads1}` \
	in2=`realpath ${reads2}` \
	out=${sample_id}.merged.fastq.gz \
	outu=${sample_id}.unmerged.fastq.gz \
	threads=${task.cpus} \
	-eoom
	"""
}
process FETCH_ACCESSION2TAXID {

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


process CONSTRUCT_TAX_TREE {

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

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	output:

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
// // 	publishDir params.results, mode: 'copy'

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

process SKETCH_FAST_DB_WITH_BBSKETCH {

	/* */

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	memory 32.GB

	input:
	path fast_db

	output:
	path "human_virus_taxa.sketch"

	when:
	params.download_only == false && params.fast_mode == true && params.bbsketch == true

	script:
	"""
	bbsketch.sh -Xmx32g \
    in=`realpath ${fast_db}` \
    out=human_virus_taxa.sketch \
    k=32,24 autosize=t depth=t minsize=300 \
    prefilter=t tossjunk=t ow unpigz files=1 \
    mode=taxa taxpath=${params.taxpath} \
    tree=${params.taxpath}/tree.taxtree.gz
	"""
}

process SKETCH_DB_WITH_BBSKETCH {

	/* */

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	memory 32.GB

	input:
	path nt_fasta

	output:
	path "taxa*.sketch"

	when:
	params.download_only == false && params.fast_mode == false && params.bbsketch == true

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
// 	publishDir params.bbsketch_classifications, mode: 'copy', overwrite: true

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
	in=stdin.fastq out=${sample_id}.bbsketch.tsv \
	tree=${params.taxpath}/tree.taxtree.gz *.sketch \
	k=32,24 mode=sequence level=1 format=3 records=1 ow \
	sortbyani=t printtaxa=t printdepth=t

	csvtk sort -t -k "3:nr" -l ${sample_id}.bbsketch.tsv \
	| csvtk grep -t --ignore-case -f "Ref" -r -p virus \
	| csvtk grep -t --ignore-case -f "Ref" -r -p human \
	-o ${sample_id}human_virus_only.bbsketch.tsv.tsv
	"""
}
