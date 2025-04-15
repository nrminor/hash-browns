process SOURMASH_DB_SKETCH {

	/* */

	storeDir params.db_cache

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries

	input:
	path fasta_db

	output:
	path "*.sig.gz"

	when:
	(params.tools && params.tools.contains("sourmash")) || params.all || params.sourmash

	script:
	def ref_prefix = file(fasta_db)
		.getBaseName()
		.replace(".gz", "")
		.replace(".fna", "")
		.replace(".fasta", "")
		.replace(".fa", "")
	"""
	sourmash sketch dna -p k=31,k=51,scaled=1000,abund --singleton -f -o ${ref_prefix}.sig.gz ${fasta_db}
	"""
}

process SOURMASH_INDEX {

	/*
	Note!

	This process will require Sourmash version >= 4.9.0. As of April 15th, 2025,
	the maximum version on PyPI and Bioconda is 4.8.14.
	*/

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries

	input:
	path fasta_db

	output:
	path "*.rocksdb"

	when:
	(params.tools && params.tools.contains("sourmash")) || params.all || params.sourmash

	script:
	def ref_prefix = file(fasta_db).getBaseName()
	"""
	sourmash index \
	--ksize 31 --dna --scaled 1000 \
	${ref_prefix}.rocksdb fasta_db -F rocksdb
	"""
}

process SOURMASH_TAX_PREPARE {

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries

	input:
	path tax_csv

	output:
	path "*"

	when:
	(params.tools && params.tools.contains("sourmash")) || params.all || params.sourmash

	script:
	"""
	sourmash tax prepare --taxonomy ${tax_csv} -o tax.db
	"""
}

process SOURMASH_SKETCH_SAMPLE {

	/* */

	tag "${sample_id}"

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}_reads.sig")

	when:
	(params.tools && params.tools.contains("sourmash")) || params.all || params.sourmash

	script:
	"""
	sourmash sketch dna -p scaled=1000,k=31,k=51 ${reads} -o ${sample_id}_reads.sig
	"""
}

process SOURMASH_GATHER {

	/* */

	tag "${sample_id}"

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	input:
	tuple val(sample_id), path(sample_sigs), path(ref_sigs)

	output:
	tuple val(sample_id), path("${sample_id}_sourmash_results.csv")

	when:
	(params.tools && params.tools.contains("sourmash")) || params.all || params.sourmash

	script:
	"""
	sourmash gather \
	--prefetch --estimate-ani-ci --create-empty-results -k 31 \
	${sample_sigs} ${ref_sigs} -o ${sample_id}_sourmash_results.csv
	"""
}

process SOURMASH_TAX_METAGENOME {

	/* */

	tag "${sample_id}"

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	input:
	tuple val(sample_id), path(gather_csv), path(taxonomy)

	output:
	path "*"

	when:
	(params.tools && params.tools.contains("sourmash")) || params.all || params.sourmash

	script:
	"""
	sourmash tax metagenome
    --gather-csv ${gather_csv} \
    --taxonomy ${taxonomy}
	"""
}

process SOURMASH_TAX_SUMMARIZE {
	input:
	path db

	script:
	"""
	sourmash tax summarize gtdb-rs202.taxonomy.v2.db -o ranks.csv
	"""
}
