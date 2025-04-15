process SOURMASH_DB_SKETCH {

	/* */


	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries

	input:
	path fasta_db

	output:
	path "*.sig.gz"

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
	script:
	"""
	sourmash index <database_name>.rocksdb <inputfile1> [ <inputfile2> ... ] -F rocksdb
	"""
}

process SOURMASH_TAX_PREPARE {
	script:
	"""
	sourmash tax prepare --taxonomy file1.csv file2.csv -o tax.db
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
	params.download_only == false && params.sourmash == true

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
	path "*"

	script:
	"""
	sourmash gather \
	--prefetch --estimate-ani-ci --create-empty-results -k 31 \
	${sample_sigs} ${ref_sigs} -o ${sample_id}_sourmash_results.csv
	"""
}

process SOURMASH_TAX_METAGENOME {
	script:
	"""
	sourmash tax metagenome
    --gather-csv HSMA33MX_gather_x_gtdbrs202_k31.csv \
    --taxonomy gtdb-rs202.taxonomy.v2.csv
	"""
}

process SOURMASH_TAX_SUMMARIZE {
	script:
	"""
	sourmash tax summarize gtdb-rs202.taxonomy.v2.db -o ranks.csv
	"""
}
