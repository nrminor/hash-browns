process FETCH_MANIFEST_DATABASES {

	errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
	maxRetries 1

	storeDir params.db_cache

	input:
	path data_manifest

	output:
	path "*.syldb", emit: sylph_db
	path "*k51.zip", emit: sourmash_k51
	path "*k31.zip", emit: sourmash_k31
	path "*k21.zip", emit: sourmash_k21
	path "*.lineages.csv.gz", emit: sourmash_taxonomy
	path "*gottcha_db.species.*", emit: gottcha_db

	when:
	params.download || params.download_only

	script:
	"""
	sdf pull --urls --overwrite && \
	sdf status > ${params.date}_status_check.txt
	"""
}
