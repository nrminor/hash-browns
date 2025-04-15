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
