process FETCH_SYLPH_DATABASES {

    errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
    maxRetries 1

    storeDir params.db_cache

    input:
    path refman_registry

    output:
    tuple val("sylph"), path("syldb.tar")

    when:
    (params.download || params.download_only) && params.sylph

    script:
    """
    refman download sylph
	"""
}

process FETCH_SOURMASH_DATABASES {

    errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
    maxRetries 1

    storeDir params.db_cache

    input:
    path refman_registry

    output:
    tuple val("sourmash"), path("sourmash-bacterial-viral.tar")

    when:
    (params.download || params.download_only) && params.sourmash

    script:
    """
    refman download sourmash
	"""
}

process FETCH_GOTTCHA2_DATABASES {

    errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
    maxRetries 1

    storeDir params.db_cache

    input:
    path refman_registry

    output:
    tuple val("gottcha2"), path("gottcha-db.tar")

    when:
    (params.download || params.download_only) && params.gottcha2

    script:
    """
    refman download gottcha2
	"""
}

process EXTRACT_DB_FILES {
    tag "${tool}"
    errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
    maxRetries 1

    storeDir params.db_cache

    input:
    tuple val(tool), path(tar_archive)

    output:
    tuple val(tool), path("*")

    script:
    """
    tar xf ${tar_archive}
    """
}
