include { FETCH_MANIFEST_DATABASES } from "../modules/scidataflow"
include {
    FETCH_SYLPH_DATABASES ;
    FETCH_SOURMASH_DATABASES ;
    FETCH_GOTTCHA2_DATABASES ;
    FETCH_NVD_DATABASES ;
    EXTRACT_DB_FILES
} from "../modules/refman"

workflow DOWNLOADS {
    take:
    _ch_data_manifest
    ch_refman_registry

    main:
    FETCH_SYLPH_DATABASES(ch_refman_registry)

    FETCH_SOURMASH_DATABASES(ch_refman_registry)

    FETCH_GOTTCHA2_DATABASES(ch_refman_registry)

    FETCH_NVD_DATABASES(ch_refman_registry)

    EXTRACT_DB_FILES(
        FETCH_SYLPH_DATABASES.out.map { dbs -> tuple("sylph", dbs) }.mix(
            FETCH_SOURMASH_DATABASES.out.map { dbs -> tuple("sourmash", dbs) }
        ).mix(
            FETCH_GOTTCHA2_DATABASES.out.map { dbs -> tuple("gottcha2", dbs) }
        ).mix(
            FETCH_NVD_DATABASES.out.map { dbs -> tuple("nvd", dbs) }
        )
    )

    ch_sylph_dbs = params.download && !params.download_only
        ? EXTRACT_DB_FILES.out.filter { tool, _files -> tool == "sylph" }.map { _tool, files -> files }.flatten()
        : Channel.empty()

    ch_sourmash_dbs = params.download && !params.download_only
        ? EXTRACT_DB_FILES.out.filter { tool, _files -> tool == "sourmash" }.map { _tool, files -> files }
        : Channel.empty()

    ch_gottcha2_db = params.download && !params.download_only
        ? EXTRACT_DB_FILES.out.filter { tool, _files -> tool == "gottcha2" }.map { _tool, files -> files }
        : Channel.empty()

    ch_nvd_db = params.download && !params.download_only
        ? EXTRACT_DB_FILES.out.filter { tool, _files -> tool == "nvd" }.map { _tool, files -> files }.flatten()
        : Channel.empty()

    emit:
    sylph_dbs    = ch_sylph_dbs
    sourmash_dbs = ch_sourmash_dbs
    gottcha2_db  = ch_gottcha2_db
    nvd_dbs      = ch_nvd_db
}
