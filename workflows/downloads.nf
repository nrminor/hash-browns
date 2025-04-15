include { FETCH_MANIFEST_DATABASES } from "../modules/scidataflow"

workflow DOWNLOADS {
    take:
    ch_data_manifest

    main:
    FETCH_MANIFEST_DATABASES(ch_data_manifest)

    emit:
    sylph_db          = FETCH_MANIFEST_DATABASES.out.sylph_db
    sourmash_k51      = FETCH_MANIFEST_DATABASES.out.sourmash_k51
    sourmash_k31      = FETCH_MANIFEST_DATABASES.out.sourmash_k31
    sourmash_k21      = FETCH_MANIFEST_DATABASES.out.sourmash_k21
    sourmash_taxonomy = FETCH_MANIFEST_DATABASES.out.sourmash_taxonomy
}
