include { FETCH_MANIFEST_DATABASES } from "../modules/scidataflow"

workflow DOWNLOADS {
    take:
    ch_data_manifest

    main:
    FETCH_MANIFEST_DATABASES(ch_data_manifest)

    ch_sylph_dbs = params.download || params.download_only
        ? FETCH_MANIFEST_DATABASES.out.sylph_db
        : Channel.empty()

    ch_sourmash_k51 = params.download || params.download_only
        ? FETCH_MANIFEST_DATABASES.out.sourmash_k51
        : Channel.empty()

    ch_sourmash_k31 = params.download || params.download_only
        ? FETCH_MANIFEST_DATABASES.out.sourmash_k31
        : Channel.empty()

    ch_sourmash_k21 = params.download || params.download_only
        ? FETCH_MANIFEST_DATABASES.out.sourmash_k21
        : Channel.empty()

    ch_sourmash_taxonomy = params.download || params.download_only
        ? FETCH_MANIFEST_DATABASES.out.sourmash_taxonomy
        : Channel.empty()

    ch_gottcha2_db = params.download || params.download_only
        ? FETCH_MANIFEST_DATABASES.out.gottcha_db
        : Channel.empty()

    emit:
    sylph_db          = ch_sylph_dbs
    sourmash_k51      = ch_sourmash_k51
    sourmash_k31      = ch_sourmash_k31
    sourmash_k21      = ch_sourmash_k21
    sourmash_taxonomy = ch_sourmash_taxonomy
    gottcha2_db       = ch_gottcha2_db
}
