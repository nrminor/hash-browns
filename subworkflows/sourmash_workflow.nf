include {
    SOURMASH_DB_SKETCH ;
    SOURMASH_SKETCH_SAMPLE ;
    SOURMASH_GATHER ;
    SOURMASH_TAX_PREPARE ;
    SOURMASH_TAX_METAGENOME
} from "../modules/sourmash"

workflow SOURMASH_WORKFLOW {
    take:
    ch_query_fastqs
    ch_custom_fa_db
    ch_sourmash_k51
    ch_sourmash_k31
    ch_sourmash_k21
    ch_sourmash_taxonomy

    main:
    // SOURMASH_TAX_PREPARE(ch_sourmash_taxonomy)

    SOURMASH_DB_SKETCH(ch_custom_fa_db)

    SOURMASH_SKETCH_SAMPLE(ch_query_fastqs)

    ch_all_dbs = SOURMASH_DB_SKETCH.out.mix(ch_sourmash_k51).mix(ch_sourmash_k31).mix(ch_sourmash_k21)

    SOURMASH_GATHER(SOURMASH_SKETCH_SAMPLE.out.combine(ch_all_dbs))

    SOURMASH_TAX_METAGENOME(
        SOURMASH_GATHER.out.combine(ch_sourmash_taxonomy)
    )
}
