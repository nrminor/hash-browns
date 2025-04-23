include {
    SYLPH_SKETCH_DB ;
    SYLPH_SKETCH_SAMPLE ;
    CLASSIFY_WITH_SYLPH ;
    SYLPH_TAX_DOWNLOAD ;
    SYLPH_TAXPROF
} from "../modules/sylph"

workflow SYLPH_WORKFLOW {
    take:
    ch_sample_fastqs
    ch_sylph_db
    ch_custom_fa

    main:
    SYLPH_SKETCH_DB(ch_custom_fa)

    SYLPH_SKETCH_SAMPLE(ch_sample_fastqs)

    ch_all_databases = SYLPH_SKETCH_DB.out.mix(ch_sylph_db)

    CLASSIFY_WITH_SYLPH(
        SYLPH_SKETCH_SAMPLE.out.combine(ch_all_databases)
    )

    SYLPH_TAX_DOWNLOAD()
}
