include { SYLPH_WORKFLOW       } from "../subworkflows/sylph_workflow"
include { SOURMASH_WORKFLOW    } from "../subworkflows/sourmash_workflow"
include { GOTTCHA2_WORKFLOW    } from "../subworkflows/gottcha2_workflow"
include { KMCP_WORKFLOW        } from "../subworkflows/kmcp_workflow"
include { STROBEALIGN_WORKFLOW } from "../subworkflows/strobealign_workflow"
include { NVD                  } from "../subworkflows/nvd"

workflow CLASSIFICATION {
    take:
    ch_query_fastqs
    ch_custom_fa_db
    ch_sylph_dbs
    _ch_sourmash_dbs
    ch_gottcha2_db
    ch_nvd_dbs

    main:
    SYLPH_WORKFLOW(
        ch_query_fastqs.map { id, _platform, fastq -> tuple(id, file(fastq)) },
        ch_sylph_dbs.flatten(),
        ch_custom_fa_db,
    )

    // SOURMASH_WORKFLOW(
    //     ch_query_fastqs.map { id, _platform, fastq -> tuple(id, file(fastq)) },
    //     ch_custom_fa_db,
    //     ch_sourmash_dbs,
    // )

    GOTTCHA2_WORKFLOW(
        ch_gottcha2_db,
        ch_query_fastqs,
    )

    NVD(
        ch_query_fastqs,
        ch_nvd_dbs,
    )
}
