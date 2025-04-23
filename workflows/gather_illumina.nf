include { MERGE_PAIRS } from "../modules/bbmap"

workflow GATHER_ILLUMINA {
    take:
    ch_paired_fastqs

    main:
    MERGE_PAIRS(ch_paired_fastqs).map { id, fastq -> tuple(id, "illumina", file(fastq)) }

    emit:
    MERGE_PAIRS.out
}
