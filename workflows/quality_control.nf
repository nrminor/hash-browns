include { VALIDATE_SEQS  } from "../modules/seqfu"
include { READ_QC        } from "../modules/nanoq"
include { FASTQC_REPORT ; MULTIQC_REPORT } from "../modules/multiqc"

workflow QUALITY_CONTROL {
    take:
    ch_fastqs

    main:
    VALIDATE_SEQS(
        ch_fastqs
    )

    READ_QC(
        VALIDATE_SEQS.out
    )

    FASTQC_REPORT(
        READ_QC.out
    )

    MULTIQC_REPORT(
        FASTQC_REPORT.out.multiqc_data.collect()
    )

    emit:
    READ_QC.out
}
