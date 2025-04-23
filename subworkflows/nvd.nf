include { PREPARE_NVD_INPUT                   } from "../modules/utils"
include {
    EXTRACT_HUMAN_VIRUS_READS ;
    CLASSIFY_CONTIGS_FIRST_PASS ;
    GENERATE_CONTIGS_TAXA_LIST ;
    CLASSIFY_CONTIGS_SECOND_PASS ;
    GENERATE_STAT_CONTIG_REPORT ;
    IDENTIFY_HUMAN_VIRUS_FAMILY_CONTIGS
} from "../modules/stat"
include { MASK_LOW_COMPLEXITY ; FILTER_SHORT_CONTIGS } from "../modules/bbmap"
include { RUN_SPADES                          } from "../modules/spades"
include {
    MEGABLAST ;
    ANNOTATE_MEGABLAST_RESULTS ;
    FILTER_NON_VIRUS_MEGABLAST_NODES ;
    REMOVE_MEGABLAST_MAPPED_CONTIGS ;
    BLASTN_CLASSIFY ;
    ANNOTATE_BLASTN_RESULTS ;
    FILTER_NON_VIRUS_BLASTN_NODES ;
    EXTRACT_UNCLASSIFIED_CONTIGS
} from "../modules/blast"
include { EXTRACT_HUMAN_VIRUS_CONTIGS         } from "../modules/seqkit"
include { MAP_READS_TO_CONTIGS ; COUNT_MAPPED_READS } from "../modules/minimap2"
include {
    LABKEY_UPLOAD_BLAST ;
    LABKEY_UPLOAD_FASTA ;
    CREATE_TARBALL ;
    UPLOAD_FILES_TO_LABKEY
} from "../modules/labkey"


/*

TODO:

NVD as a workflow declaration is largely recapitulated here, but a few core refactorings
remain:

    1) Refactor Python scripts so that they are not hardcoded for Snakemake. Some of the
       scripts from earlier in NVD are a nice demonstration how to support execution
       within Snakemake *and* execution through the command line.
    2) Re-implement handling of the global taxa list sqlite database.
    3) Split the workflow from below, which is too long for my taste, into smaller
       conceptual chunks in subworkflows. These workflows could be something to the tune
       of the following:
        - CONTIG_PREP
        - FIND_HUMAN_VIRUS_CONTIGS
        - COARSE_CLASSIFY_WITH_MEGABLAST
        - FINE_CLASSIFY_WITH_BLASTN
        - BUNDLE_FOR_LABKEY
    4) Runtime testing

*/

workflow NVD {
    take:
    ch_sample_fastqs // Queue channel of sample IDs and (merged) FASTQ files: tuple val(sample_id), path(fastq)
    ch_nvd_dbs

    main:
    // SPLIT OUT THE VARIOUS DATABASES FOR NVD
    // -------------------------------------------------------------------------------//
    // queue channel of a single file with a blastdb in it
    ch_blast_db = ch_nvd_dbs
        .flatten()
        .filter { db -> file(db).getPath().contains("core_nt") }

    // Queue channel of a single file, the STAT dbss
    ch_stat_dbss = ch_nvd_dbs
        .flatten()
        .filter { db -> file(db).getBaseName().endsWith(".dbss") }

    // Queue channel of a single file, the STAT index
    ch_stat_index = ch_nvd_dbs
        .flatten()
        .filter { db -> file(db).getBaseName().endsWith("dbs") }

    _ch_stat_annotation = ch_nvd_dbs
        .flatten()
        .filter { db -> file(db).getBaseName().endsWith(".annotation") }

    // Queue channel of a single file, the list of human virus taxa
    ch_human_virus_taxlist = ch_nvd_dbs
        .flatten()
        .filter { db -> file(db).getBaseName().contains("human_viruses_taxlist.txt") }

    // -------------------------------------------------------------------------------//

    // PREPARE CONTIGS
    // -------------------------------------------------------------------------------//

    PREPARE_NVD_INPUT(
        ch_sample_fastqs
    )

    EXTRACT_HUMAN_VIRUS_READS(
        PREPARE_NVD_INPUT.out.combine(ch_stat_dbss).combine(ch_human_virus_taxlist)
    )

    RUN_SPADES(
        EXTRACT_HUMAN_VIRUS_READS.out
    )

    MASK_LOW_COMPLEXITY(
        RUN_SPADES.out
    )

    FILTER_SHORT_CONTIGS(
        MASK_LOW_COMPLEXITY.out
    )

    // -------------------------------------------------------------------------------//

    // EXTRACT HUMAN VIRUSES WITH TWO PASSES OF STAT
    // -------------------------------------------------------------------------------//

    CLASSIFY_CONTIGS_FIRST_PASS(
        FILTER_SHORT_CONTIGS.out.combine(ch_stat_index)
    )

    GENERATE_CONTIGS_TAXA_LIST(
        CLASSIFY_CONTIGS_FIRST_PASS.out
    )

    CLASSIFY_CONTIGS_SECOND_PASS(
        GENERATE_CONTIGS_TAXA_LIST.out.combine(ch_stat_dbss)
    )

    GENERATE_STAT_CONTIG_REPORT(
        CLASSIFY_CONTIGS_SECOND_PASS.out
    )

    IDENTIFY_HUMAN_VIRUS_FAMILY_CONTIGS(
        CLASSIFY_CONTIGS_SECOND_PASS.out
    )

    EXTRACT_HUMAN_VIRUS_CONTIGS(
        IDENTIFY_HUMAN_VIRUS_FAMILY_CONTIGS.out.mix(
            FILTER_SHORT_CONTIGS.out
        ).groupTuple(by: 0)
    )

    // -------------------------------------------------------------------------------//

    // RUN COARSE-GRAINED CLASSIFICATION/FUZZY MATCHING WITH MEGABLAST
    // -------------------------------------------------------------------------------//

    MEGABLAST(
        EXTRACT_HUMAN_VIRUS_CONTIGS.out.combine(ch_blast_db)
    )

    ANNOTATE_MEGABLAST_RESULTS(
        MEGABLAST.out
    )

    ch_gettax_sqlite_path = ANNOTATE_MEGABLAST_RESULTS.out.taxa_sqlite

    FILTER_NON_VIRUS_MEGABLAST_NODES(
        ANNOTATE_MEGABLAST_RESULTS.out.hits
    )

    REMOVE_MEGABLAST_MAPPED_CONTIGS(
        MEGABLAST.out.mix(
            EXTRACT_HUMAN_VIRUS_CONTIGS.out
        ).groupTuple(by: 0)
    )

    // -------------------------------------------------------------------------------//

    // RUN FINE-GRAINED CLASSIFICATION/EXACT(ISH) MATCHING WITH BLASTN
    // -------------------------------------------------------------------------------//

    BLASTN_CLASSIFY(
        REMOVE_MEGABLAST_MAPPED_CONTIGS.out.combine(ch_blast_db)
    )

    ANNOTATE_BLASTN_RESULTS(
        BLASTN_CLASSIFY.out.combine(ch_gettax_sqlite_path)
    )

    FILTER_NON_VIRUS_BLASTN_NODES(
        ANNOTATE_BLASTN_RESULTS.out
    )

    ch_merged_blast_results = FILTER_NON_VIRUS_MEGABLAST_NODES.out
        .mix(
            FILTER_NON_VIRUS_BLASTN_NODES
        )
        .groupTuple(by: 0)
        .collectFile { sample_id, blast_file ->
            ["${sample_id}.txt", file(blast_file).text + '\n']
        }

    EXTRACT_UNCLASSIFIED_CONTIGS(
        EXTRACT_HUMAN_VIRUS_CONTIGS.out.mix(
            MEGABLAST.out,
            BLASTN_CLASSIFY.out,
        ).groupTuple(by: 0)
    )
    // -------------------------------------------------------------------------------//

    // BUNDLE FOR LABKEY
    // -------------------------------------------------------------------------------//

    MAP_READS_TO_CONTIGS(
        EXTRACT_HUMAN_VIRUS_CONTIGS.out.mix(
            EXTRACT_HUMAN_VIRUS_READS.out
        )
    )

    COUNT_MAPPED_READS(
        MAP_READS_TO_CONTIGS.out
    )

    LABKEY_UPLOAD_BLAST(
        ch_merged_blast_results.mix(
            MAP_READS_TO_CONTIGS.out,
            ch_sample_fastqs,
        ).groupTuple(by: 0)
    )

    LABKEY_UPLOAD_FASTA(
        EXTRACT_HUMAN_VIRUS_CONTIGS.out
    )

    CREATE_TARBALL(
        EXTRACT_HUMAN_VIRUS_CONTIGS.out.mix(
            EXTRACT_UNCLASSIFIED_CONTIGS.out,
            MAP_READS_TO_CONTIGS.out,
        )
    )

    UPLOAD_FILES_TO_LABKEY(
        CREATE_TARBALL.out
    )
}
