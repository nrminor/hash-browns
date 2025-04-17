#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DOWNLOADS       } from "./workflows/downloads"
include { CLASSIFICATION  } from "./workflows/classification"
include { QUALITY_CONTROL } from "./workflows/quality_control"
include { GATHER_ILLUMINA } from "./workflows/gather_illumina"


// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //

workflow {

	// assert params.taxpath : "Please set path on your system for storing taxonomy files."
	// assert params.nt_storedir : "Please set path on your system for storing the NCBI nt database."

	// prints to the screen and to the log
	// https://patorjk.com/software/taag/#p=testall&t=HASH-BROWNS
	log.info(
		"""
			=======================================================================================================================
			=  ====  =====  ======      ===  ====  ============      ====       =====    ====  ====  ====  ==  =======  ===      ==
			=  ====  ====    ====  ====  ==  ====  ============  ===  ===  ====  ===  ==  ===  ====  ====  ==   ======  ==  ====  =
			=  ====  ===  ==  ===  ====  ==  ====  ============  ====  ==  ====  ==  ====  ==  ====  ====  ==    =====  ==  ====  =
			=  ====  ==  ====  ===  =======  ====  ============  ===  ===  ===   ==  ====  ==  ====  ====  ==  ==  ===  ===  ======
			=        ==  ====  =====  =====        ==        ==      ====      ====  ====  ==   ==    ==  ===  ===  ==  =====  ====
			=  ====  ==        =======  ===  ====  ============  ===  ===  ====  ==  ====  ===  ==    ==  ===  ====  =  =======  ==
			=  ====  ==  ====  ==  ====  ==  ====  ============  ====  ==  ====  ==  ====  ===  ==    ==  ===  =====    ==  ====  =
			=  ====  ==  ====  ==  ====  ==  ====  ============  ===  ===  ====  ===  ==  =====    ==    ====  ======   ==  ====  =
			=  ====  ==  ====  ===      ===  ====  ============      ====  ====  ====    =======  ====  =====  =======  ===      ==
			=======================================================================================================================

			HASH-BROWNS: Metagenomic read classification well-done.
			(version 0.1.0)
			===================================
			fastq directory : ${params.fastq_dir ?: ""}
			ref fasta       : ${params.ref_fasta ?: ""}
			results dir     : ${params.results ?: ""}

			Storage directories:
			-----------------------------------
			Database cache  : ${params.db_cache}

			Chosen tools:
			-----------------------------------
			Sylph           : ${params.sylph ? "activated" : "off"}
			Sourmash        : ${params.sourmash ? "activated" : "off"}
			GOTTCHA2        : ${params.gottcha2 ? "activated" : "off"}
			Strobealign     : Coming soon!
			KMCP            : Coming soon!
			Mmseqs2         : Coming soon!
			CLARK           : Coming soon!
			MetaPHlAn4      : Coming soon!
			Centrifuge      : Coming soon!
			Megan           : Coming soon!
			mOTUs           : Coming soon!

			Run settings:
			-----------------------------------
			realtime_dir   : ${params.realtime_dir ?: ""}
			cleanup mode   : ${params.cleanup ? "on" : "off"}
			download only? : ${params.download_only ?: ""}
			available cpus : ${params.available_cpus}
			run date       : ${params.date}

			""".stripIndent()
	)

	// fastq channels
	if (params.realtime_dir) {
		ch_single_fastqs = Channel.watchPath("${params.realtime_dir}/**/*.fastq*", 'create,modify')
			.map { fastq -> tuple(file(fastq), file(fastq).countFastq()) }
			.filter { it[0].name.contains('_R1') == false && it[0].name.contains('_R2') == false }
			.filter { it[1] >= 20 }
			.map { fastq, _count -> tuple(file(fastq).getSimpleName(), "nanopore", file(fastq)) }
		ch_paired_fastqs = Channel.fromFilePairs("${params.realtime_dir}/*_{1,2}.fastq.gz", flat: true)
	}
	else {
		ch_single_fastqs = Channel.fromPath("${params.fastq_dir}/*.fastq*")
			.map { fastq -> tuple(file(fastq), file(fastq).countFastq()) }
			.filter { it[0].name.contains('_R1') == false && it[0].name.contains('_R2') == false }
			.filter { it[1] >= 20 }
			.map { fastq, _count -> tuple(file(fastq).getSimpleName(), "nanopore", file(fastq)) }
		ch_paired_fastqs = Channel.fromFilePairs("${params.fastq_dir}/*_{1,2}.fastq.gz", flat: true)
	}

	GATHER_ILLUMINA(ch_paired_fastqs)

	ch_fastqs = ch_single_fastqs.mix(GATHER_ILLUMINA.out)

	ch_custom_fa_db = params.ref_fasta
		? Channel.fromPath(params.ref_fasta)
		: Channel.empty()

	ch_data_manifest = params.download || params.download_only
		? Channel.fromPath(params.data_manifest)
		: Channel.empty()

	ch_refman_registry = params.download || params.download_only
		? Channel.fromPath(params.refman_registry)
		: Channel.empty()

	DOWNLOADS(
		ch_data_manifest,
		ch_refman_registry,
	)

	QUALITY_CONTROL(ch_fastqs)

	CLASSIFICATION(
		QUALITY_CONTROL.out,
		ch_custom_fa_db,
		DOWNLOADS.out.sylph_dbs,
		DOWNLOADS.out.sourmash_dbs,
		DOWNLOADS.out.gottcha2_db,
	)
}
