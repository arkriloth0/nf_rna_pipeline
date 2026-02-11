/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// base
include { FASTQC                } from '../modules/nf-core/fastqc/main'
include { MULTIQC               } from '../modules/nf-core/multiqc/main'
// trimming
include { FASTP                 } from '../modules/nf-core/fastp/main'
// alignment
include { STAR_GENOMEGENERATE   } from '../modules/nf-core/star/genomegenerate/main'
include { STAR_ALIGN            } from '../modules/nf-core/star/align/main'
// quantification
include { SALMON_INDEX          } from '../modules/nf-core/salmon/index/main'
include { SALMON_QUANT          } from '../modules/nf-core/salmon/quant/main'
// genotyping
include { GATK4_CREATESEQUENCEDICTIONARY } from '../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_HAPLOTYPECALLER          } from '../modules/nf-core/gatk4/haplotypecaller/main'
// utils
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_nf_rna_pipeline_pipeline'
include { SAMTOOLS_FAIDX         } from '../modules/nf-core/samtools/faidx/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NF_RNA_PIPELINE {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()
    //
    // MODULE: Run FastQC -default is skipped with `skip_tools`
    //
    FASTQC(
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { it[1] })
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // PARSE INPUT: Parse samplesheet to separate FASTQ and BAM inputs
    //
    ch_samplesheet
        .branch { meta, fastq_bam ->
            fastq: meta.data_type == 'fastq'
            bam: meta.data_type == 'bam'
        }
        .set { ch_input }

    //
    // GROUP FASTQ FILES BY SAMPLE ID (handles multi-lane samples)
    //
    ch_input.fastq
        .map { meta, reads ->
            def group_key = meta.id
            [ group_key, meta, reads ]
        }
        .groupTuple(by: 0)
        .map { group_key, metas, reads_list ->
            // Validate consistency across lanes
            def single_end_values = metas.collect { it.single_end }.unique()
            if (single_end_values.size() > 1) {
                error("Sample '${group_key}' has inconsistent single_end values across lanes. All lanes must be either single-end or paired-end.")
            }

            def meta = metas[0]
            def all_reads = reads_list.flatten()

            if (!meta.single_end) {
                all_reads = all_reads.sort { it.name }
            }

            log.info "Sample '${meta.id}': merged ${metas.size()} lane(s), ${all_reads.size()} file(s)"

            [ meta, all_reads ]
        }
        .set { ch_fastq_grouped }

    //
    // TRIMMING: fastp for poly-X tail trimming, 3' quality trimming, quality filtering (Q>=15)
    // The 'trim' column in the samplesheet controls per-sample trimming (default: true)
    //

    // Branch into samples that need trimming vs those that don't
    ch_fastq_grouped
        .branch { meta, reads ->
            trim:    meta.trim == true || meta.trim == 'TRUE' || meta.trim == 'true' || meta.trim == 'True'
            no_trim: true
        }
        .set { ch_trim_branch }

    // Run FASTP on samples marked for trimming (with conditional version/multiqc collection)
    ch_trim_count = ch_trim_branch.trim.count()
    
    FASTP(
        ch_trim_branch.trim.map { meta, reads -> [ meta, reads, [] ] },  // tuple val(meta), path(reads), path(adapter_fasta)
        false,                  // val discard_trimmed_pass
        false,                  // val save_trimmed_fail
        false                   // val save_merged
    )
    
    // Only collect multiqc files if FASTP actually produced output
    ch_multiqc_files = ch_multiqc_files.mix(
        FASTP.out.json.ifEmpty([]).collect { it[1] }
    )

    // Merge trimmed and untrimmed channels back together for downstream steps
    ch_fastq_for_alignment = FASTP.out.reads
        .mix(ch_trim_branch.no_trim)

    //
    // VALIDATE BAM INPUT
    //
    ch_input.bam
        .count()
        .subscribe { count ->
            if (count > 0) {
                if (!params.bam_type) {
                    error "BAM input detected but --bam_type not specified. Please use --bam_type 'genome' or --bam_type 'transcriptome'"
                }
                if (!['genome', 'transcriptome'].contains(params.bam_type)) {
                    error "Invalid --bam_type '${params.bam_type}'. Must be 'genome' or 'transcriptome'"
                }
                if (params.bam_type == 'genome') {
                    log.info "BAM input with --bam_type 'genome': Will run variant calling, skipping Salmon quantification"
                } else {
                    log.info "BAM input with --bam_type 'transcriptome': Will run Salmon quantification, skipping variant calling"
                }
            }
        }

    //
    // ROUTE BAM INPUT BASED ON TYPE
    //
    ch_input.bam
        .branch { meta, bam ->
            genome: params.bam_type == 'genome'
            transcriptome: params.bam_type == 'transcriptome'
        }
        .set { ch_bam_typed }

    //
    // REFERENCE MANAGEMENT
    //

    // Compute reference directory
    // def reference_dir = params.reference_dir ?: (params.fasta ? file(params.fasta).parent : "${params.outdir}/references")

    if (params.fasta) {
        ch_fasta = Channel.fromPath(params.fasta).map { [ [:], it ] }.first()

        // Check if .fai exists, generate if not
        def fai_path = params.fasta + '.fai'
        if (file(fai_path).exists()) {
            ch_fai = Channel.fromPath(fai_path).map { [ [:], it ] }.first()
        } else {
            log.info 'FASTA index (.fai) not found, generating from reference genome'
                SAMTOOLS_FAIDX(
                ch_fasta,           // tuple val(meta), path(fasta)
                [[], []], // tuple val(meta2), path(fai) - empty since we're generating it
                false               // val get_sizes - set to false (we just want the .fai file)
            )
            ch_fai = SAMTOOLS_FAIDX.out.fai.first()
            ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())
        }

        // Check if .dict exists, generate if not
        def dict_path = params.fasta.replaceAll(/\.fa(sta)?$/, '.dict')
        if (file(dict_path).exists()) {
            ch_dict = Channel.fromPath(dict_path).map { [ [:], it ] }.first()
        } else {
            log.info 'Sequence dictionary (.dict) not found, generating from reference genome'
            GATK4_CREATESEQUENCEDICTIONARY(ch_fasta)
            ch_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict.first()
            ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions.first())
        }
    }

    //
    // ALIGNMENT (if FQs provided as input)
    //
    // STAR index
    ch_gtf = Channel.fromPath(params.gtf).map { [ [:], it ] }.first()
    if (!params.star_index && params.fasta && params.gtf) {
        ch_fasta = Channel.fromPath(params.fasta).map { [ [:], it ] }.first()

        STAR_GENOMEGENERATE(
            ch_fasta,
            ch_gtf
        )

        ch_star_index = STAR_GENOMEGENERATE.out.index.first()
        ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions.first())
    } else if (params.star_index) {
        ch_star_index = Channel.fromPath(params.star_index).map { [ [:], it ] }.first()
    }

    STAR_ALIGN(
        ch_fastq_for_alignment,                  // tuple val(meta), path(reads)
        ch_star_index,                       // tuple val(meta2), path(index)
        ch_gtf,                              // tuple val(meta3), path(gtf)
        params.salmon_star_ignore_sjdbgtf,  // val star_ignore_sjdbgtf
        params.salmon_seq_platform ?: '',           // val seq_platform
        params.salmon_seq_center ?: ''              // val seq_center
    )

    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())

    // Transcriptome BAMs: STAR output + user-provided transcriptome BAMs
    ch_transcriptome_bam = STAR_ALIGN.out.bam_transcript
        .mix(ch_bam_typed.transcriptome)

    // Genome BAMs: STAR output + user-provided genome BAMs
    ch_genome_bam = STAR_ALIGN.out.bam_sorted_aligned
        .mix(ch_bam_typed.genome)

    //
    // QUANTIFICATION with salmon
    //
    // Create salmon index if needed
    ch_salmon_index = Channel.empty()
    if (params.transcriptome && !params.salmon_index) {
        ch_transcriptome = Channel.fromPath(params.transcriptome)
        SALMON_INDEX(ch_fasta.map { meta, fa -> [fa] }, ch_transcriptome)
        ch_salmon_index = SALMON_INDEX.out.index.first()
        ch_versions = ch_versions.mix(SALMON_INDEX.out.versions.first())
    } else if (params.salmon_index) {
        ch_salmon_index = Channel.fromPath(params.salmon_index).first()
    }

    ch_transcript_fasta = params.transcriptome ?
        Channel.fromPath(params.transcriptome).first() : Channel.empty()

    if (params.salmon_quant_mode.contains('alignment') && ch_transcriptome_bam) {
        // Alignment mode with BAM files
        SALMON_QUANT(
            ch_transcriptome_bam.map { meta, bam -> [ meta, [bam] ] }, // tuple val(meta), path(reads)
            ch_salmon_index,                                   // path index
            ch_gtf.map { meta, gtf -> gtf },                  // path gtf
            ch_transcript_fasta,                              // path transcript_fasta
            true,                                             // val alignment_mode
            'A'                                               // val lib_type
        )
    } else {
        // Mapping mode with FASTQ files
        SALMON_QUANT(
            ch_fastq_for_alignment,                        // tuple val(meta), path(reads)
            ch_salmon_index,                             // path index
            ch_gtf.map { meta, gtf -> gtf },             // path gtf
            ch_transcript_fasta,                         // path transcript_fasta
            false,                                       // val alignment_mode
            'A'                                          // val lib_type
        )
    }
    ch_versions = ch_versions.mix(SALMON_QUANT.out.versions.first())

    //
    // GENOTYPING with GATK's HaplotypeCaller
    //

    // Prepare input for GATK4_HAPLOTYPECALLER
    ch_bam_for_hc = ch_genome_bam.map { meta, bam ->
        // nf-core module expects: tuple val(meta), path(input), path(input_index), path(intervals), path(dragstr_model)
        [ meta, bam, [], [], [] ]
    }

    GATK4_HAPLOTYPECALLER(
        ch_bam_for_hc,                           // tuple val(meta), path(input), path(input_index), path(intervals), path(dragstr_model)
        ch_fasta,                                // tuple val(meta2), path(fasta)
        ch_fai,                                  // tuple val(meta3), path(fai)
        ch_dict,                                 // tuple val(meta4), path(dict)
        Channel.empty().map { [ [:], [] ] },     // tuple val(meta5), path(dbsnp)
        Channel.empty().map { [ [:], [] ] }      // tuple val(meta6), path(dbsnp_tbi)
    )

    ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions.first())

    //
    // Collate and save software versions
    //
    def topic_versions = Channel.topic('versions')
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':') + 1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'nf_rna_pipeline_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: 'nextflow_schema.json')
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
