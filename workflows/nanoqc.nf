/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { UCSC_GTFTOGENEPRED          } from '../modules/nf-core/ucsc/gtftogenepred/main'
include { SAMTOOLS_STATS              } from '../modules/nf-core/samtools/stats/main'
include { NANOPLOT                    } from '../modules/nf-core/nanoplot/main'
include { SAMTOOLS_VIEW               } from '../modules/nf-core/samtools/view/main'
include { PICARD_COLLECTRNASEQMETRICS } from '../modules/nf-core/picard/collectrnaseqmetrics/main'
include { RSEQC_BAMSTAT               } from '../modules/nf-core/rseqc/bamstat/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap            } from 'plugin/nf-schema'
include { paramsSummaryMultiqc        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText      } from '../subworkflows/local/utils_nfcore_nanoqc_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NANOQC {
    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    // read filtering on mapq and raw read length
    if (params.filter_reads) {
        ch_filtered = SAMTOOLS_VIEW(
            ch_samplesheet,
            [[], []],
            [],
            "bai",
        )
        ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions)
        ch_bam = ch_filtered.bam
    }
    else {
        ch_bam = ch_samplesheet.map { meta, bam, bai -> tuple(meta, bam) }
    }
    // collect stats on bam files
    SAMTOOLS_STATS(ch_samplesheet, [[id: 'genome'], params.fasta])
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats.collect { it[1] })
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)
    NANOPLOT(ch_bam)
    ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT.out.txt.collect { it[1] })
    ch_versions = ch_versions.mix(NANOPLOT.out.versions)
    RSEQC_BAMSTAT(ch_bam)
    ch_multiqc_files = ch_multiqc_files.mix(RSEQC_BAMSTAT.out.txt.collect { it[1] })
    ch_versions = ch_versions.mix(RSEQC_BAMSTAT.out.versions)
    if (params.protocol == 'RNA') {
        UCSC_GTFTOGENEPRED([[id: 'genome'], params.gtf])
        ref_flat = UCSC_GTFTOGENEPRED.out.refflat.collect { it[1] }
        PICARD_COLLECTRNASEQMETRICS(ch_bam, ref_flat, params.fasta, [])
        ch_multiqc_files = ch_multiqc_files.mix(PICARD_COLLECTRNASEQMETRICS.out.metrics.collect { it[1] })
        ch_versions = ch_versions.mix(PICARD_COLLECTRNASEQMETRICS.out.versions)
    }
    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nanoqc_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = Channel.fromPath(
        "${projectDir}/assets/multiqc_config.yml",
        checkIfExists: true
    )
    ch_multiqc_custom_config = params.multiqc_config
        ? Channel.fromPath(params.multiqc_config, checkIfExists: true)
        : Channel.empty()
    ch_multiqc_logo = params.multiqc_logo
        ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : Channel.empty()

    summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )
    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true,
        )
    )

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
