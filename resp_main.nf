#!/usr/bin/env nextflow 
/*
 * pipeline input parameters
 */
params.run = "22_08_25"
params.shared = "/export/home/public/agletdinov_shared"
params.results_project = "/export/home/agletdinov/work/nextflow_projects/resp_panel"
params.reads = "${params.results_project}/fastq/${params.run}/*R{1,2}*.fastq.gz"
params.adapters = "${params.shared}/adapters/adapters.fasta"
//params.primer_a = "${params.shared}/primers/resp_panel/resp_primers.fasta"
//params.primer_g = "${params.shared}/primers/resp_panel/resp_primers.fasta"
//params.kraken2db = "/export/home/public/agletdinov_shared/kraken2db/minikraken2_v2_8GB_201904_UPDATE"
params.outdir = "${params.results_project}/results/${params.run}"

//params.maxForks = 50  // Задайте необходимое максимальное число процессов
params.old_panel_fasta = "${params.shared}/genomes/resp_panels/02_2025/resp_panel-02_2025.fasta"
params.new_panel_fasta = "${params.shared}/genomes/resp_panels/07_2025/resp_panel-07_2025.fasta"

singleEnd_TRIM_ADAPT = false
params.singleEnd_TRIM_4_NUCL = false
params.singleEnd_KRAKEN2 = false
params.singleEnd_BWA_MEM_BAM_SORT = false


log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    reads          : ${params.reads}
    adapters       : ${params.adapters}
    primer_a       : ${params.primer_a}
    primer_g       : ${params.primer_g}
    kraken2db      : ${params.kraken2db}
    outdir         : ${params.outdir}
    maxForks       : ${params.maxForks}
    old_panel_fasta: ${params.old_panel_fasta}
    new_panel_fasta: ${params.new_panel_fasta}
    """
    .stripIndent()

include { TAXONOMY_ANALYSIS } from './modules/resp_panel/resp_modules.nf'
include { TAXONOMY_ANALYSIS_WITH_REF_ALIGN } from './modules/resp_panel/resp_modules.nf'
include { TAXONOMY_ANALYSIS_WITH_REF_ALIGN_BASE } from './modules/resp_panel/resp_modules.nf'
include { MULTIQC } from './modules/multiqc.nf'

workflow taxonomy_analysis{
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }
    TAXONOMY_ANALYSIS(read_pairs_ch)
    MULTIQC(TAXONOMY_ANALYSIS.out)
}

workflow taxonomy_analysis_with_ref_align_base{
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

    TAXONOMY_ANALYSIS_WITH_REF_ALIGN_BASE(read_pairs_ch, params.new_panel_fasta)
    MULTIQC(TAXONOMY_ANALYSIS_WITH_REF_ALIGN_BASE.out)
}

workflow taxonomy_analysis_with_ref_align{
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

    TAXONOMY_ANALYSIS_WITH_REF_ALIGN(read_pairs_ch, params.old_panel_fasta, params.new_panel_fasta)
    MULTIQC(TAXONOMY_ANALYSIS_WITH_REF_ALIGN.out)
}

//workflow.onComplete {
//    log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
//

