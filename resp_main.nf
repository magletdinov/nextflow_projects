#!/usr/bin/env nextflow 
/*
 * pipeline input parameters
 */
params.shared = "/export/home/public/agletdinov_shared"
params.results_project = "/export/home/agletdinov/work/nextflow_projects/resp_panel"
params.reads = "${params.results_project}/fastq/25_04_24/*R{1,2}*.fastq.gz"
params.adapters = "${params.shared}/adapters/adapters.fasta"
params.primer_a = "${params.shared}/primers/resp_panel/resp_primers.fasta"
params.primer_g = "${params.shared}/primers/resp_panel/resp_primers.fasta"
params.kraken2db = "/export/home/public/agletdinov_shared/kraken2db/minikraken2_v2_8GB_201904_UPDATE"
params.outdir = "${params.results_project}/results_kraken2/25_04_24"

//params.maxForks = 50  // Задайте необходимое максимальное число процессов

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
    """
    .stripIndent()

include { TAXONOMY_ANALYSIS } from './modules/resp_panel/resp_modules.nf'
include { MULTIQC } from './modules/multiqc.nf'

workflow taxonomy_analysis{
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }
    TAXONOMY_ANALYSIS(read_pairs_ch)
    MULTIQC(TAXONOMY_ANALYSIS.out)
}
//workflow.onComplete {
//    log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
//

