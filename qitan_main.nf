#!/usr/bin/env nextflow 
/*
 * pipeline input parameters
 */
params.run = "21_01_25"
params.shared = "/export/home/public/agletdinov_shared"
params.results_project = "/export/home/agletdinov/work/nextflow_projects/qitan"
params.reads = "${params.results_project}/fastq/${params.run}/*.fastq"
params.outdir = "${params.results_project}/results/${params.run}"
//params.maxForks = 50  // Задайте необходимое максимальное число процессов
params.singleEnd = true

params.genome = "${params.shared}/genomes/vir/human_respiratory_syncytial_virus.fasta"

log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    reads          : ${params.reads}
    outdir         : ${params.outdir}
    maxForks       : ${params.maxForks}
    kraken2db      : ${params.kraken2db}
    """
    .stripIndent()

include { QITAN } from './modules/qitan/qitan_modules.nf'

include { MULTIQC } from './modules/multiqc.nf'


workflow chitan{
    Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .set { reads_ch }
    genome = file(params.genome)
    reads_ch.view()
    QITAN(reads_ch, genome)
    MULTIQC(QITAN.out)
    //SENDMAIL_PY(MULTIQC.out)
}


//workflow.onComplete {
//    log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
//