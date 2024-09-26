#!/usr/bin/env nextflow 
/*
 * pipeline input parameters
 */
params.run = "06_09_24"
params.shared = "/export/home/public/agletdinov_shared"
params.results_project = "/export/home/agletdinov/work/nextflow_projects/total_seq"
params.reads = '/export/data/public/lgi/seq/nanopore160/2024/2024_09_06_Nanopore_160_Midnight/2/pod5_pass/barcode*/*.pod5'
params.outdir = "${params.shared}/nextflow_projects/total_seq/results/${params.run}"
//params.maxForks = 50  // Задайте необходимое максимальное число процессов
params.dorado_corrected_models = "/export/home/public/agletdinov_shared/dorado_corrected_models/herro-v1"


log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    reads          : ${params.reads}
    outdir         : ${params.outdir}
    maxForks       : ${params.maxForks}
    """
    .stripIndent()

include { NANOPORE } from './modules/nanopore/nanopore_modules.nf'

include { MULTIQC } from './modules/multiqc.nf'
include { SENDMAIL; SENDMAIL_PY } from './modules/sendMail.nf'


workflow nanopore{
    //Channel
    //    .fromPath(params.reads)
    //    .map { file -> [file.parent.name, file] }
    //    .groupTuple()
    //    .set { nanopore_pod5_ch }
    Channel
        .fromPath(params.reads)
        .map { file -> [file.parent.name, file.parent] }  // Получаем имя директории (например, barcodeN) и её полный путь
        .distinct()  // Убираем дублирующиеся записи для одной директории
        .set { nanopore_pod5_ch }
    NANOPORE(nanopore_pod5_ch)
    MULTIQC(NANOPORE.out)
    SENDMAIL_PY(MULTIQC.out)
}


//workflow.onComplete {
//    log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
//