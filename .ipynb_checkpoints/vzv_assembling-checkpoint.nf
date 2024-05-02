#!/usr/bin/env nextflow 
/*
 * pipeline input parameters
 */
//params.reads = "$projectDir/test/vzv/fastq/146*R{1,2}*.fastq.gz"
params.reads = "/export/home/agletdinov/work/nextflow_projects/VZV/fastq/25_04_24/*R{1,2}*.fastq.gz"
params.adapters = "$projectDir/test/vzv/adapters/adapters.fasta"
params.genome = "/export/home/public/agletdinov_shared/genomes/vzv/NC_001348.1.fasta"
params.primer_a = "$projectDir/test/vzv/primers/aADAPTERX_v2.fasta"
params.primer_g = "$projectDir/test/vzv/primers/gXADAPTER_v2.fasta"
params.key_areas = "$projectDir/test/vzv/key_areas/vzv_amplicons.fasta"
params.script = "$projectDir/bin/rename_fasta_id.py"
params.kraken2db = "/export/home/public/agletdinov_shared/kraken2db/minikraken2_v2_8GB_201904_UPDATE"
params.outdir = "/export/home/agletdinov/work/nextflow_projects/VZV/results/25_04_24"
//params.outdir = "check_raw_data/"

//params.maxForks = 50  // Задайте необходимое максимальное число процессов

log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    reads          : ${params.reads}
    adapters       : ${params.adapters}
    primer_a       : ${params.primer_a}
    primer_g       : ${params.primer_g}
    key_areas      : ${params.key_areas}
    script         : ${params.script}
    kraken2db      : ${params.kraken2db}
    outdir         : ${params.outdir}
    maxForks       : ${params.maxForks}
    """
    .stripIndent()

include { VZV } from './vzv_modules.nf'
include { VZV_CHECK } from './vzv_modules.nf'
include { MULTIQC } from './multiqc.nf'

workflow vzv_full_run{
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }
    key_areas = file(params.key_areas)
    script = params.script
    VZV(read_pairs_ch, key_areas, script)
    MULTIQC(VZV.out)
}

workflow vzv_check_raw{
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }
    key_areas = file(params.genome)
    script = params.script
    VZV_CHECK(read_pairs_ch, key_areas, script)
    MULTIQC(VZV_CHECK.out)
}
//workflow.onComplete {
//    log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
//

