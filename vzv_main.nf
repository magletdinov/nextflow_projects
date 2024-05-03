#!/usr/bin/env nextflow 
/*
 * pipeline input parameters
 */
params.shared = "/export/home/public/agletdinov_shared"
params.results_project = "/export/home/agletdinov/work/nextflow_projects/VZV"
params.reads = "${params.results_project}/fastq/25_04_24/*R{1,2}*.fastq.gz"
params.adapters = "${params.shared}/adapters/adapters.fasta"
params.genome = "${params.shared}/genomes/vzv/NC_001348.1.fasta"
params.primer_a = "${params.shared}/primers/vzv/aADAPTERX_v2.fasta"
params.primer_g = "${params.shared}/primers/vzv/gXADAPTER_v2.fasta"
params.key_areas = "${params.shared}/genomes/vzv/vzv_amplicons.fasta"
params.script1 = "${projectDir}/bin/rename_fasta_id.py"
params.script2 = "/export/home/agletdinov/work/git_projects/gitlab/vzv-classifier/vzv_analysis.py"
params.script3 = "${projectDir}/bin/vzv_run.py"
params.script4 = "${projectDir}/bin/clade_aggregate.py"
params.kraken2db = "/export/home/public/agletdinov_shared/kraken2db/minikraken2_v2_8GB_201904_UPDATE"
params.outdir = "${params.results_project}/results/25_04_24"
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
    script1        : ${params.script1}
    script2        : ${params.script2}
    script3        : ${params.script3}
    script4        : ${params.script4}
    kraken2db      : ${params.kraken2db}
    outdir         : ${params.outdir}
    maxForks       : ${params.maxForks}
    """
    .stripIndent()

include { VZV } from './modules/vzv/vzv_modules.nf'
include { VZV_CHECK } from './modules/vzv/vzv_modules.nf'
include { MULTIQC } from './modules/multiqc.nf'
include { AGG_CLADES } from './modules/vzv/vzv_clades.nf'

workflow vzv_full_run{
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }
    key_areas = file(params.key_areas)
    script1 = params.script1
    script2 = params.script2
    script3 = params.script3
    script4 = params.script4
    VZV(read_pairs_ch, key_areas, script1, script2, script3)
    AGG_CLADES(script4, VZV.out)
    MULTIQC(VZV.out | concat(AGG_CLADES.out) | collect)
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

