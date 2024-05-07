#!/usr/bin/env nextflow 
/*
 * pipeline input parameters
 */
params.shared = "/export/home/public/agletdinov_shared"
params.results_project = "/export/home/agletdinov/work/nextflow_projects/total_seq"
params.reads = "${params.results_project}/fastq/02_05_24/*R{1,2}*.fastq.gz"
params.adapters = "${params.shared}/adapters/adapters.fasta"
//params.kraken2db = "/export/home/public/agletdinov_shared/kraken2db/k2_pluspf_16gb"
params.kraken2db = "/export/home/public/agletdinov_shared/kraken2db/minikraken2_v2_8GB_201904_UPDATE"
//params.kraken2db = "/export/home/public/agletdinov_shared/kraken2db/k2_eupathdb48"
params.metaphlandb = "/export/home/public/agletdinov_shared/metaphlandb"
def myDict = [
    '2': 'dva',
    '4': 'cheture',
]
params.myDict = myDict
params.methods = ["4"]
params.bact_genome_dir = "/export/home/public/agletdinov_shared/genomes/bacterias"
params.genomes = ["cp", "sp", "va", "ec"]
params.outdir = "${params.results_project}/results/02_05_24"
params.bwa_index = "${params.outdir}/bwa_index"
//params.maxForks = 50  // Задайте необходимое максимальное число процессов

def genome_dict = [
    'cp': "${params.bact_genome_dir}/chlamydia_pneumoniae/chlamydia_pneumoniae.fasta",
    'sp': "${params.bact_genome_dir}/streptococcus_pneumoniae/streptococcus_pneumoniae.fasta",
    'va': "${params.bact_genome_dir}/veillonella_atypica/veillonella_atypica.fasta",
    'ec': "${params.bact_genome_dir}/escherichia_coli/escherichia_coli.fasta"
]
params.genome_dict = genome_dict

log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    reads          : ${params.reads}
    adapters       : ${params.adapters}
    kraken2db      : ${params.kraken2db}
    outdir         : ${params.outdir}
    maxForks       : ${params.maxForks}
    """
    .stripIndent()

include { TAXONOMY_ANALYSIS } from './modules/total_seq/total_modules.nf'
include { MULTIQC } from './modules/multiqc.nf'

workflow taxonomy_analysis{
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }
    methods = params.methods
    genomes = params.genomes
    bwa_index = params.bwa_index
    TAXONOMY_ANALYSIS(read_pairs_ch, methods, genomes, bwa_index)
    MULTIQC(TAXONOMY_ANALYSIS.out)
}
//workflow.onComplete {
//    log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
//

