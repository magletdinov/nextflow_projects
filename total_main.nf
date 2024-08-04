#!/usr/bin/env nextflow 
/*
 * pipeline input parameters
 */
params.run = "20_07_24"
params.shared = "/export/home/public/agletdinov_shared"
params.results_project = "/export/home/agletdinov/work/nextflow_projects/total_seq"
params.reads = "${params.results_project}/fastq/${params.run}/*R{1,2}*.fastq.gz"
//params.reads = "${params.results_project}/fastq/${params.run}/k18*R{1,2}*.fastq.gz"
params.adapters = "${params.shared}/adapters/adapters.fasta"
params.metaphlandb = "/export/home/public/agletdinov_shared/metaphlandb"
def bracken_settings_dict = [
    'S': 'species',
    'G': 'genus',
]
params.bracken_settings_dict = bracken_settings_dict
params.bracken_settings = ['S', 'G']
params.krakentools_flag = true
def taxid_dict = [
    '3050337': ["k10_bird_S5"],
    '694014':   ["k18_bird_S13", "k16_bird_S11", "k24_bird_S19"],
]
params.taxid_dict = taxid_dict
//taxid_dict_v2 = params.taxid_dict.collectEntries{ taxid, names ->
//  [taxid, names.collect{ name -> "${name}-${taxid}-megahit" }]
//}
//params.taxid_dict_v2 = taxid_dict_v2

def taxid_dict_v2 = [
    '3050337': ["k10_bird_S5-3050337-megahit-3050337"],
    '694014':   ["k18_bird_S13-694014-megahit-694014", "k16_bird_S11-694014-megahit-694014", "k24_bird_S19-694014-megahit-694014"],
]
params.taxid_dict_v2 = taxid_dict_v2

params.taxid = ['3050337', '694014']

params.methods = ["4"]
params.bact_genome_dir = "/export/home/public/agletdinov_shared/genomes/bacterias"
params.vir_genome_dir = "/export/home/public/agletdinov_shared/genomes/vir"
params.genomes = ["cp", "sp", "va", "ec"]
params.genome = "${params.shared}/genomes/sars_cov_2/NC_045512.2.fasta"
params.outdir = "${params.results_project}/results/${params.run}"
params.bwa_index = "${params.outdir}/bwa_index"
//params.maxForks = 50  // Задайте необходимое максимальное число процессов

def genome_dict = [
    'cp': "${params.bact_genome_dir}/chlamydia_pneumoniae/chlamydia_pneumoniae.fasta",
    'sp': "${params.bact_genome_dir}/streptococcus_pneumoniae/streptococcus_pneumoniae.fasta",
    'va': "${params.bact_genome_dir}/veillonella_atypica/veillonella_atypica.fasta",
    'ec': "${params.bact_genome_dir}/escherichia_coli/escherichia_coli.fasta"
]
params.genome_dict = genome_dict
def genome_dict_taxid = [
    '3050337': "${params.vir_genome_dir}/papiine_betaherpesvirus_3.fasta",
    '694014': "${params.vir_genome_dir}/infectious_bronchitis_virus.fasta",
]
params.genome_dict_taxid = genome_dict_taxid
log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    reads          : ${params.reads}
    adapters       : ${params.adapters}
    kraken2db      : ${params.kraken2db}
    outdir         : ${params.outdir}
    maxForks       : ${params.maxForks}
    taxid_dict     :${params.taxid_dict}
    taxid_dict_v2  :${params.taxid_dict_v2}
    """
    .stripIndent()

include { TAXONOMY_ANALYSIS } from './modules/total_seq/total_modules.nf'
include { TAXONOMY_ANALYSIS_SIMPLE } from './modules/total_seq/total_modules.nf'
include { TAXONOMY_ANALYSIS_SARS } from './modules/total_seq/total_modules.nf'
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

workflow taxonomy_analysis_simple{
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }
    methods = params.methods
    bracken_settings = params.bracken_settings
    taxid = params.taxid
    TAXONOMY_ANALYSIS_SIMPLE(read_pairs_ch, methods, bracken_settings, taxid)
    MULTIQC(TAXONOMY_ANALYSIS_SIMPLE.out)
}

workflow taxonomy_analysis_sars{
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }
    key_areas = file(params.genome)
    methods = params.methods
    bracken_settings = params.bracken_settings
    TAXONOMY_ANALYSIS_SARS(read_pairs_ch, key_areas, methods, bracken_settings)
    MULTIQC(TAXONOMY_ANALYSIS_SARS.out)
}
//workflow.onComplete {
//    log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
//