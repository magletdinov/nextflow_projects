#!/usr/bin/env nextflow 
/*
 * pipeline input parameters
 */
//params.run = "07_08_24_rerun_02_08_24"
params.run = "04_10_24_nextseq"
params.shared = "/export/home/public/agletdinov_shared"
params.results_project = "/export/home/agletdinov/work/nextflow_projects/total_seq"
params.reads = "${params.results_project}/fastq/${params.run}/*R{1,2}*.fastq.gz"
//params.reads = "${params.results_project}/fastq/02_08_24/*R{1,2}*.fastq.gz"
params.singleEnd = false

//params.reads = "${params.results_project}/fastq/${params.run}/k18*R{1,2}*.fastq.gz"
params.adapters = "${params.shared}/adapters/adapters.fasta"
params.metaphlandb = "/export/home/public/agletdinov_shared/metaphlandb"
def bracken_settings_dict = [
    'S': 'species',
    'G': 'genus',
]
params.bracken_settings_dict = bracken_settings_dict
params.bracken_settings = ['S', 'G']
params.krakentools_flag = false
params.extract_taxid = false
//def taxid_dict = [
//    '3050337': ["k10_bird_S5"],
//    '694014':   ["k18_bird_S13", "k16_bird_S11", "k24_bird_S19"],
//]

def taxid_dict = [
    '1286':   ["PL-02-08_S1_L001", "LR-02-08_S3_L001", "LL-02-08_S2_L001"], 
    '649161': ["PL-02-08_S1_L001", "LR-02-08_S3_L001", "LL-02-08_S2_L001"],
    '28901':  ["PL-02-08_S1_L001", "LR-02-08_S3_L001", "LL-02-08_S2_L001"],
    '1767':   ["PL-02-08_S1_L001", "LR-02-08_S3_L001"],
    '2267275':["PL-02-08_S1_L001", "LR-02-08_S3_L001", "LL-02-08_S2_L001"],
    '1764':   ["PL-02-08_S1_L001", "LR-02-08_S3_L001", "LL-02-08_S2_L001"],
    
    '6035':   ["PL-02-08_S1_L001"],
    
    '3050299':["LR-02-08_S3_L001", "LL-02-08_S2_L001"],
    '2163996':["LR-02-08_S3_L001", "LL-02-08_S2_L001"],
    '68416'  :["LR-02-08_S3_L001", "LL-02-08_S2_L001"],
] 

params.taxid_dict = taxid_dict

def taxid_dict_v2 = [
    '1286':   ["PL-02-08_S1_L001-1286-megahit-1286", "LR-02-08_S3_L001-1286-megahit-1286", "LL-02-08_S2_L001-1286-megahit-1286"], 
    '649161': ["PL-02-08_S1_L001-649161-megahit-649161", "LR-02-08_S3_L001-649161-megahit-649161", "LL-02-08_S2_L001-649161-megahit-649161"],
    '28901':  ["PL-02-08_S1_L001-28901-megahit-28901", "LR-02-08_S3_L001-28901-megahit-28901", "LL-02-08_S2_L001-28901-megahit-28901"],
    '1767':   ["PL-02-08_S1_L001-1767-megahit-1767", "LR-02-1767"],
    '2267275':["PL-02-08_S1_L001-2267275-megahit-2267275", "LR-02-08_S3_L001-2267275-megahit-2267275", "LL-02-08_S2_L001-2267275-megahit-2267275"],
    '1764':   ["PL-02-08_S1_L001-1767-megahit-1767", "LR-02-08_S3_L001-1767-megahit-1767", "LL-02-08_S2_L001-1767-megahit-1767"],
    
    '6035':   ["PL-02-08_S1_L001-6035-megahit-6035"],
    
    '3050299':["LR-02-08_S3_L001-3050299-megahit-3050299", "LL-02-08_S2_L001-3050299-megahit-3050299"],
    '2163996':["LR-02-08_S3_L001-2163996-megahit-2163996", "LL-02-08_S2_L001-2163996-megahit-2163996"],
    '68416'  :["LR-02-08_S3_L001-68416-megahit-68416", "LL-02-08_S2_L001-68416-megahit-68416"],
] 

params.taxid_dict_v2 = taxid_dict_v2



//def taxid_dict = [
//    '3049954': ["K_S2_L001", "L_S3_L001"],
//    '40324'  : ["L_S3_L001"]
//]
//params.taxid_dict = taxid_dict

//taxid_dict_v2 = params.taxid_dict.collectEntries{ taxid, names ->
//  [taxid, names.collect{ name -> "${name}-${taxid}-megahit-${taxid}" }]
//}
//params.taxid_dict_v2 = taxid_dict_v2

//def taxid_dict_v2 = [
//    '3050337': ['k10_bird_S5-3050337-megahit-3050337'],
//    '694014':   ['k18_bird_S13-694014-megahit-694014', 'k16_bird_S11-694014-megahit-694014', 'k24_bird_S19-694014-megahit-694014'],
//]
//params.taxid_dict_v2 = taxid_dict_v2

//params.taxid = ['1286', '649161', '28901', '1767', '2267275', '1764', '6035', '3050299', '2163996', '68416']
//params.taxid = ['3049954', '40324']
params.taxid = 1491

params.methods = ["4"]
params.bact_genome_dir = "/export/home/public/agletdinov_shared/genomes/bacterias"
params.vir_genome_dir = "/export/home/public/agletdinov_shared/genomes/vir"
params.genomes = ["ticks"]
params.genome = "${params.shared}/genomes/ticks2/ticks.fasta"
params.outdir = "${params.results_project}/results/${params.run}"
//params.outdir = "${params.shared}/nextflow_projects/total_seq/results/${params.run}"
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
    'ticks': "${params.shared}/genomes/ticks2/ticks.fasta",
]
params.genome_dict_taxid = genome_dict_taxid

params.taxid_list = file('/export/home/public/agletdinov_shared/kraken2db/k2_eupathdb48/list_of_eupath_taxid.tsv')
    .text
    .split('\t')
    .findAll { it } // Убираем пустые строки
    .join(',')

params.blastnDB = "/export/home/public/tools/database/nt"
params.db = "/export/home/public/tools/database/nt"
params.to_nodes = "/export/home/public/agletdinov_shared/ncbi_taxonomy/update/nodes.dmp"
params.to_names = "/export/home/public/agletdinov_shared/ncbi_taxonomy/update/names.dmp"
params.bowtie2db = "/export/home/public/agletdinov_shared/bowtie2db/"
def bowtie2_index = [
    '10_tick_S14': 'GRCh38_noalt_as',
    '11_tick_S15': 'GRCh38_noalt_as',
    '12_tick_S16': 'GRCh38_noalt_as',
    '13_tick_S17': 'GRCh38_noalt_as',
    '14_tick_S18': 'GRCh38_noalt_as',
    '15_tick_S19': 'GRCh38_noalt_as',
    '16_tick_S20': 'GRCh38_noalt_as',
    '17_tick_S21': 'GRCh38_noalt_as',
    '18_tick_S22': 'GRCh38_noalt_as',
    '19_tick_S23': 'GRCh38_noalt_as',
    '1_tick_S8': 'GRCh38_noalt_as',
    '20_tick_S24': 'GRCh38_noalt_as',
    '21_tick_S25': 'GRCh38_noalt_as',
    '22_tick_S26': 'GRCh38_noalt_as',
    '23_tick_S27': 'GRCh38_noalt_as',
    '24_tick_S28': 'GRCh38_noalt_as',
    '25_tick_S29': 'GRCh38_noalt_as',
    '26_tick_S30': 'GRCh38_noalt_as',
    '2_tick_S9': 'GRCh38_noalt_as',
    '6_tick_S10': 'GRCh38_noalt_as',
    '7_tick_S11': 'GRCh38_noalt_as',
    '8_tick_S12': 'GRCh38_noalt_as',
    '9_tick_S13': 'GRCh38_noalt_as',
    'C-_1-26_S31': 'GRCh38_noalt_as'
]
params.bowtie2_index = bowtie2_index
params.chunkSize = 100
params.blastn_report = "${params.outdir}/blastn"
log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    reads          : ${params.reads}
    adapters       : ${params.adapters}
    kraken2db      : ${params.kraken2db}
    outdir         : ${params.outdir}
    maxForks       : ${params.maxForks}
    blastnDB       : ${params.blastnDB}
    bowtie2db      : ${params.bowtie2db}
    to_nodes       : ${params.to_nodes}
    to_names       : ${params.to_names}
    """
    .stripIndent()

include { TAXONOMY_ANALYSIS } from './modules/total_seq/total_modules.nf'
include { TAXONOMY_ANALYSIS_SIMPLE } from './modules/total_seq/total_modules.nf'
include { TAXONOMY_ANALYSIS_BAD_R2 } from './modules/total_seq/total_modules.nf'
include { TAXONOMY_ANALYSIS_CONTIGS } from './modules/total_seq/total_modules.nf'
include { TAXONOMY_ANALYSIS_READS } from './modules/total_seq/total_modules.nf'
include { TAXONOMY_ANALYSIS_SARS } from './modules/total_seq/total_modules.nf'
include { TAXONOMY_ANALYSIS_READS_EUPATH } from './modules/total_seq/total_modules.nf'
include { TAXONOMY_ANALYSIS_TYSIA } from './modules/total_seq/total_modules.nf'
include { TAXONOMY_ANALYSIS_TICKS } from './modules/total_seq/total_modules.nf'
include { TAXONOMY_ANALYSIS_TYSIA_LITE } from './modules/total_seq/total_modules.nf'


include { MULTIQC } from './modules/multiqc.nf'
include { SENDMAIL; SENDMAIL_PY } from './modules/sendMail.nf'

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
    taxid_dict_v2 = params.taxid_dict_v2
    taxid = params.taxid
    //blastnDB = file(params.blastnDB)
    db = file( params.db )
    to_nodes = params.to_nodes
    TAXONOMY_ANALYSIS_SIMPLE(read_pairs_ch, methods, bracken_settings, taxid_dict_v2, taxid, db, to_nodes)
    MULTIQC(TAXONOMY_ANALYSIS_SIMPLE.out)
}


workflow taxonomy_analysis_bad_r2{
    Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .set { read_pairs_ch }
    methods = params.methods
    bracken_settings = params.bracken_settings
    taxid = params.taxid
    TAXONOMY_ANALYSIS_BAD_R2(read_pairs_ch, methods, bracken_settings, taxid)
    MULTIQC(TAXONOMY_ANALYSIS_BAD_R2.out)
}

workflow taxonomy_analysis_contigs{
    Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .set { read_pairs_ch }
    methods = params.methods
    bracken_settings = params.bracken_settings
    TAXONOMY_ANALYSIS_CONTIGS(read_pairs_ch, methods, bracken_settings)
    MULTIQC(TAXONOMY_ANALYSIS_CONTIGS.out)
}


workflow taxonomy_analysis_reads{
    Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .set { read_pairs_ch }
    methods = params.methods
    bracken_settings = params.bracken_settings
    db = file( params.db )
    TAXONOMY_ANALYSIS_READS_EUPATH(read_pairs_ch, methods, bracken_settings, db)
    MULTIQC(TAXONOMY_ANALYSIS_READS_EUPATH.out)
}


workflow taxonomy_analysis_reads_eupath{
    Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .set { read_pairs_ch }
    bowtie2db = params.bowtie2db
    taxid_list = params.taxid_list
    bracken_settings = params.bracken_settings
    taxid = "5833"
    TAXONOMY_ANALYSIS_READS_EUPATH(read_pairs_ch, bowtie2db, taxid_list, bracken_settings, taxid)
    MULTIQC(TAXONOMY_ANALYSIS_READS_EUPATH.out)
}

workflow taxonomy_analysis_reads_tysia{
    Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .set { read_pairs_ch }
    bowtie2db = params.bowtie2db
    bowtie2_index = params.bowtie2_index
    taxid = params.taxid
    bracken_settings = params.bracken_settings
    db = file( params.db )
    to_nodes = params.to_nodes
    to_names = params.to_names
    TAXONOMY_ANALYSIS_TYSIA(read_pairs_ch, bowtie2db, taxid, bracken_settings, db, to_nodes, to_names)
    MULTIQC(TAXONOMY_ANALYSIS_TYSIA.out)
    //SENDMAIL_PY(MULTIQC.out)
}

workflow taxonomy_analysis_reads_ticks{
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }
    genome = file(params.genome)
    //genomes = params.genomes
    bracken_settings = params.bracken_settings

    bowtie2db = params.bowtie2db
    bowtie2_index = params.bowtie2_index
    db = file( params.db )
    to_nodes = params.to_nodes
    to_names = params.to_names
    TAXONOMY_ANALYSIS_TICKS(read_pairs_ch, genome, bracken_settings, bowtie2db, db, to_nodes, to_names)
    MULTIQC(TAXONOMY_ANALYSIS_TICKS.out)
}

workflow taxonomy_analysis_reads_tysia_lite{
    Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .set { read_pairs_ch }
    bowtie2db = params.bowtie2db
    bowtie2_index = params.bowtie2_index
    taxid = params.taxid
    bracken_settings = params.bracken_settings
    TAXONOMY_ANALYSIS_TYSIA_LITE(read_pairs_ch, bowtie2db, taxid, bracken_settings)
    MULTIQC(TAXONOMY_ANALYSIS_TYSIA_LITE.out)
    //SENDMAIL_PY(MULTIQC.out)
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