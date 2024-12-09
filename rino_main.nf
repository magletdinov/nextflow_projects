#!/usr/bin/env nextflow 
/*
 * pipeline input parameters
 */
params.run = "02_11_24"
params.shared = "/export/home/public/agletdinov_shared"
params.results_project = "/export/home/agletdinov/work/nextflow_projects/rino"
params.reads = "${params.results_project}/fastq/${params.run}/*R{1,2}*.fastq.gz"
params.singleEnd = false

params.adapters = "${params.shared}/adapters/adapters.fasta"
params.primer_a = "${params.shared}/primers/vzv/aADAPTERX_v2.fasta"
params.primer_g = "${params.shared}/primers/vzv/gXADAPTER_v2.fasta"

params.metaphlandb = "/export/home/public/agletdinov_shared/metaphlandb"
def bracken_settings_dict = [
    'S': 'species',
    'G': 'genus',
]
params.bracken_settings_dict = bracken_settings_dict
params.bracken_settings = ['S', 'G']
params.krakentools_flag = false
params.extract_taxid = false

params.taxid = 1491

params.methods = ["4"]
params.outdir = "${params.results_project}/results/${params.run}"
params.bwa_index = "${params.outdir}/bwa_index"
//params.maxForks = 50  // Задайте необходимое максимальное число процессов

params.blastnDB = "/export/home/public/tools/database/nt"
params.db = "/export/home/public/tools/database/nt"
params.to_nodes = "/export/home/public/agletdinov_shared/ncbi_taxonomy/update/nodes.dmp"
params.to_names = "/export/home/public/agletdinov_shared/ncbi_taxonomy/update/names.dmp"
params.bowtie2db = "/export/home/public/agletdinov_shared/bowtie2db/"
def bowtie2_index = [
    '37-BRICS7_S10_L001': 'GRCh38_noalt_as',
    '38-BRICS8_S11_L001': 'GRCh38_noalt_as',
    '39-BRICS9_S12_L001': 'GRCh38_noalt_as',
    '40-BRICS10_S13_L001': 'GRCh38_noalt_as',
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

include { RINO } from './modules/rino/rino_modules.nf'


include { MULTIQC } from './modules/multiqc.nf'
include { SENDMAIL; SENDMAIL_PY } from './modules/sendMail.nf'

workflow rino{
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
    RINO(read_pairs_ch, bowtie2db, db, to_nodes, to_names)
    MULTIQC(RINO.out)
    //SENDMAIL_PY(MULTIQC.out)
}