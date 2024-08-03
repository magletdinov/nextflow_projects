include { EXTRACT_KRAKEN_READS_TEST } from 'modules/modules.nf'

params.run = "20_07_24"
params.results_project = "/export/home/agletdinov/work/nextflow_projects/total_seq"
params.outdir = "${params.results_project}/results/${params.run}"
params.reads = "${params.results_project}/results/${params.run}/cutadapt_trim/*R{1,2}*.fastq.gz"
params.kraken2_res = "${params.results_project}/results/${params.run}"

workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }
    read_pairs_ch
    EXTRACT_KRAKEN_READS_TEST(read_pairs_ch)


}