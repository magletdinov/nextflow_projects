params.run = "13_06_26-and-28_04_25-gonor"
params.shared = "/export/home/public/agletdinov_shared"
params.results_project = "/export/home/agletdinov/work/nextflow_projects/bact_cultures"
params.reads = "${params.results_project}/fastq/${params.run}/*R{1,2}*.fastq.gz"
params.outdir = "${params.results_project}/results/${params.run}/reads_preprocessing"

params.adapters = "${params.shared}/adapters"
params.taxid_genus = 482 //Извлекаем не только Neisseria gonorrhoeae, но и всех потомков рода ${БАКТЕРИЯ}. В идеале после сборки проверить, что собрались только БАКТЕРИЯ (taxid)
params.taxid_species = 485
//params.kraken2db = "/dev/shm/nt_2023"
//params.kraken2db = "/dev/shm/standardPFP"

params.ref_genome = "${params.shared}/genomes/bacterias/neisseria_gonorrhoeae/neisseria_gonorrhoeae.fasta"
params.ref_annotation = "${params.shared}/annotations/bacterias/neisseria_gonorrhoeae/neisseria_gonorrhoeae.gtf"

params.singleEnd_TRIM_4_NUCL = false
params.singleEnd_KRAKEN2 = false
params.singleEnd_EXTRACT_KRAKEN_READS_TAXID = false
params.singleEnd_EXTRACT_KRAKEN_READS_TAXID_FASTA = true

include { TRIMMOMATIC; TRIM_4_NUCL; EXTRACT_KRAKEN_READS_TAXID; isolateSPAdes; EXTRACT_KRAKEN_READS_TAXID_FASTA; QUAST_SINGLE_INPUT; } from '../modules/modules.nf'
include { KRAKEN2 as KRAKEN2_1; KRAKEN2_FASTA as KRAKEN2_2; } from '../modules/modules.nf'
include { FASTQC as FASTQC_1; } from '../modules/modules.nf'
include { FASTQC as FASTQC_2; } from '../modules/modules.nf'

include { MULTIQC } from '../modules/multiqc.nf'

workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }
    dir_name_0 = "raw_data"
    FASTQC_1(read_pairs_ch, dir_name_0)
    TRIMMOMATIC(read_pairs_ch)
    trim_4_nucl_files = TRIM_4_NUCL(TRIMMOMATIC.out.sample_id__p_fastq)
    dir_name_1 = "after_filtering"
    FASTQC_2(trim_4_nucl_files.fastq, dir_name_1)
    dir_name = "before_extract"
    KRAKEN2_1(trim_4_nucl_files.sample_id__fastq, dir_name)
    EXTRACT_KRAKEN_READS_TAXID(trim_4_nucl_files.sample_id__fastq.join(KRAKEN2_1.out.sample_id__output).join(KRAKEN2_1.out.sample_id__report), params.taxid_genus)
    isolateSPAdes(EXTRACT_KRAKEN_READS_TAXID.out.sample_id__fastq)
    dir_name = "after_extract"
    KRAKEN2_2(isolateSPAdes.out.sample_id__scaffolds, dir_name)
    EXTRACT_KRAKEN_READS_TAXID_FASTA(isolateSPAdes.out.sample_id__scaffolds.join(KRAKEN2_2.out.sample_id__output).join(KRAKEN2_2.out.sample_id__report), params.taxid_species)
    ref_genome = file(params.ref_genome)
    ref_annotation = file(params.ref_annotation)
    QUAST_SINGLE_INPUT(EXTRACT_KRAKEN_READS_TAXID_FASTA.out.sample_id__fasta.join(EXTRACT_KRAKEN_READS_TAXID.out.sample_id__fastq), ref_genome, ref_annotation)
    MULTIQC(KRAKEN2_1.out.both_report_output | concat(KRAKEN2_2.out.both_report_output) | concat(QUAST_SINGLE_INPUT.out.report) | concat(FASTQC_1.out) | concat(FASTQC_2.out) | collect)
}