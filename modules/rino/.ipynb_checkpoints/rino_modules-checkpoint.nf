include { FASTQC as FQ1 ; FASTQC as FQ2 } from '../modules.nf'
include { TRIM_ADAPT; TRIM_4_NUCL; TRIM_PRIMERS; REMOVE_HOST; BWA_INDEX; BWA_MEM_BAM_SORT; metaSPAdes; shuffling_fasta; blastn; blastn_parse_50_hits; SAMTOOLS_INDEX; SAMTOOLS_CONSENSUS; SAMTOOLS_STATS } from '../modules.nf'
include { KRAKEN2 as KRAKEN2_1 } from '../modules.nf'

workflow RINO {
  take:
    read_pairs_ch
    bowtie2_db
    db
    to_nodes
    to_names
  main:
    FQ1(read_pairs_ch)
    TRIM_ADAPT(read_pairs_ch)
    TRIM_4_NUCL(TRIM_ADAPT.out)
    TRIM_PRIMERS(TRIM_4_NUCL.out)
    REMOVE_HOST(TRIM_PRIMERS.out, bowtie2_db)
    dir_name = "all_data"
    KRAKEN2_1(REMOVE_HOST.out, dir_name)
    metaSPAdes(REMOVE_HOST.out)
    shuffling_fasta(metaSPAdes.out.id_scaffolds)    
    ch_fasta = shuffling_fasta.out.splitFasta(by: params.chunkSize, file:true)
    ch_hits = blastn(ch_fasta, db).id_report
    ch_collected = ch_hits.collectFile(storeDir:params.blastn_report).map{ [it.name, it] }
    blastn_parse_50_hits(ch_collected, to_nodes, to_names)
  emit: 
     FQ1.out | concat(KRAKEN2_1.out.report) | concat(blastn_parse_50_hits.out) | collect
}