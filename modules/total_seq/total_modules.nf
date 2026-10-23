include { FASTQC as FQ1 ; FASTQC as FQ2 } from '../modules.nf'
include { METAPHLAN as METAPHLAN1 ; METAPHLAN as METAPHLAN2 } from '../modules.nf'
include { TRIM_ADAPT; TRIM_4_NUCL; KRAKEN2; BRACKEN; BRACKEN_EACH; EXTRACT_KRAKEN_READS; EXTRACT_KRAKEN_READS_FASTA; EXTRACT_KRAKEN_READS_TAXID_LIST; EXTRACT_KRAKEN_READS_TAXID; EXTRACT_KRAKEN_READS_FROM_FASTA_TO_FASTA; BWA_INDEX; BWA_INDEX_FULL_GENOME; BWA_MEM_BAM_SORT; BWA_MEM_BAM_SORT_FULL_GENOME; BWA_INDEX_MEM_BAM_SORT; SAMTOOLS_INDEX; SAMTOOLS_CONSENSUS_LITE; SAMTOOLS_STATS; SAMTOOLS_EXTRACT_FASTQ; MEGAHIT; MEGAHIT_MS; QUAST; UNICYCLER; BLASTN; blastn; blastn_parse_1_hit; blastn_parse_50_hits; DIAMOND; REMOVE_HOST; REMOVE_HUMAN; metaSPAdes; shuffling_fasta; QUAST_NO_REFERENCE } from '../modules.nf'
//include { KRAKEN2 as KRAKEN2_1 ; KRAKEN2 as KRAKEN2_2 } from '../modules.nf'
include { KRAKEN2 as KRAKEN2_1 ; KRAKEN2 as KRAKEN2_2 } from '../modules.nf'
include { KRAKEN2_FASTA as KRAKEN2_FASTA_1 ; KRAKEN2_FASTA as KRAKEN2_FASTA_2 } from '../modules.nf'
include { QUAST_NO_REFERENCE as QUAST_NO_REFERENCE_PRIMARY; QUAST_NO_REFERENCE as QUAST_NO_REFERENCE_SECONDARY; QUAST_NO_REFERENCE as QUAST_NO_REFERENCE_SECONDARY_FILT_LENGTH } from '../modules.nf'

include { FILTER_FASTA_BY_LENGTH as FILTER_FASTA_BY_LENGTH_1 ; FILTER_FASTA_BY_LENGTH as FILTER_FASTA_BY_LENGTH_2 } from '../modules.nf'
include { blastn as blastn_1 ; blastn as blastn_2 } from '../modules.nf'
include { blastn_parse_50_hits as blastn_parse_50_hits_1 ; blastn_parse_50_hits as blastn_parse_50_hits_2 } from '../modules.nf'
workflow TAXONOMY_ANALYSIS_ALHOVSKIJ {
  take:
    read_pairs_ch
    primary_assembles_ch
    bracken_settings
    db
    to_nodes
    to_names
  main:
    read_pairs_ch//.view()
    FQ1(read_pairs_ch)
    dir_name = "raw_data"
    KRAKEN2_1(read_pairs_ch, dir_name)
    dir1_name = "primary_assembles"
    FILTER_FASTA_BY_LENGTH_1(primary_assembles_ch, '1000', dir1_name)
    blastn_1(FILTER_FASTA_BY_LENGTH_1.out.sample_id__filtered_fasta, db, dir1_name)
    blastn_parse_50_hits_1(blastn_1.out.sample_id__report.join(FILTER_FASTA_BY_LENGTH_1.out.sample_id__filtered_fasta), to_nodes, to_names, dir1_name)
    
    KRAKEN2_FASTA_1(primary_assembles_ch, dir_name)
    dir8_name = "primary_assembles_before_length_filtering"
    KRAKEN2_FASTA_2(primary_assembles_ch, dir8_name)
    dir2_name = "viral_and_unclassified_from_primary_assembles"
    EXTRACT_KRAKEN_READS_FROM_FASTA_TO_FASTA(primary_assembles_ch.join(KRAKEN2_FASTA_1.out.sample_id__output).join(KRAKEN2_FASTA_1.out.sample_id__report), dir2_name)
    BWA_INDEX(EXTRACT_KRAKEN_READS_FROM_FASTA_TO_FASTA.out.sample_id__fasta
        .map { i, j -> tuple(j, i)})//.view()
    combine_ch = read_pairs_ch.combine(BWA_INDEX.out)//.view()
    //combine_ch.view()
    BWA_MEM_BAM_SORT(combine_ch)
    SAMTOOLS_EXTRACT_FASTQ(BWA_MEM_BAM_SORT.out.sample_id__bam)
    
    origin_ch = SAMTOOLS_EXTRACT_FASTQ.out.sample_id__fastq
    alh_ch = origin_ch.filter { sample_id, files -> sample_id.contains('alhovscij') } //на этих последовательностях запустил CLC
    megahit_ch = origin_ch.filter { sample_id, files -> sample_id.endsWith('-megahit')  }
    megahit_ms_ch = origin_ch.filter { sample_id, files -> sample_id.contains('megahit_metasens') }
    metaspades_ch = origin_ch.filter { sample_id, files -> sample_id.contains('spades_meta') }
    
    //alh_ch.view()
    //megahit_ch.view()
    //megahit_ms_ch.view()
    //metaspades_ch.view()
    
    MEGAHIT(megahit_ch)
    MEGAHIT_MS(megahit_ms_ch)
    
    megahit_res_ch = MEGAHIT.out.sample_id__contigs
    megahit_res_ch = megahit_res_ch.map{ sample_id, file -> tuple(sample_id.split("-megahit_simple")[0], file) }

    megahit_ms_res_ch = MEGAHIT_MS.out.sample_id__contigs
    megahit_ms_res_ch = megahit_ms_res_ch.map{ sample_id, file -> tuple(sample_id.split("-megahit_ms")[0], file) }
   
    Channel
        .fromPath("/export/home/agletdinov/work/nextflow_projects/total_seq/results/novaseq_alhovscij/secondary_assembles/clc/2d-2d_contigs-alhovscij.contigs.fa")
        .map { path -> tuple(path.getName().split(".contigs.fa")[0], path) }
        .set { alh_ch_res_ch }

    dir3_name = "before_length_filtering"
    assembler_reads_res_ch = megahit_res_ch.join(megahit_ch).mix(megahit_ms_res_ch.join(megahit_ms_ch)).mix(alh_ch_res_ch.join(alh_ch))
    QUAST_NO_REFERENCE_SECONDARY(assembler_reads_res_ch, dir3_name)
    
    assembler_res_ch = assembler_reads_res_ch.map { sample_id, fasta, reads -> tuple(sample_id, fasta) }//.view()
    dir6_name = "secondary_assembles"
    FILTER_FASTA_BY_LENGTH_2(assembler_res_ch, '1000', dir6_name)
    //FILTER_FASTA_BY_LENGTH.out.sample_id__filtered_fasta.view()
    
    dir4_name = "after_length_filtering"
    //alh_ch.mix(megahit_ch).mix(megahit_ms_ch).view()
    assembler_filt_reads_res_ch = FILTER_FASTA_BY_LENGTH_2.out.sample_id__filtered_fasta.join(alh_ch.mix(megahit_ch).mix(megahit_ms_ch))
    QUAST_NO_REFERENCE_SECONDARY_FILT_LENGTH(assembler_filt_reads_res_ch, dir4_name)
   
    
    dir7_name = "secondary_assembles"
    blastn_2(FILTER_FASTA_BY_LENGTH_2.out.sample_id__filtered_fasta, db, dir7_name)
    blastn_parse_50_hits_2(blastn_2.out.sample_id__report.join(FILTER_FASTA_BY_LENGTH_2.out.sample_id__filtered_fasta), to_nodes, to_names, dir7_name)
    
  emit: 
     FQ1.out | concat(KRAKEN2_FASTA_1.out.both_report_output) | collect
}


workflow TAXONOMY_ANALYSIS {
  take:
    read_pairs_ch
    methods
    genomes
    bwa_index
  main:
    //FASTQC_1 = FASTQC(read_pairs_ch)
    FQ1(read_pairs_ch)
    TRIM_ADAPT(read_pairs_ch)
    TRIM_4_NUCL(TRIM_ADAPT.out)
    //FASTQC_2 = FASTQC(TRIM_PRIMERS.out)
    //FQ2(TRIM_PRIMERS.out)
    BWA_INDEX_FULL_GENOME(genomes)
    BWA_MEM_BAM_SORT_FULL_GENOME(TRIM_4_NUCL.out, genomes)
    SAMTOOLS_INDEX(BWA_MEM_BAM_SORT_FULL_GENOME.out)
    SAMTOOLS_CONSENSUS_LITE(BWA_MEM_BAM_SORT_FULL_GENOME.out)
    SAMTOOLS_STATS(BWA_MEM_BAM_SORT_FULL_GENOME.out)    
    
    KRAKEN2(TRIM_4_NUCL.out)
    BRACKEN(KRAKEN2.out.id_report)
    MEGAHIT(TRIM_4_NUCL.out)
    UNICYCLER(TRIM_4_NUCL.out)
    METAPHLAN1(MEGAHIT.out.id_contigs)
    METAPHLAN2(UNICYCLER.out.id_contigs)
    //METAPHLAN_AGG(METAPHLAN.out.collect())
    
  emit: 
     FQ1.out | concat(KRAKEN2.out.report) | concat(BRACKEN.out) | concat(SAMTOOLS_STATS.out) | concat(METAPHLAN1.out.report) | concat(METAPHLAN2.out.report) | collect
     //FASTQC.out | concat(KRAKEN2.out) | collect
}

workflow TAXONOMY_ANALYSIS_SIMPLE {
  take:
    read_pairs_ch
    methods
    bracken_settings
    taxid_dict_v2
    taxid
    //blastnDB
    db
    to_nodes
  main:
    //FASTQC_1 = FASTQC(read_pairs_ch)
    FQ1(read_pairs_ch)
    TRIM_ADAPT(read_pairs_ch)
    TRIM_4_NUCL(TRIM_ADAPT.out)
    //FASTQC_2 = FASTQC(TRIM_PRIMERS.out)
    //FQ2(TRIM_PRIMERS.out)
    KRAKEN2(TRIM_4_NUCL.out)
    BRACKEN_EACH(KRAKEN2.out.id_report, bracken_settings)
    EXTRACT_KRAKEN_READS(TRIM_4_NUCL.out.join(KRAKEN2.out.id_output).join(KRAKEN2.out.id_report), taxid)
    MEGAHIT(EXTRACT_KRAKEN_READS.out.id_fasta)
    KRAKEN2_FASTA(MEGAHIT.out.id_contigs)
    EXTRACT_KRAKEN_READS_FASTA(MEGAHIT.out.id_contigs.join(KRAKEN2_FASTA.out.id_output).join(KRAKEN2_FASTA.out.id_report), taxid_dict_v2, taxid)
    blastn(MEGAHIT.out.simple_id_contigs, db)
    blastn_parse_1_hit(blastn.out.id_report, to_nodes).view()
    //BLASTN(MEGAHIT.out.simple_id_contigs, blastnDB)
    
    //def (value1, value2) = '1128-2'.tokenize( '-' )
    //QUAST(EXTRACT_KRAKEN_READS_FASTA.out)
    //DIAMOND(MEGAHIT.out.id_contigs)
    //METAPHLAN1(MEGAHIT.out.id_contigs)
    //METAPHLAN_AGG(METAPHLAN.out.collect())
    
  emit: 
     FQ1.out | concat(KRAKEN2.out.report) | concat(BRACKEN_EACH.out) | collect
     //FASTQC.out | concat(KRAKEN2.out) | collect
}

workflow TAXONOMY_ANALYSIS_CONTIGS {
  take:
    read_pairs_ch
    methods
    bracken_settings
  main:
    FQ1(read_pairs_ch)
    TRIM_ADAPT(read_pairs_ch)
    TRIM_4_NUCL(TRIM_ADAPT.out)
    MEGAHIT(TRIM_4_NUCL.out)
    KRAKEN2_FASTA(MEGAHIT.out.id_contigs)
    BRACKEN_EACH(KRAKEN2_FASTA.out.id_report, bracken_settings)
    
  emit: 
     FQ1.out | concat(KRAKEN2_FASTA.out.report) | concat(BRACKEN_EACH.out) | collect
     //FASTQC.out | concat(KRAKEN2.out) | collect
}

workflow TAXONOMY_ANALYSIS_READS {
  take:
    read_pairs_ch
    methods
    bracken_settings
    db
    to_nodes
    to_names
  main:
    FQ1(read_pairs_ch)
    TRIM_ADAPT(read_pairs_ch)
    TRIM_4_NUCL(TRIM_ADAPT.out)
    KRAKEN2(TRIM_4_NUCL.out)
    REMOVE_HUMAN(TRIM_4_NUCL.out.join(KRAKEN2.out.id_output).join(KRAKEN2.out.id_report))
    BRACKEN_EACH(KRAKEN2.out.id_report, bracken_settings)
    blastn(REMOVE_HUMAN.out.id_onefasta, db)
    blastn_parse_50_hits(blastn.out.sample_id_report, to_nodes, to_names)
  emit: 
     FQ1.out | concat(KRAKEN2.out.report) | concat(BRACKEN_EACH.out) | collect
     //FASTQC.out | concat(KRAKEN2.out) | collect
}


workflow TAXONOMY_ANALYSIS_READS_EUPATH {
  take:
    read_pairs_ch
    bowtie2_db
    taxid_list
    bracken_settings
    taxid
  main:
    FQ1(read_pairs_ch)
    TRIM_ADAPT(read_pairs_ch)
    TRIM_4_NUCL(TRIM_ADAPT.out)
    REMOVE_HOST(TRIM_ADAPT.out, bowtie2_db)
    dir_name = "all_data"
    KRAKEN2_1(REMOVE_HOST.out, dir_name)
    EXTRACT_KRAKEN_READS_TAXID_LIST(REMOVE_HOST.out.join(KRAKEN2_1.out.id_output).join(KRAKEN2_1.out.id_report), taxid_list)
    dir_name = "after_filtration"
    KRAKEN2_2(EXTRACT_KRAKEN_READS_TAXID_LIST.out.id_fastq, dir_name)
    BRACKEN_EACH(KRAKEN2_2.out.id_report, bracken_settings)
    EXTRACT_KRAKEN_READS_TAXID(EXTRACT_KRAKEN_READS_TAXID_LIST.out.id_fastq.join(KRAKEN2_2.out.id_output).join(KRAKEN2_2.out.id_report), taxid)

    
  emit: 
     FQ1.out | concat(KRAKEN2_1.out.report) | concat(BRACKEN_EACH.out) | collect
     //FASTQC.out | concat(KRAKEN2_FASTA_2.out) | collect
}


workflow TAXONOMY_ANALYSIS_TYSIA {
  take:
    read_pairs_ch
    bowtie2_db
    //bowtie2_index
    taxid
    bracken_settings
    db
    to_nodes
    to_names
  main:
    FQ1(read_pairs_ch)
    TRIM_ADAPT(read_pairs_ch)
    TRIM_4_NUCL(TRIM_ADAPT.out)
    REMOVE_HOST(TRIM_4_NUCL.out, bowtie2_db)
    dir_name = "all_data"
    KRAKEN2_1(TRIM_4_NUCL.out, dir_name)
    EXTRACT_KRAKEN_READS_TAXID(TRIM_4_NUCL.out.join(KRAKEN2_1.out.id_output).join(KRAKEN2_1.out.id_report), taxid)
    BRACKEN_EACH(KRAKEN2_1.out.id_report, bracken_settings)
    metaSPAdes(REMOVE_HOST.out)
    shuffling_fasta(metaSPAdes.out.id_scaffolds)
    ch_fasta = shuffling_fasta.out.splitFasta(by: params.chunkSize, file:true)
    //ch_fasta.view()
    ch_hits = blastn(ch_fasta, db).id_report
    //ch_hits.view()
    ch_collected = ch_hits.collectFile(storeDir:params.blastn_report).map{ [it.name, it] }
    //ch_collected.view()
    blastn_parse_50_hits(ch_collected, to_nodes, to_names)
  emit: 
     FQ1.out | concat(KRAKEN2_1.out.report) | concat(blastn_parse_50_hits.out) | collect
     //FQ1.out | concat(KRAKEN2_1.out.report) | collect
}

workflow TAXONOMY_ANALYSIS_TYSIA_LITE {
  take:
    read_pairs_ch
    bowtie2_db
    taxid
    bracken_settings
  main:
    FQ1(read_pairs_ch)
    TRIM_ADAPT(read_pairs_ch)
    TRIM_4_NUCL(TRIM_ADAPT.out)
    REMOVE_HOST(TRIM_4_NUCL.out, bowtie2_db)
    dir_name = "all_data"
    KRAKEN2_1(TRIM_4_NUCL.out, dir_name)
    EXTRACT_KRAKEN_READS_TAXID(TRIM_4_NUCL.out.join(KRAKEN2_1.out.id_output).join(KRAKEN2_1.out.id_report), taxid)
    BRACKEN_EACH(KRAKEN2_1.out.id_report, bracken_settings)
  emit: 
     FQ1.out | concat(KRAKEN2_1.out.report) | concat(BRACKEN_EACH.out) | collect
     //FQ1.out | concat(KRAKEN2_1.out.report) | collect
}

workflow TAXONOMY_ANALYSIS_TICKS {
  take:
    read_pairs_ch
    genome
    //genomes
    bracken_settings
    bowtie2_db
    db
    to_nodes
    to_names
  main:
    FQ1(read_pairs_ch)
    TRIM_ADAPT(read_pairs_ch)
    TRIM_4_NUCL(TRIM_ADAPT.out)
    BWA_INDEX(genome, genome.name)
    //BWA_INDEX_FULL_GENOME(genomes)
    BWA_MEM_BAM_SORT(TRIM_4_NUCL.out, BWA_INDEX.out)
    SAMTOOLS_INDEX(BWA_MEM_BAM_SORT.out)
    SAMTOOLS_STATS(BWA_MEM_BAM_SORT.out)  
    dir_name = "all_data"
    KRAKEN2_1(TRIM_4_NUCL.out, dir_name)
    BRACKEN_EACH(KRAKEN2_1.out.id_report, bracken_settings)

    REMOVE_HOST(TRIM_4_NUCL.out, bowtie2_db)
    metaSPAdes(REMOVE_HOST.out)
    shuffling_fasta(metaSPAdes.out.id_scaffolds)
    ch_fasta = shuffling_fasta.out.splitFasta(by: params.chunkSize, file:true)
    //ch_fasta.view()
    ch_hits = blastn(ch_fasta, db).id_report
    //ch_hits.view()
    ch_collected = ch_hits.collectFile(storeDir:params.blastn_report).map{ [it.name, it] }
    //ch_collected.view()
    blastn_parse_50_hits(ch_collected, to_nodes, to_names)
    
  emit: 
     FQ1.out | concat(SAMTOOLS_STATS.out) | concat(KRAKEN2_1.out.report) | concat(BRACKEN_EACH.out) | concat(blastn_parse_50_hits.out)| collect
     //FQ1.out | concat(KRAKEN2_1.out.report) | collect
}


workflow TAXONOMY_ANALYSIS_SARS {
  take:
    read_pairs_ch
    key_areas
    methods
    bracken_settings
  main:
    FQ1(read_pairs_ch)
    TRIM_ADAPT(read_pairs_ch)
    TRIM_4_NUCL(TRIM_ADAPT.out)
    BWA_INDEX(key_areas, key_areas.name)
    BWA_MEM_BAM_SORT(TRIM_4_NUCL.out, BWA_INDEX.out)
    SAMTOOLS_INDEX(BWA_MEM_BAM_SORT.out)
    SAMTOOLS_CONSENSUS_LITE(BWA_MEM_BAM_SORT.out)
    SAMTOOLS_STATS(BWA_MEM_BAM_SORT.out)    
    KRAKEN2(TRIM_4_NUCL.out)
    BRACKEN_EACH(KRAKEN2.out.id_report, bracken_settings)
    
  emit: 
     FQ1.out | concat(SAMTOOLS_STATS.out) | concat(KRAKEN2.out.report) | concat(BRACKEN_EACH.out) | collect
     //FASTQC.out | concat(KRAKEN2.out) | collect
}

workflow TAXONOMY_ANALYSIS_BAD_R2{
  take:
    read_pairs_ch
    methods
    bracken_settings
    taxid
  main:
    //FASTQC_1 = FASTQC(read_pairs_ch)
    FQ1(read_pairs_ch)
    TRIM_ADAPT(read_pairs_ch)
    TRIM_4_NUCL(TRIM_ADAPT.out)
    //FASTQC_2 = FASTQC(TRIM_PRIMERS.out)
    //FQ2(TRIM_PRIMERS.out)
    KRAKEN2(TRIM_4_NUCL.out)
    BRACKEN_EACH(KRAKEN2.out.id_report, bracken_settings)
    EXTRACT_KRAKEN_READS(TRIM_4_NUCL.out.join(KRAKEN2.out.id_output).join(KRAKEN2.out.id_report), taxid)
  emit: 
     FQ1.out | concat(KRAKEN2.out.report) | concat(BRACKEN_EACH.out) | collect
     //FASTQC.out | concat(KRAKEN2.out) | collect
}