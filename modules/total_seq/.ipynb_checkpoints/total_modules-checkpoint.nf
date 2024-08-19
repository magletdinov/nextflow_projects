include { FASTQC as FQ1 ; FASTQC as FQ2 } from '../modules.nf'
include { METAPHLAN as METAPHLAN1 ; METAPHLAN as METAPHLAN2 } from '../modules.nf'
include { TRIM_ADAPT; TRIM_4_NUCL; KRAKEN2; KRAKEN2_FASTA; BRACKEN; BRACKEN_EACH; EXTRACT_KRAKEN_READS; EXTRACT_KRAKEN_READS_FASTA; EXTRACT_KRAKEN_READS_TAXID_LIST; EXTRACT_KRAKEN_READS_TAXID; BWA_INDEX; BWA_INDEX_FULL_GENOME; BWA_MEM_BAM_SORT; BWA_MEM_BAM_SORT_FULL_GENOME; SAMTOOLS_INDEX; SAMTOOLS_CONSENSUS_LITE; SAMTOOLS_STATS; MEGAHIT; QUAST; UNICYCLER; BLASTN; blastn; blastn_parse_1_hit; blastn_parse_50_hits; DIAMOND; REMOVE_HOST; REMOVE_HUMAN; metaSPAdes; } from '../modules.nf'
//include { KRAKEN2 as KRAKEN2_1 ; KRAKEN2 as KRAKEN2_2 } from '../modules.nf'
include { KRAKEN2 as KRAKEN2_1 ; KRAKEN2 as KRAKEN2_2 } from '../modules.nf'

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
  main:
    FQ1(read_pairs_ch)
    TRIM_ADAPT(read_pairs_ch)
    TRIM_4_NUCL(TRIM_ADAPT.out)
    KRAKEN2(TRIM_4_NUCL.out)
    REMOVE_HUMAN(TRIM_4_NUCL.out.join(KRAKEN2.out.id_output).join(KRAKEN2.out.id_report))
    BRACKEN_EACH(KRAKEN2.out.id_report, bracken_settings)
    blastn(REMOVE_HUMAN.out.id_onefasta, db)
    
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
    bracken_settings
    db
    to_nodes
    to_names
  main:
    FQ1(read_pairs_ch)
    TRIM_ADAPT(read_pairs_ch)
    TRIM_4_NUCL(TRIM_ADAPT.out)
    REMOVE_HOST(TRIM_ADAPT.out, bowtie2_db)
    dir_name = "all_data"
    KRAKEN2_1(REMOVE_HOST.out, dir_name)
    BRACKEN_EACH(KRAKEN2_1.out.id_report, bracken_settings)
    metaSPAdes(REMOVE_HOST.out)
    blastn(metaSPAdes.out.id_scaffolds, db)
    blastn_parse_50_hits(blastn.out.id_report, to_nodes, to_names)

  emit: 
     FQ1.out | concat(KRAKEN2_1.out.report) | concat(BRACKEN_EACH.out) | concat(blastn_parse_50_hits.out) | collect
     //FASTQC.out | concat(KRAKEN2_FASTA_2.out) | collect
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