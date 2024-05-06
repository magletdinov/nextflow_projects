include { FASTQC as FQ1 ; FASTQC as FQ2 } from '../modules.nf'
include { TRIM_ADAPT; TRIM_4_NUCL; KRAKEN2; BRACKEN; BWA_INDEX_FULL_GENOME; BWA_MEM_BAM_SORT_FULL_GENOME; SAMTOOLS_INDEX; SAMTOOLS_CONSENSUS_LITE; SAMTOOLS_STATS; MEGAHIT; METAPHLAN } from '../modules.nf'


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
    TRIM_4_NUCL(TRIM_ADAPT.out, methods)
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
    METAPHLAN(MEGAHIT.out.id_contigs)
    //METAPHLAN_AGG(METAPHLAN.out.collect())
    
  emit: 
     FQ1.out | concat(KRAKEN2.out.report) | concat(BRACKEN.out) | concat(SAMTOOLS_STATS.out) | concat(METAPHLAN.out.report) | collect
     //FASTQC.out | concat(KRAKEN2.out) | collect
}