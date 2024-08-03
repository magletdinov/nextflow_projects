include { FASTQC as FQ1 ; FASTQC as FQ2 } from '../modules.nf'
include { METAPHLAN as METAPHLAN1 ; METAPHLAN as METAPHLAN2 } from '../modules.nf'
include { TRIM_ADAPT; TRIM_4_NUCL; KRAKEN2; BRACKEN; BRACKEN_EACH; EXTRACT_KRAKEN_READS; BWA_INDEX; BWA_INDEX_FULL_GENOME; BWA_MEM_BAM_SORT; BWA_MEM_BAM_SORT_FULL_GENOME; SAMTOOLS_INDEX; SAMTOOLS_CONSENSUS_LITE; SAMTOOLS_STATS; MEGAHIT; UNICYCLER; DIAMOND } from '../modules.nf'


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
    //MEGAHIT(TRIM_4_NUCL.out)
    //DIAMOND(MEGAHIT.out.id_contigs)
    //METAPHLAN1(MEGAHIT.out.id_contigs)
    //METAPHLAN_AGG(METAPHLAN.out.collect())
    
  emit: 
     FQ1.out | concat(KRAKEN2.out.report) | concat(BRACKEN_EACH.out) | collect
     //FASTQC.out | concat(KRAKEN2.out) | collect
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