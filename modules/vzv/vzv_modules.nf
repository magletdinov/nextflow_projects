include { FASTQC; TRIM_ADAPT; TRIM_4_NUCL; TRIM_PRIMERS; BWA_INDEX; BWA_MEM_BAM_SORT; SAMTOOLS_INDEX; SAMTOOLS_CONSENSUS; SAMTOOLS_STATS; RENAME_FASTA_ID; KRAKEN2; IDENTIFY_CLADE } from '../modules.nf'


workflow VZV {
  take:
    read_pairs_ch
    key_areas
    script1
    script2
    script3
    //methods
  main:
    FASTQC(read_pairs_ch)
    TRIM_ADAPT(read_pairs_ch)
    TRIM_4_NUCL(TRIM_ADAPT.out)
    TRIM_PRIMERS(TRIM_4_NUCL.out)
    KRAKEN2(TRIM_PRIMERS.out)
    BWA_INDEX(key_areas, key_areas.name)
    BWA_MEM_BAM_SORT(TRIM_PRIMERS.out, BWA_INDEX.out)
    SAMTOOLS_INDEX(BWA_MEM_BAM_SORT.out)
    SAMTOOLS_CONSENSUS(BWA_MEM_BAM_SORT.out)
    SAMTOOLS_STATS(BWA_MEM_BAM_SORT.out)
    RENAME_FASTA_ID(SAMTOOLS_CONSENSUS.out, script1)
    IDENTIFY_CLADE(RENAME_FASTA_ID.out, script2, script3)
  emit: 
     //BWA_MEM_BAM_SORT.out | concat(SAMTOOLS_STATS.out) | concat(KRAKEN2.out) | collect
     FASTQC.out | concat(SAMTOOLS_STATS.out) | concat(KRAKEN2.out.report) | concat(IDENTIFY_CLADE.out) | collect
     //SAMTOOLS_STATS.out | concat(KRAKEN2.out.report) | collect
}

workflow VZV_CHECK {
  take:
    read_pairs_ch
    key_areas
    script
  main:
    FASTQC(read_pairs_ch)
    TRIM_ADAPT(read_pairs_ch)
    TRIM_4_NUCL(TRIM_ADAPT.out)
    TRIM_PRIMERS(TRIM_4_NUCL.out)
    KRAKEN2(TRIM_PRIMERS.out)
    BWA_INDEX(key_areas, key_areas.name)
    BWA_MEM_BAM_SORT(TRIM_PRIMERS.out, BWA_INDEX.out)
    SAMTOOLS_INDEX(BWA_MEM_BAM_SORT.out)
    SAMTOOLS_CONSENSUS(BWA_MEM_BAM_SORT.out)
    SAMTOOLS_STATS(BWA_MEM_BAM_SORT.out)
  emit: 
     //BWA_MEM_BAM_SORT.out | concat(SAMTOOLS_STATS.out) | concat(KRAKEN2.out) | collect
     FASTQC.out | concat(SAMTOOLS_STATS.out) | concat(KRAKEN2.out) | collect
}