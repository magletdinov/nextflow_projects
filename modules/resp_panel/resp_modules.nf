include { FASTQC; TRIM_ADAPT; TRIM_4_NUCL; TRIM_PRIMERS; KRAKEN2 } from '../modules.nf'


workflow TAXONOMY_ANALYSIS {
  take:
    read_pairs_ch
  main:
    //FASTQC_1 = FASTQC(read_pairs_ch)
    FASTQC(read_pairs_ch)
    TRIM_ADAPT(read_pairs_ch)
    TRIM_4_NUCL(TRIM_ADAPT.out)
    TRIM_PRIMERS(TRIM_4_NUCL.out)
    //FASTQC_2 = FASTQC(TRIM_PRIMERS.out)
    KRAKEN2(TRIM_PRIMERS.out)
  emit: 
     //FASTQC_1.out | concat(FASTQC_2.out) | concat(KRAKEN2.out) | collect
     FASTQC.out | concat(KRAKEN2.out) | collect
}