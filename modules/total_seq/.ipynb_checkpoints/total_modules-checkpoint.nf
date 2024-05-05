include { FASTQC as FQ1 ; FASTQC as FQ2 } from '../modules.nf'
include { TRIM_ADAPT; TRIM_4_NUCL; KRAKEN2; BRACKEN } from '../modules.nf'


workflow TAXONOMY_ANALYSIS {
  take:
    read_pairs_ch
  main:
    //FASTQC_1 = FASTQC(read_pairs_ch)
    FQ1(read_pairs_ch)
    TRIM_ADAPT(read_pairs_ch)
    TRIM_4_NUCL(TRIM_ADAPT.out)
    //FASTQC_2 = FASTQC(TRIM_PRIMERS.out)
    //FQ2(TRIM_PRIMERS.out)
    KRAKEN2(TRIM_4_NUCL.out)
    BRACKEN(KRAKEN2.out.id_report)
  emit: 
     FQ1.out | concat(KRAKEN2.out.report) | concat(BRACKEN.out) | collect
     //FASTQC.out | concat(KRAKEN2.out) | collect
}