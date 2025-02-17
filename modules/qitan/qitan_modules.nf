include { FASTQC as FQ1; FASTQC as FQ2 } from '../modules.nf'
include { TRIM_PRIMERS; KRAKEN2; BWA_INDEX; BWA_MEM_BAM_SORT } from '../modules.nf'
workflow QITAN {
  take:
    reads_ch
    genome
  main:
    FQ1(reads_ch)
    //TRIM_PRIMERS(reads_ch)
    //FQ2(TRIM_PRIMERS.out)
    dir_name = "all_data"
    KRAKEN2(reads_ch, dir_name)
    BWA_INDEX(genome, genome.name)
    BWA_MEM_BAM_SORT(reads_ch, BWA_INDEX.out)
    
  emit: 
     //FQ1.out | collect
     FQ1.out | concat(KRAKEN2.out.report) | collect
}