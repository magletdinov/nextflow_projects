include { FASTQC as FQ1 } from '../modules.nf'
include { DORADO_GPU } from '../modules.nf'

workflow NANOPORE {
  take:
    read_pairs_ch
  main:
    DORADO_GPU(read_pairs_ch)
    FQ1(DORADO_GPU.out)
  emit: 
     FQ1.out
     //FQ1.out | concat(KRAKEN2_1.out.report) | collect
}