include { FASTQC as FQ1; FASTQC as FQ2 } from '../modules.nf'
include { DORADO_GPU; DORADO_CORRECTION_GPU } from '../modules.nf'

workflow NANOPORE {
  take:
    read_pairs_ch
  main:
    DORADO_GPU(read_pairs_ch).view()
    FQ1(DORADO_GPU.out)
    DORADO_CORRECTION_GPU(DORADO_GPU.out).view()
    FQ2(DORADO_CORRECTION_GPU.out)
  emit: 
     FQ1.out | concat(FQ2.out) | collect
     //FQ1.out | concat(KRAKEN2_1.out.report) | collect
}