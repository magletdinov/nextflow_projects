include { FASTQC as FQ1 ; FASTQC as FQ2 } from '../modules.nf'
include { BWA_INDEX as BWA_INDEX1 ; BWA_INDEX as BWA_INDEX2 } from '../modules.nf'
include { TRIM_ADAPT; TRIM_4_NUCL; TRIM_PRIMERS; KRAKEN2; BRACKEN_UPDATED; BWA_MEM_BAM_SORT; SAMTOOLS_STATS } from '../modules.nf'

workflow TAXONOMY_ANALYSIS {
  take:
    read_pairs_ch
  main:
    //FASTQC_1 = FASTQC(read_pairs_ch)
    dir_name_0 = "raw_data"
    FQ1(read_pairs_ch, dir_name_0)
    TRIM_ADAPT(read_pairs_ch)
    TRIM_4_NUCL(TRIM_ADAPT.out)
    //TRIM_PRIMERS(TRIM_4_NUCL.out)
    //FASTQC_2 = FASTQC(TRIM_PRIMERS.out)
    //FQ2(TRIM_PRIMERS.out)
    dir_name = "before_extract"
    KRAKEN2(TRIM_4_NUCL.out, dir_name)
    BRACKEN_UPDATED(KRAKEN2.out.sample_id__report, dir_name)
  emit: 
     FQ1.out | concat(KRAKEN2.out.both_report_output) | concat(BRACKEN_UPDATED.out.tsv) | collect
     //FQ1.out | concat(FQ2.out) | collect
     //FQ1.out | concat(FQ2.out) | concat(KRAKEN2.out) | collect
     //FASTQC.out | concat(KRAKEN2.out) | collect
}

workflow TAXONOMY_ANALYSIS_WITH_REF_ALIGN_BASE {
  take:
    read_pairs_ch
    new_panel_fasta
  main:
    dir_name_0 = "raw_data"
    FQ1(read_pairs_ch, dir_name_0)
    TRIM_ADAPT(read_pairs_ch)
    TRIM_4_NUCL(TRIM_ADAPT.out)
    BWA_INDEX1(tuple(
      file(params.new_panel_fasta), "new_panel")
      )
    combine_ch = TRIM_4_NUCL.out.combine(BWA_INDEX1.out)
    //combine_ch.view()
    BWA_MEM_BAM_SORT(combine_ch)
    SAMTOOLS_STATS(BWA_MEM_BAM_SORT.out.sample_id__bam)     

    dir_name = "before_extract"
    KRAKEN2(TRIM_4_NUCL.out, dir_name)
    BRACKEN_UPDATED(KRAKEN2.out.sample_id__report, dir_name)
  emit: 
     FQ1.out | concat(KRAKEN2.out.both_report_output) | concat(BRACKEN_UPDATED.out.tsv) | collect
}

workflow TAXONOMY_ANALYSIS_WITH_REF_ALIGN {
  take:
    read_pairs_ch
    old_panel_fasta
    new_panel_fasta
  main:
    dir_name_0 = "raw_data"
    FQ1(read_pairs_ch, dir_name_0)
    TRIM_ADAPT(read_pairs_ch)
    TRIM_4_NUCL(TRIM_ADAPT.out)
    read_pairs_ch_old_panel = TRIM_4_NUCL.out.filter { sample_id, files -> sample_id.contains('old')  }
    read_pairs_ch_new_panel = TRIM_4_NUCL.out.filter { sample_id, files -> !sample_id.contains('old')  }
    BWA_INDEX1(tuple(
      file(params.old_panel_fasta), "old_panel")
      )

    BWA_INDEX2(tuple(
      file(params.new_panel_fasta), "new_panel")
      )
    combine_ch_old = read_pairs_ch_old_panel.combine(BWA_INDEX1.out)//.view()
    combine_ch_new = read_pairs_ch_new_panel.combine(BWA_INDEX2.out)//.view()
    combine_ch = combine_ch_old.mix(combine_ch_new)
    //combine_ch.view()
    BWA_MEM_BAM_SORT(combine_ch)
    SAMTOOLS_STATS(BWA_MEM_BAM_SORT.out.sample_id__bam)     

    dir_name = "before_extract"
    KRAKEN2(TRIM_4_NUCL.out, dir_name)
    BRACKEN_UPDATED(KRAKEN2.out.sample_id__report, dir_name)
    

 
    
    
    //old_panel_fasta_with_name = old_panel_fasta.map { f -> 
    //                                                      prefix = f.baseName  // Имя файла без расширения
    //                                                     tuple(f, prefix) }.view()
    //new_panel_fasta_with_name = new_panel_fasta.map { f -> tuple(f, f.name) }.view()

    
    //BWA_INDEX2(new_panel_fasta_with_name)

  emit: 
     FQ1.out | concat(KRAKEN2.out.both_report_output) | concat(BRACKEN_UPDATED.out.tsv) | collect
     //FQ1.out | concat(FQ2.out) | collect
     //FQ1.out | concat(FQ2.out) | concat(KRAKEN2.out) | collect
     //FASTQC.out | concat(KRAKEN2.out) | collect
}