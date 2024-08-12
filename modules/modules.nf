process FASTQC {
    //conda = 'bioconda::fastqc'
    conda "/export/home/agletdinov/mambaforge/envs/multiqc"

    tag "FastQC on ${sample_id}"
    publishDir "${params.outdir}/fastqc", mode: "copy"

    input:
    tuple val(sample_id), path(reads)

    output:
    path("*")

    script:
    """
    fastqc ${reads}
    multiqc . --filename ${sample_id}
    """
}

process TRIM_ADAPT {
    conda = '/export/home/agletdinov/mambaforge/envs/fastp'

    tag "Fastp on ${sample_id}"
    publishDir "${params.outdir}/fastp_adapt", mode: "copy"
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path('*fastq.gz')
    
    cpus 20
    //maxForks params.maxForks
    //errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    //maxRetries 5
    
    script:
    if (params.singleEnd) {
        """
        fastp -i ${reads} -o ${sample_id}_t_a.fastq.gz --thread ${task.cpus} --adapter_fasta ${params.adapters}
        """
    } else {
        """
        fastp -i ${reads[0]} -I ${reads[1]} -o ${sample_id}_t_a_R1.fastq.gz -O  ${sample_id}_t_a_R2.fastq.gz --thread ${task.cpus} --adapter_fasta ${params.adapters}
        """
    }
}
  
process TRIM_4_NUCL {
    conda = 'bioconda::cutadapt'

    tag "Cutadapt on ${sample_id}"
    //publishDir "${params.outdir}/cutadapt_trim/${mode}_n", mode: "copy"
    publishDir "${params.outdir}/cutadapt_trim", mode: "copy"
    
    input:
    tuple val(sample_id), path(reads)
    //each mode

    output:
    //tuple val("${sample_id}_${mode}"), path('*fastq.gz')
    tuple val("${sample_id}"), path('*fastq.gz')
    
    cpus 20

    //maxForks params.maxForks
    //errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    //maxRetries 5

    script:
    if (params.singleEnd) {
        """
        cutadapt -u 4 -u -4 -o ${sample_id}_t_4.fastq.gz \
        ${reads} -j ${task.cpus}
        """
    } else {
        """
        cutadapt -u 4 -u -4 -U 4 -U -4 -o ${sample_id}_t_4_R1.fastq.gz -p ${sample_id}_t_4_R2.fastq.gz \
        ${reads[0]} ${reads[1]} \
        -j ${task.cpus}
        """
    }
}


process TRIM_PRIMERS {
    conda = 'bioconda::cutadapt'
    
    tag "Cutadapt (primers) on ${sample_id}"
    publishDir "${params.outdir}/cutadapt_primers", mode: "copy"
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path('*fastq.gz')
    
    cpus 20
    //maxForks params.maxForks
    //errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    //maxRetries 5

    script:
    """
    cutadapt -a file:${params.primer_a} -A file:${params.primer_a} -g file:${params.primer_g} -G file:${params.primer_g} \
    -o ${sample_id}_t_pr_R1.fastq.gz -p ${sample_id}_t_pr_R2.fastq.gz \
        ${reads[0]} ${reads[1]} \
        -j ${task.cpus}
    """
}

process REMOVE_HOST {
    conda = '/export/home/agletdinov/mambaforge/envs/bowtie2'

    tag "Remove host from ${sample_id}"
    publishDir "${params.outdir}/remove_host", mode:'copy'
    cpus 40
    
    input:
    tuple val(sample_id), path(reads)
    path(db)

    output:
    tuple val(sample_id), path('*fastq.gz')
    
    script:
    """
    idx_base=\$(find ${db}/ -name '*.bt2' | awk -F \".\" '{print \$1 | \"sort -u\"}')

    bowtie2 -p ${task.cpus} -x \${idx_base} \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        -S ${sample_id}_mapped_and_unmapped.sam
    samtools view -bS --threads ${task.cpus} ${sample_id}_mapped_and_unmapped.sam > ${sample_id}_mapped_and_unmapped.bam

    samtools view -b -f 12 -F 256 ${sample_id}_mapped_and_unmapped.bam > ${sample_id}_bothReadsUnmapped.bam

    samtools sort -n -m 10G -@ 2 ${sample_id}_bothReadsUnmapped.bam -o ${sample_id}_bothReadsUnmapped_sorted.bam

    samtools fastq -@ 10 ${sample_id}_bothReadsUnmapped_sorted.bam \
                    -1 ${sample_id}_host_removed_R1.fastq.gz \
                    -2 ${sample_id}_host_removed_R2.fastq.gz \
                    -0 /dev/null -s /dev/null -n
    """
}


process BWA_INDEX {
    conda = 'bioconda::bwa'

    tag "Bwa index for ${params.key_areas}"
    publishDir "${params.outdir}/bwa_index", mode:'copy'
    
    input:
    path(ref_fasta)
    val(prefix)
    
    //errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    //maxRetries 5
    
    output:
    tuple val(prefix), path("${prefix}.{ann,amb,sa,bwt,pac}")

    """
    bwa index \\
        -p "${prefix}" \\
        "${ref_fasta}"
    """
}

process BWA_INDEX_FULL_GENOME {
    conda = '/export/home/agletdinov/mambaforge/envs/bwa'

    tag "Bwa index for ${genome}"
    publishDir "${params.outdir}/bwa_index/${genome}", mode:'copy'
    
    input:
    each genome
    //errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    //maxRetries 5
    
    output:
    tuple val("${genome}"), path("${genome}.{ann,amb,sa,bwt,pac}")

    """
    bwa index \\
        -p "${genome}" \\
        "${params.genome_dict_taxid[genome]}"
    """
}

process BWA_MEM_BAM_SORT {
    conda 'bioconda::bwa bioconda::samtools'
    
    tag "Bwa mem on ${sample_id}"
    publishDir "${params.outdir}/bam", mode:'copy'
    
    input:
    tuple val(sample_id), path(reads)
    tuple val(idxbase), path("bwa_index/*")
    
    //errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    //maxRetries 5
    
    output:
    tuple val(sample_id), path("${sample_id}.aln.sorted.bam")
   
    cpus 20
    //maxForks params.maxForks
    
    script:
    """
    bwa mem -t ${task.cpus} "bwa_index/${idxbase}" ${reads[0]} ${reads[1]}|
    samtools view --threads ${task.cpus} -1|
    samtools sort --threads ${task.cpus} -o "${sample_id}.aln.sorted.bam" 
    """
}

process BWA_MEM_BAM_SORT_FULL_GENOME {
    conda 'bioconda::bwa bioconda::samtools'
    
    tag "Bwa mem on ${sample_id}_${idxbase}"
    publishDir "${params.outdir}/bam", mode:'copy'
    
    input:
    tuple val(sample_id), path(reads)
    //tuple val(idxbase), path("bwa_index/${idxbase}/*")
    each genome
    //errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    //maxRetries 5
    
    output:
    tuple val("${sample_id}_${genome}"), path("${sample_id}_${genome}.aln.sorted.bam")
   
    cpus 20
    //maxForks params.maxForks
    
    script:
    """
    bwa mem -t ${task.cpus} "${params.bwa_index}/${genome}/${genome}" ${reads[0]} ${reads[1]} |
    samtools view --threads ${task.cpus} -1 |
    samtools sort --threads ${task.cpus} -o "${sample_id}_${genome}.aln.sorted.bam" 
    """
}


process SAMTOOLS_INDEX {
    conda 'bioconda::samtools'
    
    tag "Samtools index on ${sample_id}"
    publishDir "${params.outdir}/bam", mode:'copy'
    
    input:
    tuple val(sample_id), path(bam)
    
    //errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    //maxRetries 5
    
    output:
    tuple val(sample_id), path("${sample_id}.aln.sorted.bam.bai")
   
    cpus 5
    //maxForks params.maxForks
    
    script:
    """
    samtools index -@ ${task.cpus} ${bam}
    """
}

process SAMTOOLS_CONSENSUS {
    conda 'bioconda::samtools'
    
    tag "Samtools consensus on ${sample_id}"
    publishDir "${params.outdir}/consensus", mode:'copy'
    
    input:
    tuple val(sample_id), path(bam)
    
    //errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    //maxRetries 5
    
    output:
    tuple val(sample_id), path("*.fasta")
   
    cpus 20
    //maxForks params.maxForks
    
    script:
    """
    samtools consensus -aa -d 20 -m simple -o ${sample_id}.fasta ${bam}
    """
}

process SAMTOOLS_CONSENSUS_LITE {
    conda 'bioconda::samtools'
    
    tag "Samtools consensus on ${sample_id}"
    publishDir "${params.outdir}/consensus", mode:'copy'
    
    input:
    tuple val(sample_id), path(bam)
    
    //errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    //maxRetries 5
    
    output:
    tuple val(sample_id), path("*.fasta")
   
    cpus 20
    //maxForks params.maxForks
    
    script:
    """
    samtools consensus -m simple -o ${sample_id}.fasta ${bam}
    """
}

process SAMTOOLS_STATS {
    conda 'bioconda::samtools'
    
    tag "Samtools stats on ${sample_id}"
    publishDir "${params.outdir}/statistics/${sample_id}", mode:'copy'
    
    input:
    tuple val(sample_id), path(bam)
    
    //errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    //maxRetries 5
    
    output:
    path("*")
   
    cpus 20
    //maxForks params.maxForks
    
    
    //samtools depth ${bam} > depth_${sample_id}.txt
    //samtools flagstat ${bam} > flagstat_${sample_id}.txt
    //samtools view -c -F 4 ${bam} > proportion_of_mapped_reads_${sample_id}.txt    
    script:
    """
    samtools stats -@ ${task.cpus} ${bam} > samtools_stats_${sample_id}.txt
    samtools idxstats ${bam} > samtools_idxstats_${sample_id}.txt

    """
}


process RENAME_FASTA_ID {
    //conda 'conda-forge::biopython'
    conda "/export/home/agletdinov/mambaforge/envs/reat"
     
    tag "Rename fasta id for ${sample_id}"
    publishDir "${params.outdir}/consensus_new_head", mode:'copy'
    
    input:
    tuple val(sample_id), path(consensus)
    path(script1)    
    //errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    //maxRetries 5
    
    output:
    tuple val(sample_id), path("*.fasta")
   
    cpus 3
    //maxForks params.maxForks
    
    script:
    """
    python3 ${script1} -i ${consensus} -o ${sample_id}_renamed.fasta
    """
}

process KRAKEN2 {
    //conda 'kraken2'
    conda "/export/home/agletdinov/mambaforge/envs/kraken2"
    maxForks 1
    cpus 45

    tag "Kraken2 on ${sample_id}"
    publishDir "${params.outdir}/kraken2/${dir_name}", mode: "copy"
    
    input:
    tuple val(sample_id), path(reads)
    val (dir_name)
    
    output:
    path('*.kraken'), emit: report
    tuple val(sample_id), path('*report.kraken'), emit: id_report
    tuple val(sample_id), path('*output.kraken'), emit: id_output
    
    script:
    if (params.singleEnd) {
        """
        kraken2 --db ${params.kraken2db} --threads ${task.cpus} --memory-mapping\
         --output ${sample_id}_output.kraken --report ${sample_id}_report.kraken ${reads} > ${sample_id}_kraken.txt
        """
    } else {
        """
        kraken2 --db ${params.kraken2db} --threads ${task.cpus} --memory-mapping\
         --output ${sample_id}_output.kraken --report ${sample_id}_report.kraken --paired ${reads[0]} ${reads[1]} > ${sample_id}_kraken.txt
        """
    }
}


process KRAKEN2_FASTA {
    //conda 'kraken2'
    conda "/export/home/agletdinov/mambaforge/envs/kraken2"
    //maxForks 1
    errorStrategy 'ignore'
    cpus 45
        
    tag "Kraken2 on ${sample_id}"
    publishDir "${params.outdir}/kraken2/${dir_name}", mode: "copy"
    
    input:
    tuple val(sample_id), path(reads)
    val (dir_name)
    output:
    path('*.kraken'), emit: report
    tuple val(sample_id), path('*report.kraken'), emit: id_report
    tuple val(sample_id), path('*output.kraken'), emit: id_output
    
    script:
    """
    kraken2 --db ${params.kraken2db} --memory-mapping --threads ${task.cpus}  --output ${sample_id}_output.kraken --report ${sample_id}_report.kraken ${reads} > ${sample_id}_kraken.txt
    """
}

process BRACKEN {
    //conda 'kraken2'
    conda "/export/home/agletdinov/mambaforge/envs/bracken"
    //maxForks 1
    cpus 20

    tag "Bracken on ${sample_id}"
    publishDir "${params.outdir}/bracken", mode: "copy"
    
    input:
    tuple val(sample_id), path(kraken_report)
    
    output:
    path('*.tsv')
    
    script:
    """
    bracken -l S -t ${task.cpus} -d ${params.kraken2db} -i ${kraken_report} -o ${sample_id}_bracken_S_mqc.tsv
    """
}

process BRACKEN_EACH {
    //conda 'kraken2'
    conda "/export/home/agletdinov/mambaforge/envs/bracken"
    //maxForks 1
    cpus 20
    

    tag "Bracken on ${sample_id}"
    publishDir "${params.outdir}/bracken/${params.bracken_settings_dict[resolution]}", mode: "copy"
    
    input:
    tuple val(sample_id), path(kraken_report)
    each resolution
    
    output:
    path('*.tsv')

    
    script:
    """
    bracken -l ${resolution} -t ${task.cpus} -d ${params.kraken2db} -i ${kraken_report} -o ${sample_id}_bracken_${resolution}_mqc.tsv
    """
}

process EXTRACT_KRAKEN_READS {
    //conda 'kraken2'
    conda "/export/home/agletdinov/mambaforge/envs/kraken2"
    //maxForks 1
    cpus 20
    
    tag "Extract kraken reads on ${sample_id}"
    //publishDir "${params.outdir}/krakentools/extract_${taxid}", mode: "copy"
    publishDir "${params.outdir}/krakentools/extract_taxids/${sample_id}", mode: "copy"
    
    input:
    tuple val(sample_id), path(reads), path(kraken_output), path(kraken_report)
    each taxid
    //тут нельзя писать println!!! (ниже when)
    when:
    params.taxid_dict.containsKey(taxid) && sample_id in params.taxid_dict[taxid] && params.krakentools_flag == true
    output:
    tuple val("${sample_id}-${taxid}"), path('*fast*'), emit: id_files    
    tuple val("${sample_id}-${taxid}"), path('*fasta'), emit: id_fasta
    tuple val("${sample_id}-${taxid}"), path('*fastq*'), emit: id_fastq
    
    script:
    if (params.singleEnd) {
        """
        extract_kraken_reads.py -k ${kraken_output} -r ${kraken_report} -t ${taxid} \
            -s ${reads} \
            -o ${sample_id}-${taxid}-withchildren_1.fastq --include-children --fastq-output
        seqtk seq -a ${sample_id}-${taxid}-withchildren_1.fastq > ${sample_id}-${taxid}-withchildren_1.fasta
        """
    } else {
        """
        extract_kraken_reads.py -k ${kraken_output} -r ${kraken_report} -t ${taxid} \
            -s ${reads[0]} -s2 ${reads[1]} \
            -o ${sample_id}-${taxid}-withchildren_1.fastq -o2 ${sample_id}-${taxid}-withchildren_2.fastq --include-children --fastq-output
        seqtk seq -a ${sample_id}-${taxid}-withchildren_1.fastq > ${sample_id}-${taxid}-withchildren_1.fasta
        seqtk seq -a ${sample_id}-${taxid}-withchildren_2.fastq > ${sample_id}-${taxid}-withchildren_2.fasta
        """
    }
}

process EXTRACT_KRAKEN_READS_FASTA {
    //conda 'kraken2'
    conda "/export/home/agletdinov/mambaforge/envs/kraken2"
    //maxForks 1
    //errorStrategy 'ignore'
    cpus 20
    

    tag "Extract kraken reads on ${sample_id}"
    //publishDir "${params.outdir}/krakentools/extract_${taxid}", mode: "copy"
    publishDir "${params.outdir}/krakentools/extract_taxids_contigs/${sample_id}", mode: "copy"
    
    input:
    tuple val(sample_id), path(reads), path(kraken_output), path(kraken_report)
    val(taxid_dict_v2)
    each taxid
    //"${sample_id}".contains(taxid) && params.taxid_dict_v2[taxid].contains("${sample_id}-${taxid}")
    //"${sample_id}-${taxid}" in params.taxid_dict_v2[taxid]
    //println "${sample_id}-${taxid}"
    //println params.taxid_dict_v2
    //println "${params.taxid_dict_v2["3050337"]}"
    //"${sample_id}-${taxid}" in params.taxid_dict_v2["${taxid}"]
    when:
    sample_id + '-' + taxid in taxid_dict_v2[taxid]
    output:
    tuple val("${sample_id}-${taxid}"), path('*fasta')
    
    script: 
    """
    extract_kraken_reads.py -k ${kraken_output} -r ${kraken_report} -t ${taxid} \
    -s ${reads} -o ${sample_id}-${taxid}-withchildren-contigs.fasta --include-children
    """
}

process EXTRACT_KRAKEN_READS_TAXID_LIST {
    conda "/export/home/agletdinov/mambaforge/envs/kraken2"
    cpus 20
    
    tag "Extract kraken reads on ${sample_id}"
    publishDir "${params.outdir}/krakentools/extract_taxids_eupath/${sample_id}", mode: "copy"
    
    input:
    tuple val(sample_id), path(reads), path(kraken_output), path(kraken_report)
    val taxid_list
    
    output:
    tuple val(sample_id), path('*fast*'), emit: id_files    
    //tuple val(sample_id), path('*[12].fasta*'), emit: id_fasta
    tuple val(sample_id), path('*[12].fastq'), emit: id_fastq
    //tuple val(sample_id), path('*eupath-withchildren.fastq*'), emit: id_onefastq
    
    
    script:
    if (params.singleEnd) {
        """
        extract_kraken_reads.py -k ${kraken_output} -r ${kraken_report} -t ${taxid_list} \
            -s ${reads} -fastq-output\
            -o ${sample_id}-eupath-withchildren_1.fastq --include-children -
        seqtk seq -a ${sample_id}-eupath-withchildren_1.fastq > ${sample_id}-eupath-withchildren_1.fasta
        """
    } else {
        """
        extract_kraken_reads.py -k ${kraken_output} -r ${kraken_report} -t 885310 5759 885312 885313 885314 414452 665079 27973 333668 36329 5858 864141 85471 5855 383379 5811 65658 5755 43142 29196 61605 196912 65661 211522 32599 65662 32600 5757 28015 294381 41668 5763 27987 237895 756076 1169474 1169540 911350 5057 5059 451804 330879 1033177 890382 28583 89421 113399 227321 380704 425011 5061 510516 166112 332648 403673 237561 294748 367775 284593 404692 454286 569365 235443 214684 178876 578454 469472 490068 490069 490065 490066 469470 469471 264361 215243 1279085 5518 1089452 426428 1089451 660027 660029 334819 559515 5037 747725 318829 76775 100951 1231661 367110 510953 1223555 1223556 763407 763924 273507 1048749 588596 1223558 639000 164328 500485 1223559 431595 1223560 559292 402676 483514 284812 999809 37727 28564 284591 4952 336722 941442 5741 1394984 2049356 723287 986730 989654 989655 284813 907965 876142 1178016 657448 1485682 5866 5865 32595 1133968 83556 353154 5874 1537102 869250 880535 5823 720590 880536 31271 208452 1120755 1274352 57266 5833 57267 137071 647221 5849 5851 5850 880534 5854 73239 88456 5801 51314 84963 5804 44415 51315 51316 5802 572307 42890 1074873 507601 508771 943120 432359 1463230 5659 420245 5661 981087 435258 347515 929439 5679 157538 5684 5689 1470209 679716 5691 185431 1068625 5693 353153 85056 5697 1055687 1257118 370354 885318 885319 885315 885311 370355 1076696 117008 93969 441375 353152 1603071 857276 110365 690307 112090 767769 1392248 602072 344612 331117 332952 1160497 157072 1137211 578462 1353008 1392255 1392256 1392250 1036612 341663 767770 1036611 1073089 1073090 498019 86049 240176 294750 1231522 1220924 1296111 1296103 1296105 1296108 794803 45357 396776 246410 306902 283643 222929 454284 443226 294747 212818 1229664 2502994 1442368 1227346 948311 447093 544712 544711 2059318 1296100 1296121 1296119 41688 1220926 425265 747676 242507 510951 482561 502780 403677 42068 502779 611791 761204 67593 246409 747089 563466 1398154 1156394 771870 695850 645134 1397361 441960 578456 431241 441959 413071 237631 336963 184922 598745 658858 5742 348837 1240240 1288291 1866961 1003232 481877 1081671 1355160 646526 1081669 1913371 578461 578460 40302 1805483 881290 935791 1354746 146866 1358809 72359 993615 948595 189622 5857 1237626 138298 54757 483139 99158 1074872 943122 1194599 943167 1130820 943118 943119 935652 943121 1130821 933077 412133 75058 1416333 85057 71804 429131 67003\
            -s ${reads[0]} -s2 ${reads[1]} --fastq-output\
            -o ${sample_id}-eupath-withchildren_1.fastq -o2 ${sample_id}-eupath-withchildren_2.fastq --include-children

        
        """
        //seqfu interleave -1 ${sample_id}-eupath-withchildren_1.fastq.gz -2 ${sample_id}-eupath-withchildren_2.fastq.gz > ${sample_id}-eupath-withchildren.fastq.gz
        //seqtk seq -a ${sample_id}-eupath-withchildren_1.fastq > ${sample_id}-eupath-withchildren_1.fasta
        //seqtk seq -a ${sample_id}-eupath-withchildren_2.fastq > ${sample_id}-eupath-withchildren_2.fasta
        //cat ${sample_id}-eupath-withchildren_1.fasta ${sample_id}-eupath-withchildren_2.fasta > ${sample_id}-eupath-withchildren.fasta
        //gzip ${sample_id}-eupath-withchildren_1.fastq ${sample_id}-eupath-withchildren_1.fastq
    }
}


process EXTRACT_KRAKEN_READS_TAXID {
    conda "/export/home/agletdinov/mambaforge/envs/kraken2"
    cpus 20
    
    tag "Extract kraken reads on ${sample_id}"
    publishDir "${params.outdir}/krakentools/extract_taxids_eupath/${sample_id}", mode: "copy"
    
    input:
    tuple val(sample_id), path(reads), path(kraken_output), path(kraken_report)
    val taxid
    
    output:
    tuple val(sample_id), path('*fast*'), emit: id_files    
    tuple val(sample_id), path('*[12].fasta*'), emit: id_fasta

    when:
    params.extract_taxid == true   
    
    script:
    """
    extract_kraken_reads.py -k ${kraken_output} -r ${kraken_report} -t ${taxid}\
        -s ${reads[0]} -s2 ${reads[1]} \
        -o ${sample_id}-${taxid}-eupath-withchildren_1.fasta -o2 ${sample_id}-${taxid}-eupath-withchildren_2.fasta --include-children
    """
}


process REMOVE_HUMAN {
    //conda 'kraken2'
    conda "/export/home/agletdinov/mambaforge/envs/kraken2"
    //maxForks 1
    cpus 20
    
    tag "Remove human reads from ${sample_id} by kraken2"
    publishDir "${params.outdir}/krakentools/remove_human", mode: "copy"
    
    input:
    tuple val(sample_id), path(reads), path(kraken_output), path(kraken_report)
    
    output:
    tuple val(sample_id), path('*fast*'), emit: id_files    
    tuple val(sample_id), path('*fasta'), emit: id_fasta
    tuple val(sample_id), path('*without_human.fasta'), emit: id_onefasta
    tuple val(sample_id), path('*fastq*'), emit: id_fastq
    
    script:
    """
    extract_kraken_reads.py -k ${kraken_output} -r ${kraken_report} -t 9606 \
        -s ${reads[0]} -s2 ${reads[1]} \
        -o ${sample_id}-without_human-1.fastq -o2 ${sample_id}-without_human-2.fastq --include-children --fastq-output --exclude
    seqtk seq -a ${sample_id}-without_human-1.fastq > ${sample_id}-without_human-1.fasta
    seqtk seq -a ${sample_id}-without_human-2.fastq > ${sample_id}-without_human-2.fasta
    cat ${sample_id}-without_human-1.fasta ${sample_id}-without_human-2.fasta > ${sample_id}-without_human.fasta
    """
}


process IDENTIFY_CLADE {
    conda "/export/home/agletdinov/mambaforge/envs/reat"
     
    tag "Clade for ${sample_id}"
    publishDir "${params.outdir}/clades", mode:'copy'
    
    input:
    tuple val(sample_id), path(consensus)
    path(script2)
    path(script3)
    //errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    //maxRetries 5
    
    output:
    //path("clade_for_${sample_id}.log")
    path("*_clade.json")
   
    cpus 1
    //maxForks params.maxForks
    script:
    """
    mkdir temp_in_for_${sample_id}
    mv ${consensus} temp_in_for_${sample_id}
    python3 ${script2} -i temp_in_for_${sample_id} -o temp_out_for_${sample_id}
    python3 ${script3} -i temp_out_for_${sample_id}/full_report.json -o ${sample_id}_clade.json
    """
}

process UNICYCLER {
    //conda 'bioconda::unicycler'
    conda "/export/home/agletdinov/mambaforge/envs/unicycler"
    cpus 40
    //memory 500.GB
    //maxForks 2
    tag "Unicycler on ${sample_id}"
    publishDir "${params.outdir}/unicycler", mode:'copy'
    
    input:
    tuple val(sample_id), path(reads)

    output:
    path('*'), emit: report
    tuple val("${sample_id}_unicycler"), path("${sample_id}/assembly.fasta"), emit: id_contigs

    script:
    """
    mkdir uni_res_for_${sample_id}
    unicycler \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        -o ${sample_id} \
        -t ${task.cpus}
    """
}

process MEGAHIT {
    //conda 'bioconda::megahit'
    conda "/export/home/agletdinov/mambaforge/envs/megahit"
    memory 300.GB
    maxForks 3
    cpus 40

    tag "Megahit on ${sample_id}"
    publishDir "${params.outdir}/megahit", mode:'copy'
    
    input:
    tuple val(sample_id), path(reads)

    output:
    path('*'), emit: report
    tuple val(sample_id), path("${sample_id}/*fa"), emit: simple_id_contigs
    tuple val("${sample_id}-megahit"), path("${sample_id}/*fa"), emit: id_contigs

    script:
    """
    megahit -1 ${reads[0]} -2 ${reads[1]} \
        -o ${sample_id} --out-prefix ${sample_id} --memory 0.3 -t ${task.cpus} \
        --presets meta-sensitive
    """
}

process QUAST {
    conda "/export/home/agletdinov/mambaforge/envs/quast"
    cpus 40

    tag "Quast for ${sample_id}"
    publishDir "${params.outdir}/quast/${sample_id}", mode:'copy'
    
    input:
    tuple val(sample_id), path(contigs)
    each taxid
    
    //errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    //maxRetries 5
    
    output:
    tuple val(sample_id), path("*")

    """
    quast.py ${contigs} /
          -t ${task.cpus} -r ${params.genome_dict_taxid[taxid]}
    """
}

process METAPHLAN {
    //conda 'bioconda::metaphlan'
    conda "/export/home/agletdinov/mambaforge/envs/metaphlan"
    //memory = '1 MB'
    //maxForks 2
    cpus 40
    tag "Metaphlan on ${sample_id}"
    publishDir "${params.outdir}/metaphlan/${sample_id}", mode:'copy'
    
    input:
    tuple val(sample_id), path(contigs)

    output:
    tuple val(sample_id), path('*.txt'), emit: txt_report
    tuple val(sample_id), path('*'), emit: id_report
    path('*.txt'), emit: report

    script:
    """
    metaphlan ${contigs} --bowtie2db ${params.metaphlandb} --bowtie2out metagenome_${sample_id}.bowtie2.bz2 --nproc ${task.cpus} --input_type fasta -o profiled_metagenome_${sample_id}.txt
    """
}

process METAPHLAN_AGG {
    //conda 'bioconda::metaphlan'
    conda "/export/home/agletdinov/mambaforge/envs/metaphlan"
    //memory 500.GB
    //maxForks 2
    tag "Metaphlan on ${sample_id}"
    publishDir "${params.outdir}/metaphlan", mode:'copy'
    
    input:
    tuple val(sample_id), path(metaphlan_txt)

    output:
    path('merged_abundance_table.txt')
    
    script:
    """
    merge_metaphlan_tables.py ${metaphlan_txt} > merged_abundance_table.txt
    """
}

process DIAMOND {
    //conda 'bioconda::metaphlan'
    conda "/export/home/agletdinov/mambaforge/envs/diamond"
    //memory = '1 MB'
    //maxForks 2
    cpus 20
    tag "Diamond on ${sample_id}"
    publishDir "${params.outdir}/diamond/${sample_id}", mode:'copy'
    
    input:
    tuple val(sample_id), path(contigs)

    output:
    path('*.tsv')
 
    script:
    """
    diamond blastx \
        --very-sensitive \
        -d /export/home/public/tools/database/virus_nr_diamond.dmnd \
        --outfmt 6 qseqid sseqid pident evalue \
        -p ${task.cpus} \
        -q ${contigs} \
        -o ${sample_id}_matches.tsv \
        --max-target-seqs 10 \
        --evalue 1e-05 \
        -v \
        --log
    """
}

process blastn {
    //conda 'bioconda::metaphlan'
    conda "/export/home/agletdinov/mambaforge/envs/blast"
    //memory = '1 MB'
    //maxForks 2
    cpus 140
    tag "Blastn on ${sample_id}"
    publishDir "${params.outdir}/blastn/${sample_id}", mode:'copy'

    input:
    tuple val(sample_id), path(contigs)
    path db 

    output:
    tuple val(sample_id), path('*.blastn'), emit: id_report

    script:
    """
    blastn -db ${params.db}/nt -num_threads ${task.cpus} -out ${sample_id}.blastn -query ${contigs} -evalue 1e-03 -max_target_seqs 1 -max_hsps 1 -task megablast -outfmt "6 qaccver saccver sskingdoms sscinames salltitles staxids pident evalue"
    """
}

process BLASTN {
    //conda 'bioconda::metaphlan'
    conda "/export/home/agletdinov/mambaforge/envs/blastn"
    //memory = '1 MB'
    //maxForks 2
    cpus 45
    tag "Blastn on ${sample_id}"
    publishDir "${params.outdir}/blastn/${sample_id}", mode:'copy'
    
    input:
    tuple val(sample_id), path(contigs)
    path(blastnDB) 

    output:
    path('*.blastn')
 
    script:
    """
    blastn -db ${params.blastnDB}/nt -num_threads ${task.cpus} -out ${sample_id}.blastn -query ${contigs} -evalue 1e-03 -max_target_seqs 1 -max_hsps 1 -task megablast -outfmt "6 qaccver saccver sskingdoms sscinames salltitles staxids pident evalue"
    """
}


process blastn_parse{
    //conda 'bioconda::metaphlan'
    conda "/export/home/agletdinov/mambaforge/envs/reat"
    //memory = '1 MB'
    //maxForks 2
    cpus 5
    tag "Parse blastn result for ${sample_id}"
    publishDir "${params.outdir}/blastn_parse/${sample_id}", mode:'copy'
    
    input:
    tuple val(sample_id), path(report)
    path(to_nodes)

    output:
    path('*.tsv')
 
    script:
    kraken_taxid = sample_id.split('-')[-1]
    """
    #!/usr/bin/env python3
    import pandas as pd
    from pathlib import Path
    from collections import Counter

    nodes = pd.read_csv('${to_nodes}', sep="\t", header=None, usecols=[0, 2, 4], names=["tax_id", "parent tax_id", "tax_name"], index_col="tax_id")
    df_blastn = pd.read_csv('${report}', sep="\t", header=None, usecols=[0, 5], names=["seq_name", "taxid"])

    def create_path_to_root(taxid, nodes):
        try:
            path_to_root = []
            path_to_root.append(taxid)
            parent = nodes.loc[taxid]["parent tax_id"]
            if parent != 1:
                path_to_root.extend(create_path_to_root(parent, nodes))
            return path_to_root
        except:
            return []

    kraken_taxid = '${kraken_taxid}'
    df_blastn["path_to_root"] = df_blastn["taxid"].apply(lambda x: create_path_to_root(int(x), nodes))
    df_blastn["kraken_correct"] = df_blastn["taxid"].apply(lambda x: kraken_taxid in create_path_to_root(int(x), nodes))
    final_result = df_blastn["kraken_correct"].value_counts().idxmax()
    df_blastn.to_csv("extend_report.tsv", sep="\t")
    print(final_result)
    """
}


process metaSPAdes {
    //conda 'bioconda::metaphlan'
    conda "/export/home/agletdinov/mambaforge/envs/spades"
    //maxForks 2
    cpus 40
    tag "MetaSPAdes on ${sample_id}"
    publishDir "${params.outdir}/metaSPAdes/${sample_id}", mode:'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path('*'), emit: id_all
    tuple val(sample_id), path('scaffolds.fasta'), emit: id_scaffolds
 
    script:
    """
    spades.py -1 ${reads[0]} -2  ${reads[1]} -o .
    """
}