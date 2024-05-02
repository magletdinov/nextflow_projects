process FASTQC {
    conda = 'bioconda::fastqc'

    tag "FastQC on ${sample_id}"
    publishDir "${params.outdir}/fastqc", mode: "copy"

    input:
    tuple val(sample_id), path(reads)

    output:
    path("fastqc_${sample_id}_logs")

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}

process TRIM_ADAPT {
    conda = 'bioconda::fastp'

    tag "Fastp on ${sample_id}"
    publishDir "${params.outdir}/fastp_adapt", mode: "copy"
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path('*fastq.gz')
    
    cpus 3
    //maxForks params.maxForks
    //errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    //maxRetries 5
    
    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} -o ${sample_id}_t_a_R1.fastq.gz -O  ${sample_id}_t_a_R2.fastq.gz --thread ${task.cpus} --adapter_fasta ${params.adapters}
    """
}
  
process TRIM_4_NUCL {
    conda = 'bioconda::cutadapt'

    tag "Cutadapt on ${sample_id}"
    publishDir "${params.outdir}/cutadapt_trim_n", mode: "copy"
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path('*fastq.gz')
    
    cpus 3

    //maxForks params.maxForks
    //errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    //maxRetries 5

    script:
    """
    cutadapt -u 4 -u -4 -U 4 -U -4 -o ${sample_id}_t_4_R1.fastq.gz -p ${sample_id}_t_4_R2.fastq.gz \
        ${reads[0]} ${reads[1]} \
        -j ${task.cpus}
    """
}


process TRIM_PRIMERS {
    conda = 'bioconda::cutadapt'
    
    tag "Cutadapt (primers) on ${sample_id}"
    publishDir "${params.outdir}/cutadapt_primers", mode: "copy"
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path('*fastq.gz')
    
    cpus 3
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
   
    cpus 5
    //maxForks params.maxForks
    
    script:
    """
    bwa mem -t ${task.cpus} "bwa_index/${idxbase}" ${reads[0]} ${reads[1]}|
    samtools view --threads ${task.cpus} -1|
    samtools sort --threads ${task.cpus} -o "${sample_id}.aln.sorted.bam" 
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
   
    cpus 1
    //maxForks params.maxForks
    
    script:
    """
    samtools index -@ ${task.cpus} ${bam}
    """
}

process SAMTOOLS_CONSENSUS {
    conda 'bioconda::samtools'
    
    tag "Samtools consensus on ${sample_id}"
    //publishDir "${params.outdir}/consensus", mode:'copy'
    
    input:
    tuple val(sample_id), path(bam)
    
    //errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    //maxRetries 5
    
    output:
    tuple val(sample_id), path("*.fasta")
   
    cpus 3
    //maxForks params.maxForks
    
    script:
    """
    samtools consensus -aa -d 20 -m simple -o ${sample_id}.fasta ${bam}
    """
}

process SAMTOOLS_STATS {
    conda 'bioconda::samtools'
    
    tag "Samtools consensus on ${sample_id}"
    publishDir "${params.outdir}/statistics/${sample_id}", mode:'copy'
    
    input:
    tuple val(sample_id), path(bam)
    
    //errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    //maxRetries 5
    
    output:
    path("*")
   
    cpus 3
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
    conda 'conda-forge::biopython'
     
    tag "Rename fasta id for ${sample_id}"
    publishDir "${params.outdir}/consensus_new_head", mode:'copy'
    
    input:
    tuple val(sample_id), path(consensus)
    path(script)    
    //errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    //maxRetries 5
    
    output:
    path("${sample_id}_renamed.fasta")
   
    cpus 1
    //maxForks params.maxForks
    
    script:
    """
    python3 ${script} -i ${consensus} -o ${sample_id}_renamed.fasta
    """
}

process KRAKEN2 {
    //conda 'kraken2'
    conda "/export/home/agletdinov/mambaforge/envs/kraken2"
    cpus 100
    memory 700.GB
    maxForks 1
    
    tag "Kraken2 on ${sample_id}"
    publishDir "${params.outdir}/kraken2", mode: "copy"
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path('*.kraken')
    
    script:
    """
    kraken2 --db ${params.kraken2db}  --threads ${task.cpus} --gzip-compressed --report ${sample_id}_report.kraken --paired ${reads[0]} ${reads[1]} > ${sample_id}_kraken.txt
    """
}