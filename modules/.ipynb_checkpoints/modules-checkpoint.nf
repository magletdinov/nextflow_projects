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
    """
    fastp -i ${reads[0]} -I ${reads[1]} -o ${sample_id}_t_a_R1.fastq.gz -O  ${sample_id}_t_a_R2.fastq.gz --thread ${task.cpus} --adapter_fasta ${params.adapters}
    """
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
    //maxForks 1
    cpus 40

    tag "Kraken2 on ${sample_id}"
    publishDir "${params.outdir}/kraken2", mode: "copy"
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path('*.kraken'), emit: report
    tuple val(sample_id), path('*report.kraken'), emit: id_report
    tuple val(sample_id), path('*output.kraken'), emit: id_output
    
    script:
    """
    kraken2 --db ${params.kraken2db}  --threads ${task.cpus} --gzip-compressed --output ${sample_id}_output.kraken --report ${sample_id}_report.kraken --paired ${reads[0]} ${reads[1]} > ${sample_id}_kraken.txt
    """
}

process KRAKEN2_FASTA {
    //conda 'kraken2'
    conda "/export/home/agletdinov/mambaforge/envs/kraken2"
    //maxForks 1
    errorStrategy 'ignore'
    cpus 40
        
    tag "Kraken2 on ${sample_id}"
    publishDir "${params.outdir}/kraken2_temp", mode: "copy"
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path('*.kraken'), emit: report
    tuple val(sample_id), path('*report.kraken'), emit: id_report
    tuple val(sample_id), path('*output.kraken'), emit: id_output
    
    script:
    """
    kraken2 --db ${params.kraken2db}  --threads ${task.cpus}  --output ${sample_id}_output.kraken --report ${sample_id}_report.kraken ${reads} > ${sample_id}_kraken.txt
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

    when:
    params.taxid_dict.containsKey(taxid) && sample_id in params.taxid_dict[taxid] && params.krakentools_flag == true

    output:
    tuple val("${sample_id}-${taxid}"), path('*fasta')
    
    script:
    """
    extract_kraken_reads.py -k ${kraken_output} -r ${kraken_report} -t ${taxid} \
    -s ${reads[0]} -s2 ${reads[1]} \
    -o ${sample_id}-${taxid}-withchildren_1.fasta -o2 ${sample_id}-${taxid}-withchildren_2.fasta --include-children
    """
}

process EXTRACT_KRAKEN_READS_FASTA {
    //conda 'kraken2'
    conda "/export/home/agletdinov/mambaforge/envs/kraken2"
    //maxForks 1
    cpus 20
    

    tag "Extract kraken reads on ${sample_id}"
    //publishDir "${params.outdir}/krakentools/extract_${taxid}", mode: "copy"
    publishDir "${params.outdir}/krakentools/extract_taxids_contigs/${sample_id}", mode: "copy"
    
    input:
    tuple val(sample_id), path(reads), path(kraken_output), path(kraken_report)
    each taxid
    
    when:
    //"${sample_id}".contains(taxid) && params.taxid_dict_v2[taxid].contains("${sample_id}-${taxid}")
    "k16_bird_S11-694014-megahit-694014" in params.taxid_dict_v2[taxid]
    //println "Checking ${sample_id}-${taxid} in ${params.taxid_dict_v2[taxid]}"
    println "Sample ID: ${sample_id}"
    //println "Reads file: ${reads}"
    //println "Kraken output: ${kraken_output}"
    //println "Kraken report: ${kraken_report}"
    //println "Tax ID: ${taxid}"
    //println "Tax ID dictionary: ${params.taxid_dict_v2}"    
    
    //println "Checking ${sample_id}-${taxid} in ${params.taxid_dict_v2[taxid]}"
    //params.taxid_dict_v2[taxid].contains(sample_id) == true
    //println "Checking ${sample_id} in ${params.taxid_dict_v2[taxid]}"
    //params.taxid_dict_v2[taxid].contains(sample_id)
    output:
    tuple val(sample_id), path('*fasta')
    
    script:
    """
    extract_kraken_reads.py -k ${kraken_output} -r ${kraken_report} -t ${taxid} \
    -s ${reads} -o ${sample_id}-${taxid}-withchildren-contigs.fasta --include-children
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
    //tuple val(sample_id), path('*'), emit: id_contigs
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