process MULTIQC {
    //conda 'multiqc'
    conda "${CONDA_PREFIX_1}/envs/multiqc"
    
    tag "MultiQC"
    publishDir "${params.outdir}/multiqc", mode:'copy'

    input:
    path("*")

    output:
    path("multiqc_report.html")

    script:
    """
    multiqc . --cl-config "{max_table_rows: 2000}"
    """
}