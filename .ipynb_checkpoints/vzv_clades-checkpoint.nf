process AGG_CLADES {
    //conda 'multiqc'
    conda "/export/home/agletdinov/mambaforge/envs/reat"
    
    tag "Clades"
    publishDir "${params.outdir}/res_clades", mode:'copy'

    input:
    path(script4)
    path("*")

    output:
    path("clades_mqc.tsv")

    script:
    """
    python3 ${script4} -i .
    """
}