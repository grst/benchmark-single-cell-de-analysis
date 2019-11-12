#!/usr/bin/env nextflow

process run_diffxpy_ncounts_0613 {
    conda "envs/diffxpy_ncounts_0.6.13.yml"
    publishDir "results/diffxpy_ncounts"
    cpus 4

    input:
        file 'notebook.Rmd' from file('diffxpy_test/diffxpy_ncounts.Rmd')
        file "*" from Channel.fromPath("diffxpy_test/*.tsv").collect()

    output:
        file "diffxpy_ncounts.html" into diffxpy_ncounts_html_0613

    """
    reportsrender notebook.Rmd diffxpy_ncounts.html --cpus=${task.cpus}
    """
}
process run_diffxpy_ncounts_071 {
    conda "envs/diffxpy_ncounts.yml"
    publishDir "results/diffxpy_ncounts"
    cpus 4

    input:
        file 'notebook.Rmd' from file('diffxpy_test/diffxpy_ncounts.Rmd')
        file "*" from Channel.fromPath("diffxpy_test/*.tsv").collect()

    output:
        file "diffxpy_ncounts.html" into diffxpy_ncounts_html_071

    """
    reportsrender notebook.Rmd diffxpy_ncounts.html --cpus=${task.cpus}
    """
}
