#!/usr/bin/env nextflow-edge
nextflow.enable.dsl = 2

include { nxfVars } from "./nxfvars.nf"

process diffxpy0613 {
    def id = "diffxpy0.6.13"
    conda "/home/sturm/.conda/envs/2021-sc-de-benchmark-diffxpy0.6"
    publishDir params.output_dir, mode: 'link'
    cpus 22

    input:
        path(adata)

    output:
        path("${id}.csv"), emit: result

    script:
    output_file = "${id}.csv"
    """
    ${nxfVars(task)}
    jupytext --to ipynb --output ${id}.ipynb -k python3 $projectDir/analyses/01_run_de_methods/${id}.py
    jupyter nbconvert ${id}.ipynb --execute --to html  
    """
}

process diffxpy074 {
    def id = "diffxpy0.7.4"
    conda "/home/sturm/.conda/envs/2021-sc-de-benchmark-diffxpy07"
    publishDir params.output_dir, mode: 'link'
    cpus 22

    input:
        path(adata)

    output:
        path("${id}.csv"), emit: result

    script:
    output_file = "${id}.csv"
    """
    ${nxfVars(task)}
    jupytext --to ipynb --output ${id}.ipynb -k python3 $projectDir/analyses/01_run_de_methods/${id}.py
    jupyter nbconvert ${id}.ipynb --execute --to html  
    """

}

process edger {
    def id = "edger"
    conda "/home/sturm/.conda/envs/2021-sc-de-benchmark-edger"
    publishDir params.output_dir, mode: 'link'
    cpus 6

    input:
        path(adata)

    output:
        path("${id}.csv"), emit: result

    script:
    output_file = "${id}.csv"
    """
    ${nxfVars(task)}
    jupytext --to ipynb --output ${id}.ipynb -k python3 $projectDir/analyses/01_run_de_methods/${id}.py
    jupyter nbconvert ${id}.ipynb --execute --to html  
    """

}

process scVI {
    def id = "scVI"
    conda "/home/sturm/.conda/envs/pircher-sc-integrate2"
    clusterOptions '-V -S /bin/bash -q all.q@apollo-15'
    publishDir params.output_dir, mode: 'link'
    cpus 2

    input:
        path(adata)

    output:
        path("${id}.csv"), emit: result

    script:
    output_file = "${id}.csv"
    """
    ${nxfVars(task)}
    jupytext --to ipynb --output ${id}.ipynb -k python3 $projectDir/analyses/01_run_de_methods/${id}.py
    jupyter nbconvert ${id}.ipynb --execute --to html  
    """

}

process statsmodels {
    def id = "statsmodels_glm"
    conda "/home/sturm/.conda/envs/2021-sc-de-benchmark-diffxpy07"
    publishDir params.output_dir, mode: 'link'
    cpus 11

    input:
        path(adata)

    output:
        path("${id}.csv"), emit: result

    script:
    output_file = "${id}.csv"
    """
    ${nxfVars(task)}
    jupytext --to ipynb --output ${id}.ipynb -k python3 $projectDir/analyses/01_run_de_methods/${id}.py
    jupyter nbconvert ${id}.ipynb --execute --to html  
    """

}

workflow {
    diffxpy0613(file("./data/adata-myeloid.h5ad"))
    diffxpy074(file("./data/adata-myeloid.h5ad"))
    edger(file("./data/adata-myeloid.h5ad"))
    scVI(file("./data/adata-myeloid.h5ad"))
    statsmodels(file("./data/adata-myeloid.h5ad"))

}