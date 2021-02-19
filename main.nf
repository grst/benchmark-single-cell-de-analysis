nextflow.enable.dsl = 2

include { nxfVars } from "./nxfvars.nf"

process diffxpy0613 {
    def id = "diffxpy0613"
    conda "/home/sturm/.conda/envs/2021-sc-de-benchmark-diffxpy0.6"
    publishDir params.output_dir, mode: 'link'
    cpus 22

    input:
        path(adata)

    output:
        file("${id}.csv"), emit: result

    script:
    output_file = "${id}.csv"
    """
    jupytext --to ipynb $projectDir/01_run_de_methods/${id}.py
    jupyter nbconvert ${id}.ipynb --execute --to html  
    """
}

// process diffxpy074 {
//     conda "/home/sturm/.conda/envs/2021-sc-de-benchmark-diffxpy07"

// }

// process edger {
//     conda "/home/sturm/.conda/envs/2021-sc-de-benchmark-edger"

// }

// process scVI {
//     conda "/home/sturm/.conda/envs/pircher-sc-integrate2"
//     clusterOptions '-V -S /bin/bash -q all.q@apollo-15'

// }

// process statmodels {
//    conda "/home/sturm/.conda/envs/2021-sc-de-benchmark-diffxpy07"

// }

workflow {
    diffxpy0613(file("./data/adata-myeloid.h5ad"))

}