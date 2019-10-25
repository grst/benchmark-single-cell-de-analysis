#!/usr/bin/env nextflow

//use dsl2 to support modules.
nextflow.preview.dsl=2


process benchmark_islam {
    input:
        file 'lib'
        file 'notebook.Rmd'

    output:
        file "benchmark_islam.html"

    publishDir "results"

    cpus 16

    """
    reportsrender rmd notebook.Rmd benchmark_islam.html --cpus=${task.cpus}
    """
}


workflow {
    benchmark_islam(file('analysis/lib'), file('analysis/benchmark_islam.Rmd'))
}
