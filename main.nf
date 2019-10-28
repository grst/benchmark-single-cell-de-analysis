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

    cpus 32

    """
    reportsrender rmd notebook.Rmd benchmark_islam.html --cpus=4
    """
}


workflow {
    benchmark_islam(file('analysis/lib'), file('analysis/benchmark_islam.Rmd'))
}
