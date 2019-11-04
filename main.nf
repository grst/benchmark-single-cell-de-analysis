#!/usr/bin/env nextflow

//use dsl2 to support modules.
/* nextflow.preview.dsl=2 */


datasets = Channel.from(['islam', 80], ['trapnell', 150])

process run_benchmark {
    input:
        set val(dataset), val(n_samples) from datasets
        file 'lib' from file('analysis/lib')
        file 'notebook.Rmd' from file('analysis/benchmark.Rmd')

    output:
        file "benchmark_${dataset}.html" into report

    publishDir "results"

    cpus 32

    """
    reportsrender rmd notebook.Rmd benchmark_${dataset}.html --cpus=4 \
        --params="dataset=${dataset} n_samples=${n_samples}"
    """
}

