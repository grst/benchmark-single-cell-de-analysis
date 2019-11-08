#!/usr/bin/env nextflow

// tuples; [dataset, number of simulated samples/cells]
datasets = Channel.value([['islam', 80], ['trapnell', 150]])

/**
 * HACK: manually install zingeR into the respective
 * conda envs. Should create r-zinger conda package
 * at some point.
 */
process install_zinger_0613 {
    conda "envs/benchmark_diffxpy_0.6.13.yml"
    executor local

    output:
        val "installed" into installed_zinger_0613

    """
    Rscript -e "remotes::install_github('statOmics/zingeR')"
    """
}
process install_zinger_071 {
    conda "envs/benchmark_diffxpy_0.7.1.yml"
    executor local

    output:
        val "installed" into installed_zinger_071

    """
    Rscript -e "remotes::install_github('statOmics/zingeR')"
    """
}


process run_benchmark {
    conda "envs/benchmark_diffxpy_0.7.1.yml"
    publishDir "results/diffxpy_0.7.1"
    //32 = 4 * 8 (each thread will fork 8 times...)
    cpus 32

    input:
        val instlled from installed_zinger_071
        set val(dataset), val(n_samples) from datasets.flatMap()
        file 'lib' from file('analysis/lib')
        file 'notebook.Rmd' from file('analysis/benchmark.Rmd')

    output:
        file "benchmark_${dataset}.html" into report

    """
    reportsrender notebook.Rmd benchmark_${dataset}.html --cpus=4 \
        --params="dataset=${dataset} n_samples=${n_samples}"
    """
}


process run_benchmark_diffxpy_0613 {
    conda "envs/benchmark_diffxpy_0.6.13.yml"
    publishDir "results/diffxpy_0.6.13"
    cpus 32

    input:
        val installed from installed_zinger_0613
        set val(dataset), val(n_samples) from datasets.flatMap()
        file 'lib' from file('analysis/lib')
        file 'notebook.Rmd' from file('analysis/benchmark.Rmd')

    output:
        file "benchmark_${dataset}.html" into report2


    """
    reportsrender notebook.Rmd benchmark_${dataset}.html --cpus=4 \
        --params="dataset=${dataset} n_samples=${n_samples}"
    """
}


process run_diffxpy_ncounts {
    conda "envs/diffxpy_ncounts.yml"
    publishDir "results/diffxpy_ncounts"
    cpus 4

    input:
        file 'notebook.Rmd' from file('diffxpy_test/diffxpy_ncounts.Rmd')
        file "*" from Channel.fromPath("diffxpy_test/*.tsv")

    output:
        file "diffxpy_ncounts.html" into diffxpy_ncounts_html

    """
    reportsrender notebook.Rmd diffxpy_ncounts.html --cpus=${task.cpus}
    """
}
