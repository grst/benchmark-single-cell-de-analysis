#!/bin/bash
nextflow run main.nf \
  -profile icbi \
  -w /data/scratch/sturm/projects/2021/sc-de-benchmark/work
  -resume
