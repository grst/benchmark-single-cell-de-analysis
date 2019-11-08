# Experimenting with zero-inflated models for single-cell differential gene expression analysis. 

The analysis is built on a simulation-benchmark provided by *Van-den-Berge, 2019, Genome Biology*. 

**Aim: fast and reliable differential gene expression analysis method for scRNA-seq data**

## State-of-the-art: 
* edgeR with sample weights factored in has repeatedly been proposed as the best 
  method for differential gene expression analysis. However, it is too slow 
  slow with many samples (several hours with 20k cells, does not scale to 100k+ cells.)
* MAST and limma-voom are scaleable, but have been shown to no always perform best. 
* diffxpy implements differential gene expression analysis in python/tensorflow. 
  It uses the same fundamental model as edgeR but can run blazing-fast (on GPU)
  
  
Ideally, diffxpy could be extended to account for zero-inflation to achieve comparable performance as edgeR. 


## What I tried: 
* vanilla diffxpy: works fast and as well as vanilla edgeR
* scqtl: implements zero inflated negative binomial model in tensorflow. 
  not very fast and I failed to obtain the right model parameters to
  conduct a wald test. 
* fit a ZINB model from python statsmodels. Significantly better than 
  NB model, but FDR inflation
* run diffxpy with number of detected genes as co-factor. Failed because
  of matrix not being full rank (weird, as the same model works in edgeR)
* run diffxpy with ZINB weights from zinber/zinb-WAVE: diffxpy currently does 
  not support sample weights. Also zinb-WAVE/zingeR are quite slow and don't scale 
  to the 100k+ cells. 
  
## Next steps: 
* check the performance of `edgeR` with `#detected genes` as covariate. 
* contact diffxpy developers
   - why does it complain about the matrix not being full-rank (even though it works with edgeR)
   - what's their take on ZINB-models? 
* Why does the FDR-inflation occor with the statsmodel ZINB model? 