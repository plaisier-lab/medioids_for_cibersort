# Converting a single-cell expression dataset into medioids that can be used in CIBERSORT deconvolution analysis

Single-cell expression data can be used to identify discretized subpopultaions of cells from complex tissues. These discretized populations can then be used to provide exemplars for deconvoluting esitmates for these expression patterns from bulk expression data. This code was written to:

1. Compute the mediods (scanpy in Python)
2. Plot the results of the deconvolution (scanpy in Python)

## Environments for running the code
These code are designed to be run inside docker vitrual environments that can be downloaded through dockerhub:

```
docker pull saoconn1/scrnaseq_latest
```

## Running the code
### 1. Analyze the RNA-seq data and prepare it for mediod calculation (Seurat in R)
First fire up the docker image and link it to your files:

```
docker run -it -v '<path_to_your_files>:/files' saoconn1/scrnaseq_latest
```

Then 

## Results

