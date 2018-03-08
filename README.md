# RNA-Seq-data-analysis

The analysis of RNA Seq data is quite straightforward. Bioconductor is a popular website providing tools for the analysis of various high-throughput data.  edgeR is a popular package that fits generalize linear model to analyze RNA seq data (count data).  Below is a simple Bioconductor article to show this with examples and codes. For our data, we start the analyses from the step 5.2.  Our goal is to look for the differentially expressed genes across 2 clusters or 3 clusters identified by the imaging features.  
http://www.bioconductor.org/packages/devel/bioc/vignettes/systemPipeR/inst/doc/systemPipeRNAseq.pdf

For the reference, below is another thorough Bioconductor article documenting all the details of analyzing RNA seq data starting from scratch.  The package is DESeq2 (another popular tool).     
http://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html


1. Differential gene expression analysis based on the negative binomial distribution
Estimate variance-mean dependence in count data from high-throughput sequencing assays (over 50,000 mRNA) and test for differential expression based on a model using the negative binomial distribution. 

2. Building and running automated end-to-end analysis workflows for a wide range of next generation sequence (NGS) applications such as RNA-Seq.

3. Obtained gene-to-GO mappings from Human Genome Assembly GRCh38.p110 and did Batch GO term enrichment analysis


Here are some result figures.

![alt text](https://github.com/yue0814/RNA-Seq-data-analysis/blob/master/DEG2.PNG)

![alt text](https://github.com/yue0814/RNA-Seq-data-analysis/blob/master/DEG3.PNG)

![alt text](https://github.com/yue0814/RNA-Seq-data-analysis/blob/master/GOslimbarplot.PNG)

![alt text](https://github.com/yue0814/RNA-Seq-data-analysis/blob/masterGOslimbarplot2.PNG)

![alt text](https://github.com/yue0814/RNA-Seq-data-analysis/blob/master/sample_tree.PNG)
