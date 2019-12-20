---
output:
  html_document: default
  pdf_document: default
---
# scClassifier

**scClassifier** is an R package designed for the simultaneous identification of common and rare cell types using single cell RNA-sequencing data. 

## Installation

* scClassifier depends on Seurat (version >=2.2.0 and <= 2.3.4). Users can retrieve the old version of Seurat from [CRAN](https://cran.r-project.org/src/contrib/Archive/Seurat/), and then run the following instruction to install Seurat

```
install.packages("/PATH-TO-SEURAT/Seurat_<VERSION>.tar.gz")
```

* When Seurat has been installed, run the following instructions or the R script [scClassifier-Install.R](scClassifier-Install.R)
```
git clone https://github.com/homopolymer/scClassifier.git
cd scClassifier
Rscript --vanilla scClassifier-Install.R
```


* To get started, check out the Quick Start Tutorial: ([View](Tutorial/scClassifier-Human-PBMC.pdf)) ([Download](Tutorial/scClassifier-Human-PBMC.Rmd))


#### 2019/12/16:
Version 0.1.0 released! 


	
