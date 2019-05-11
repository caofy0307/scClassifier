# Using scRNA-seq as reference
#
#
library(Seurat)
library(SingleCellExperiment)
library(scClassifier)
library(dplyr)
library(mclust)
library(ggpubr)

dir <- "/Users/zengbio/Documents/workspace/2018/classification/experiment/evaluation/10x_Genomics/bulk"

# data.dir <- "/Users/zengbio/Project/SingleCell/SingleCellClassifier_Experiments/Human_Purified_PBMC"
# load(file.path(data.dir,"results_20180518.rda"))

# human purified PBMCs
load("/Users/zengbio/Project/SingleCell/SingleCellClassifier_DB/purified_PBMC_data.Robj")
use.cells <- purified_PBMC_obj@cell.names[!purified_PBMC_obj@meta.data$type %in% "Jurkat"]
purified_PBMC_obj0 <- SubsetData(purified_PBMC_obj, cells.use = use.cells)

# human immune cell examples
load("/Users/zengbio/Project/SingleCell/SingleCellClassifier_DB/immuno_navigator_human_expression.Robj")
bulk.genemap <- fData(immuno_navigator_human_expression)
bulk.cellinfo <- pData(immuno_navigator_human_expression)
bulk.profile <- exprs(immuno_navigator_human_expression)
bulk.profile[bulk.profile<0] <- 0

results <- data.frame(Iter = numeric(0), N = numeric(0), ARI = double(0), Method = character(0), stringsAsFactors = FALSE)

niter <- 10
train.n <- 1000
use.cells.old <- c()
iter <- 1
  use.cells <- sample(setdiff(purified_PBMC_obj0@cell.names, use.cells.old), train.n)
  purified_PBMC_obj <- SubsetData(purified_PBMC_obj0, cells.use = use.cells)

  # filtering
  purified_PBMC_obj <- FilterCells(purified_PBMC_obj, subset.names = c("nGene"), low.thresholds = c(90), high.thresholds = c(Inf))
  # Normalization
  purified_PBMC_obj <- NormalizeData(purified_PBMC_obj)
  # find HVG
  purified_PBMC_obj <- FindVariableGenes(purified_PBMC_obj, x.low.cutoff = 0.1, x.high.cutoff = 3, y.cutoff = 0.15)
  var.genes <- purified_PBMC_obj@var.genes
  # scaling
  purified_PBMC_obj <- ScaleData(purified_PBMC_obj)

  ### Perform PCA and JackStraw
  # PCA
  purified_PBMC_obj <- RunPCA(purified_PBMC_obj, do.print = F)
  # JackStraw
  purified_PBMC_obj <- JackStraw(purified_PBMC_obj, num.replicate = 100, num.cores = 7, do.par = TRUE)
  # # projection
  purified_PBMC_obj <- ProjectPCA(purified_PBMC_obj, do.print = F)
  # genes that are significantly correlated with PCs
  sig.genes <- PCASigGenes(purified_PBMC_obj, pcs.use = 1:5, pval.cut = 1e-10, max.per.pc = 200)
  # Run PCA again but use significant genes
  purified_PBMC_obj <- RunPCA(purified_PBMC_obj, pc.genes = sig.genes, do.print = F)


  # unsupervised clustering
  purified_PBMC_obj <- FindClusters(object = purified_PBMC_obj,
                                    reduction.type = "pca",
                                    dims.use = 1:10,
                                    k.param = 30,
                                    print.output = 0,
                                    save.SNN = FALSE,
                                    force.recalc = TRUE)

  # cluster-specific markers
  markers <- FindAllMarkers(object = purified_PBMC_obj,
                            genes.use = unique(c(var.genes, sig.genes)),
                            only.pos = TRUE,
                            min.pct = 0.25,
                            thresh.use = 0.25)
  markers %>% group_by(cluster) %>% top_n(30, avg_logFC) -> topMarkers
  # topMarkers <- topMarkers[topMarkers$p_val_adj<0.1,]
  topGenes <- unique(topMarkers$gene)
  de.genes <- topGenes

  # label normalization
  levels <- c("B cell", "T cell", "Monocyte", "Nature killer cell", "Stem cell", "Others")
  labels <- plyr::mapvalues(purified_PBMC_obj@meta.data$type,
                            c("B", "Cytotoxic T", "Helper T", "Memory T", "Naive Cytotoxic T", "Naive T", "NK", "Regulatory T", "Stem"),
                            c("B cell", "T cell", "T cell", "T cell", "T cell", "T cell", "Nature killer cell", "T cell", "Stem cell"))
  names(labels) <- purified_PBMC_obj@meta.data$type
  labels <- factor(labels, levels = levels)
  purified_PBMC_obj@meta.data$orig.type <- purified_PBMC_obj@meta.data$type
  purified_PBMC_obj@meta.data$type <- labels



  # cell type catalogue
  catalogue.genes <- rownames(bulk.profile)
  genes.use <- intersect(de.genes, catalogue.genes)
  catalogue.profile <- bulk.profile[catalogue.genes,]
  catalogue.celltype <- bulk.cellinfo$type
  sc.profile <- as.matrix(purified_PBMC_obj@data[genes.use, ])


  z <- DirichletMultinomialClassifier(sc.profile,
                                      catalogue.profile,
                                      catalogue.celltype,
                                      do.training = T,
                                      genes.use = genes.use)
  k <- z$celltype
  scc.ari <- adjustedRandIndex(k, purified_PBMC_obj@meta.data$type)
  print(paste(iter, train.n, scc.ari, "DirichletMultinomialClassifier"))
