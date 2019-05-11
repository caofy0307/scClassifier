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

niter <- 1
train.n <- 500
use.cells.old <- c()
iter <- 1

  use.cells <- sample(setdiff(purified_PBMC_obj0@cell.names, use.cells.old), train.n)
  purified_PBMC_obj <- SubsetData(purified_PBMC_obj0, cells.use = use.cells)

  # label normalization
  use.levels <- c("B cell", "T cell", "Monocyte", "Nature killer cell", "Stem cell", "Others")
  labels <- plyr::mapvalues(purified_PBMC_obj@meta.data$type,
                            c("B", "Cytotoxic T", "Helper T", "Memory T", "Naive Cytotoxic T", "Naive T", "NK", "Regulatory T", "Stem"),
                            c("B cell", "T cell", "T cell", "T cell", "T cell", "T cell", "Nature killer cell", "T cell", "Stem cell"))
  names(labels) <- purified_PBMC_obj@meta.data$type
  labels <- factor(labels, levels = use.levels)
  purified_PBMC_obj@meta.data$orig.type <- purified_PBMC_obj@meta.data$type
  purified_PBMC_obj@meta.data$type <- labels



  # cell type catalogue
  # genes.use <- SingleRSDGenes(bulk.profile, bulk.cellinfo$type, size = 300)
  genes.use <- SeuratDEGenes(as.matrix(purified_PBMC_obj@data), resolution = 1.2)
  genes.use <- intersect(genes.use, intersect(rownames(purified_PBMC_obj@data),
                                              rownames(bulk.profile)))
  catalogue.profile <- bulk.profile[genes.use, ]
  catalogue.celltype <- bulk.cellinfo$type
  sc.profile <- as.matrix(purified_PBMC_obj@data[genes.use, ])

  # z <- apply(sc.profile,
  #            2,
  #            JaitinClassifier,
  #            catalogue.profile = catalogue.profile,
  #            catalogue.celltype = catalogue.celltype)
  # k <- lapply(z, function(x) x$celltype) %>% unlist
  #
  # scc.ari <- adjustedRandIndex(k, purified_PBMC_obj@meta.data$type)
  # print(paste(iter, train.n, scc.ari, "JaitinClassifier"))


  # z <- MultinomialClassifier(sc.profile,
  #                            catalogue.profile,
  #                            catalogue.celltype,
  #                            genes.use = genes.use)
  # k <- z$celltype
  # scc.ari <- adjustedRandIndex(k, purified_PBMC_obj@meta.data$type)
  # print(paste(iter, train.n, scc.ari, "MultinomialClassifier"))
  #
  z <- DirichletMultinomialClassifier(sc.profile,
                                      catalogue.profile,
                                      catalogue.celltype,
                                      # do.training = F,
                                      genes.use = genes.use)
  k <- z$celltype
  scc.ari <- adjustedRandIndex(k, purified_PBMC_obj@meta.data$type)
  print(paste(iter, train.n, scc.ari, "DirichletMultinomialClassifier"))

  z2 <- DiffusionKernelAveraging(sc.profile, t(z$score.matrix),
                                 genes.use = genes.use,
                                 diffuse.thresh = NULL,
                                 gamma = 0.01,
                                 tau = 5,
                                 do.pca = F,
                                 pca.params = 30)
  k2 <- z2$celltype
  scc.ari <- adjustedRandIndex(k2, purified_PBMC_obj@meta.data$type)
  print(paste(iter, train.n, scc.ari, "DirichletMultinomialClassifier"))

  z2 <- RandomWalkAveraging(sc.profile, t(z$score.matrix),
                            genes.use = genes.use,
                            diffuse.thresh = NULL,
                            tau = Inf,
                            alpha = 0.4,
                            do.pca = F,
                            pca.params = 30)
  k2 <- z2$celltype
  scc.ari <- adjustedRandIndex(k2, purified_PBMC_obj@meta.data$type)
  print(paste(iter, train.n, scc.ari, "DirichletMultinomialClassifier"))

  z2 <- DanaRandomWalk(sc.profile, t(z$score.matrix),
                       genes.use = genes.use,
                       # alpha = 0.4,
                       thresh = 0.25,
                       neighbor.size = 1,
                       do.pca = F,
                       pca.params = 30)
  k2 <- z2$celltype
  scc.ari <- adjustedRandIndex(k2, purified_PBMC_obj@meta.data$type)
  print(paste(iter, train.n, scc.ari, "DirichletMultinomialClassifier"))

  # z <- WeightedClassifier(sc.profile,
  #                         catalogue.profile,
  #                         catalogue.celltype,
  #                         latent.classifier = "MultinomialClassifier",
  #                         genes.use = genes.use,
  #                         do.pca = F)
  # k <- z$celltype
  # scc.ari <- adjustedRandIndex(k, purified_PBMC_obj@meta.data$type)
  # print(paste(iter, train.n, scc.ari, "WeightedClassifier"))
