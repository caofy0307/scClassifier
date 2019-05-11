library(Seurat)
library(SingleCellExperiment)
library(scClassifier)
library(dplyr)
library(mclust)
library(ggpubr)


# human purified PBMCs
# load("/home/zengbio/Project/SingleCell/scClassifier/test/purified_PBMC_data.Robj")
load("/Users/zengbio/Project/SingleCell/SingleCellClassifier_DB/purified_PBMC_data.Robj")
use.cells <- purified_PBMC_obj@cell.names[(!purified_PBMC_obj@meta.data$type %in% c("Jurkat", "Stem")) & (!is.na(purified_PBMC_obj@meta.data$type))]
purified_PBMC_obj0 <- SubsetData(purified_PBMC_obj, cells.use = use.cells)

# label normalization
use.levels <- c("B cell", "T cell", "Monocyte", "Nature killer cell")
labels <- plyr::mapvalues(purified_PBMC_obj0@meta.data$type,
                          c("B", "Cytotoxic T", "Helper T", "Memory T", "Naive Cytotoxic T", "Naive T", "NK", "Regulatory T"),
                          c("B cell", "T cell", "T cell", "T cell", "T cell", "T cell", "Nature killer cell", "T cell"))
names(labels) <- purified_PBMC_obj0@meta.data$type
labels <- factor(labels, levels = use.levels)
purified_PBMC_obj0@meta.data$orig.type <- purified_PBMC_obj0@meta.data$type
purified_PBMC_obj0@meta.data$type <- labels


# human immune cell examples
# load("/home/zengbio/Project/SingleCell/scClassifier/test/immuno_navigator_human_expression.Robj")
load("/Users/zengbio/Project/SingleCell/SingleCellClassifier_DB/immuno_navigator_human_expression.Robj")
bulk.genemap <- fData(immuno_navigator_human_expression)
bulk.cellinfo <- pData(immuno_navigator_human_expression)
bulk.profile <- exprs(immuno_navigator_human_expression)
bulk.profile[bulk.profile<0] <- 0
bulk.profile <- exp(bulk.profile) - 1
bulk.profile <- bulk.profile[, bulk.cellinfo$type %in% c("T cell", "B cell", "Monocyte", "Nature killer cell")]
bulk.cellinfo <- bulk.cellinfo[bulk.cellinfo$type %in% c("T cell", "B cell", "Monocyte", "Nature killer cell"),]


# classification experiment
results <- data.frame(Iter = numeric(0), N = numeric(0), ARI = double(0), Method = character(0), stringsAsFactors = FALSE)

niter <- 1
train.n <- 2000
use.cells.old <- c()
iter <- 1

use.cells <- sample(setdiff(purified_PBMC_obj0@cell.names, use.cells.old), train.n)
purified_PBMC_obj <- SubsetData(purified_PBMC_obj0, cells.use = use.cells)

genes.use = SeuratDEGenes(purified_PBMC_obj@raw.data[, purified_PBMC_obj@cell.names], min.pct = 0.2, thresh.use = 0.2)
# genes.use = SeuratHVGenes(purified_PBMC_obj@raw.data[, purified_PBMC_obj@cell.names], x.high.cutoff = 5)
# genes.use = SingleRSDGenes(bulk.profile, bulk.cellinfo$type, size = 500)

# simulation
pseudo.catalogue = PseudoSingleCellGenerator(as.matrix(purified_PBMC_obj@raw.data[, purified_PBMC_obj@cell.names]),
                                             bulk.profile,
                                             bulk.cellinfo$type,
                                             genes.use = NULL,
                                             lambda = 0.01,
                                             sizes = 1000,
                                             BPPARAM = SnowParam())

y = SVMClassifier(as.matrix(purified_PBMC_obj@raw.data[, purified_PBMC_obj@cell.names]),
                  pseudo.catalogue$scdata,
                  pseudo.catalogue$celltype,
                  genes.use = genes.use,
                  do.scale = F,
                  importance.thresh = 100,
                  cross = 5)
table(y$celltype, purified_PBMC_obj@meta.data$type)
sum(y$celltype == purified_PBMC_obj@meta.data$type) / train.n


y4 = XGBoostClassifier(as.matrix(purified_PBMC_obj@raw.data[, purified_PBMC_obj@cell.names]),
                       pseudo.catalogue$scdata,
                       pseudo.catalogue$celltype,
                       eta = 1,
                       max.depth = 2,
                       nrounds = 1000,
                       booster = "gbtree",
                       genes.use = genes.use,
                       importance.thresh = 100)
table(y4$celltype, purified_PBMC_obj@meta.data$type)
sum(y4$celltype == purified_PBMC_obj@meta.data$type) / train.n



# multinomial simulation
pseudo.catalogue2 = PseudoSingleCellGenerator(as.matrix(purified_PBMC_obj@raw.data[, purified_PBMC_obj@cell.names]),
                                              bulk.profile,
                                              bulk.cellinfo$type,
                                              simulator = "Multinomial",
                                              genes.use = NULL,
                                              sizes = 1000,
                                              BPPARAM = SerialParam())

y2 = SVMClassifier(as.matrix(purified_PBMC_obj@raw.data[, purified_PBMC_obj@cell.names]),
                   pseudo.catalogue2$scdata,
                   pseudo.catalogue2$celltype,
                   genes.use = genes.use,
                   do.scale = F,
                   importance.thresh = 100,
                   cross = 5)
table(y2$celltype, purified_PBMC_obj@meta.data$type)
sum(y2$celltype == purified_PBMC_obj@meta.data$type) / train.n


y3 = XGBoostClassifier(as.matrix(purified_PBMC_obj@raw.data[, purified_PBMC_obj@cell.names]),
                       pseudo.catalogue2$scdata,
                       pseudo.catalogue2$celltype,
                       eta = 1,
                       max.depth = 2,
                       nrounds = 1000,
                       booster = "gbtree",
                       genes.use = genes.use,
                       importance.thresh = 100)
table(y3$celltype, purified_PBMC_obj@meta.data$type)
sum(y3$celltype == purified_PBMC_obj@meta.data$type) / train.n

