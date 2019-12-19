# Load Bioconductor
message("Loading Bioconductor to install required packages.")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install Bioconductor requirements: BiocParallel
message("Installing required packages from Bioconductor.")
if (!requireNamespace("BiocParallel", quietly = TRUE))
  BiocManager::install('BiocParallel', suppressUpdates=T)

# Install CRAN requirements: rsvd, expm, MCMCprecision
message("Installing required packages from CRAN")
if (!requireNamespace("rsvd", quietly = TRUE))
  BiocManager::install('rsvd', suppressUpdates=T)
if (!requireNamespace("expm", quietly = TRUE))
  BiocManager::install('expm', suppressUpdates=T)
if (!requireNamespace("MCMCprecision", quietly = TRUE))
  BiocManager::install("MCMCprecision", suppressUpdates=T)

# # Install GitHub requirements: SingleR
# message("Installing required packages from GitHub")
# if (!requireNamespace("SingleR", quietly = TRUE))
#   BiocManager::install("dviraran/SingleR", suppressUpdates=T)


# Install scClassifier (and devtools, if necessary)
# BiocManager::install("homopolymer/scClassifier", suppressUpdates=T)
install.packages(".", repos = NULL, type = "source")

# Check that scClassifier installed.
if (requireNamespace("scClassifier", quietly = TRUE)) {
  message("scClassifier installed successfully!")
  message('You can load it by typing: library("scClassifier")')
  message('Try "?scClassifier" for starting tips.')
} else {
  message("Something went wrong. It doesn't seem that scClassifier installed correctly.")
}
