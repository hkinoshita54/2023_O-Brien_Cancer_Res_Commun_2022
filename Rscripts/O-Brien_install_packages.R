# install velocyto.R

# required

# install boost by brew in terminal

# install gfortran, see below
# https://github.com/fxcoudert/gfortran-for-macOS/releases/tag/11.2-bigsur-intel

# install gcc by brew in terminal -> failed, probably d/t gfortran installed above

# the below may be necessary
# https://medium.com/biosyntax/following-up-library-dependency-when-compiling-r-packages-89f191b9f227
# modified Makeconf in Library/Frameworks/R.framework/Resources/etc
# FLIBS = ...

BiocManager::install("pcaMethods") # required to install velocyto



# set github personal access token to avoid errors related to "API rate limit"
usethis::create_github_token()
gitcreds::gitcreds_set()
usethis::edit_r_environ() # GITHUB_PAT=XXXXXXX...(git hub personal access token)

# install velocyto
library(devtools)
install_github("velocyto-team/velocyto.R")

# install seurat disk
remotes::install_github("mojaveazure/seurat-disk")

# compile Seurat packages from souse to have previous versions
install.packages("~/Downloads/Seurat_4.9.9.9058.tar.gz", repos = NULL, type="source")
install.packages("~/Downloads/SeuratObject_4.9.9.9091.tar.gz", repos = NULL, type="source")

# CytoTrace
# https://cytotrace.stanford.edu/
BiocManager::install("sva") # this is dependency
devtools::install_local("~/Downloads/CytoTRACE_0.3.3.tar.gz")

# Monocle3
# https://cole-trapnell-lab.github.io/monocle3/docs/installation/
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))
devtools::install_github('cole-trapnell-lab/monocle3')
