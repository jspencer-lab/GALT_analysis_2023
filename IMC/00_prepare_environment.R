Sys.unsetenv('GITHUB_PAT')

install.packages("BiocManager", repos='http://cran.us.r-project.org')

BiocManager::install(c("cytomapper", "imcRtools", "harmony", "randomForest", "cowplot", "bluster", "viridis",
                       "CATALYST", 'dittoSeq', 'mclust', 'scran'))
install.packages(c('patchwork', 'caret','tidyverse', 'Hmisc', 'devtools', 'dichromat'), repos='http://cran.us.r-project.org')
devtools::install_github("i-cyto/Rphenograph")

renv::snapshot()
