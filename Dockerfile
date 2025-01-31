
## BIOCONDUCTOR
## Pull the Bioconductor v3.20 Docker image
FROM bioconductor/bioconductor_docker:RELEASE_3_20

## Change home directory
WORKDIR /home/rstudio/

## Install dependencies (2 passes)
RUN R -e 'BiocManager::install(c("BiocParallel", "BiocStyle", "Biostrings", "curl", "DEqMS", "dplyr", "GOfuncR", "QFeatures", "lme4", "msqrob2", "MSstatsTMT", "Mus.musculus", "pryr", "quarto", "readr", "stringr", "tidyverse", "cowplot", "ggVennDiagram", "ggridges", "gt", "flextable", "ggh4x", "gprofiler2", "UpSetR", "ggborderline", "viridisLite", "impute", "ExploreModelMatrix"))'
RUN R -e 'BiocManager::install(c("BiocParallel", "BiocStyle", "Biostrings", "curl", "DEqMS", "dplyr", "GOfuncR", "QFeatures", "lme4", "msqrob2", "MSstatsTMT", "Mus.musculus", "pryr", "quarto", "readr", "stringr", "tidyverse", "cowplot", "ggVennDiagram", "ggridges", "gt", "flextable", "ggh4x", "gprofiler2", "UpSetR", "ggborderline", "viridisLite", "impute", "ExploreModelMatrix"))'
## Install msTrawler (timestamped, 2 passes)
RUN R -e 'BiocManager::install("calico/msTrawler", ref = "078eb79")'
RUN R -e 'BiocManager::install("calico/msTrawler", ref = "078eb79")'

## Install Authors-block Extension For Quarto
RUN quarto add --no-prompt kapsner/authors-block

## Install tinytex through Quarto
RUN quarto install --no-prompt tinytex
