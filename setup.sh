#!/bin/bash

source ./env.sh

# check dorado
# https://github.com/nanoporetech/dorado

# get R
# https://posit.co/download/rstudio-desktop/

# get R packages
Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="http://cran.us.r-project.org")'
Rscript -e 'BiocManager::install("edgeR")'
Rscript -e 'BiocManager::install("GO.db")'

Rscript -e 'BiocManager::install("clusterProfiler")'


Rscript -e "install.packages('ggplot2', repos='https://cloud.r-project.org')"
Rscript -e "install.packages('ggridges', repos='https://cloud.r-project.org')"
Rscript -e "install.packages('pheatmap', repos='https://cloud.r-project.org')"
Rscript -e "install.packages('compositions', repos='https://cloud.r-project.org')"

Rscript -e "install.packages('devtools', repos='https://cloud.r-project.org')"
Rscript -e 'devtools::install_github("psyteachr/introdataviz")'

Rscript -e "install.packages('ggalluvial', repos='https://cloud.r-project.org')"

Rscript -e "install.packages('ggsignif', repos='https://cloud.r-project.org')"


Rscript -e "install.packages('effectsize', repos='https://cloud.r-project.org')"
Rscript -e "install.packages('eulerr', repos='https://cloud.r-project.org')"


# check modkit
# https://github.com/nanoporetech/modkit/releases/download/v0.6.0/modkit_v0.6.0_u16_x86_64.tar.gz
# cargo build
# ./target/debug/modkit


# brew install bedtools
# brew install samtools


# check featureCounts
# https://sourceforge.net/projects/subread/

# get rqc
git clone https://github.com/levend1s/rqc.git ${RQC_DIR}
cd ${RQC_DIR}
# pip3 install virtualenv
python3 -m virtualenv env
source env/bin/activate
pip3 install -r requirements.txt
deactivate

# get ONT RNA-seq


# get GLORI-seq data
