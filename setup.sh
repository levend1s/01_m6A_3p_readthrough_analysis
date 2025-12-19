#!/bin/bash

source ./env.sh

# check dorado
# https://github.com/nanoporetech/dorado

# get R
# https://posit.co/download/rstudio-desktop/

# check modkit
# https://github.com/nanoporetech/modkit/releases/download/v0.6.0/modkit_v0.6.0_u16_x86_64.tar.gz

# check bedtools
# https://formulae.brew.sh/formula/samtools

# check featureCounts
# https://sourceforge.net/projects/subread/

# get rqc
git clone https://github.com/levend1s/rqc.git ${RQC_DIR}
cd ${RQC_DIR}
# pip3 install virtualenv
python3 -m virtualenv env
source env/bin/activate
pip install -r requirements.txt
deactivate

# get ONT RNA-seq


# get GLORI-seq data
