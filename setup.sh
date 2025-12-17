#!/bin/bash

source ./env.sh

# check dorado

# check bedtools

# check featureCounts

git clone https://github.com/levend1s/rqc.git ${RQC_DIR}
cd ${RQC_DIR}
python3 -m virtualenv env
source env/bin/activate
pip install -r requirements.txt
deactivate

# get ONT RNA-seq


# get GLORI-seq data
