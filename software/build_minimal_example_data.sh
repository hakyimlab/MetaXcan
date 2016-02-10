#!/usr/bin/env bash

bash build_minimal_data_requirements.sh

if [[ ! -f data/GWAS.tar.gz ]]; then
    wget https://s3.amazonaws.com/imlab-open/Data/MetaXcan/1000G-WB/data/GWAS.tar.gz -O data/GWAS.tar.gz
    cd data
    tar -xzvpf GWAS.tar.gz
    cd ..
fi

