#!/usr/bin/env bash

#aws s3 cp s3://imlab-open/Data/MetaXcan/1000G-WB/ . --recursive
if [[ ! -d "data" ]] ; then
    mkdir data
fi

if [ ! -f data/DGN-WB_0.5.db ]; then
    wget https://s3.amazonaws.com/imlab-open/Data/MetaXcan/1000G-WB/data/DGN-WB_0.5.db -O data/DGN-WB_0.5.db
fi

if [[ ! -d intermediate ]] ; then
    mkdir intermediate
fi

if [[ ! -d intermediate/cov ]] ; then
    mkdir intermediate/cov
fi

if [[ ! -f intermediate/cov/covariance.txt.gz ]]; then
    wget https://s3.amazonaws.com/imlab-open/Data/MetaXcan/1000G-WB/intermediate/covariance.txt.gz -O intermediate/cov/covariance.txt.gz
fi

