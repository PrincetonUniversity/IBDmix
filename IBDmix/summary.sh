#!/bin/bash

if [[ $# -lt 2 ]]; then
    echo "USAGE: $0 LENGTH_CUTOFF LOD_CUTOFF [input] [output]"
    exit 1
fi

if [[ $# -eq 3 ]]; then
    echo "unable to determine input or output, specify stdin/stdout with '-'"
    exit 1
fi

length=$1
lod=$2
infile=/dev/stdin
outfile=/dev/stdout
if [[ $# -eq 4 ]]; then  # use inputs
    if [[ $3 != "-" ]]; then
        infile=$3
    fi
    if [[ $4 != "-" ]]; then
        outfile=$4
    fi
fi

awk -v length_cutoff=$length -v lod_cutoff=$lod -F '\t' '{
    if($1 == "ID"){
        print $0, "length"
    }
    else{
        len = $4 - $3
        lod = $6
        if(len > length_cutoff && lod > lod_cutoff){
            print $0, len
        }
    }
}' $infile | sort \
    --key=1,1 \
    --key=3n,3 \
    > $outfile
