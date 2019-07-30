#!/bin/bash

if [[ $# -lt 3 ]]; then
    echo "USAGE: $0 LENGTH_CUTOFF LOD_CUTOFF POPULATION [input] [output]"
    exit 1
fi

if [[ $# -eq 4 ]]; then
    echo "unable to determine input or output, specify stdin/stdout with '-'"
    exit 1
fi

length=$1
lod=$2
pop=$3
infile=/dev/stdin
outfile=/dev/stdout
if [[ $# -eq 5 ]]; then  # use inputs
    if [[ $4 != "-" ]]; then
        infile=$4
    fi
    if [[ $5 != "-" ]]; then
        outfile=$5
    fi
fi

declare -A pop2anc=(
    [CHB]=EAS [JPT]=EAS [CHS]=EAS [CDX]=EAS [KHV]=EAS
    [CEU]=EUR [TSI]=EUR [FIN]=EUR [GBR]=EUR [IBS]=EUR
    [YRI]=AFR [LWK]=AFR [GWD]=AFR [MSL]=AFR [ESN]=AFR [ASW]=AFR [ACB]=AFR
    [MXL]=AMR [PUR]=AMR [CLM]=AMR [PEL]=AMR
    [GIH]=SAS [PJL]=SAS [BEB]=SAS [STU]=SAS [ITU]=SAS
)

awk -v length_cutoff=$length -v lod_cutoff=$lod \
    -v pop=$pop -v anc=${pop2anc[$pop]} -v OFS='\t' '
NR == 1{
    if($1 == "ID"){
        print $0, "length", "pop", "anc"
        next
    }
    else{
        print "ID", "chr", "start", "end", "LOD", "length", "pop", "anc"
    }
}
{
    len = $4 - $3
    lod = $5
    if(len > length_cutoff && lod > lod_cutoff){
        print $0, len, pop, anc | "sort --key=1,1 --key=3n,3"
    }
}
' $infile > $outfile
