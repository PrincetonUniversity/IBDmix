# IBDmix

1. For empirical data, first step is to merge original archaic vcf and modern human vcf
with merge.sh, which depends on 01\_GenotypeInputFromVCF.cpp

2. Second step use IBDmix to get the original calls.
There are different bash codes for calling sequence for just one population vs
multiple populations (IBDmix works as calling one population at once but
for multiple populations we can make multiple sample files and submit one
array job) and there are here:
FUIntroSegcalling\_onePop.sh
FUIntroSegcalling\_wholePop.sh
Which depends on 
FUIntroSeg_v2.1.cpp;

3. Merge the original calls for overlapping sequences from one individual.
FUIntroSeg_Summary.py
