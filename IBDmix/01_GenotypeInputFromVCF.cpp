#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

typedef unsigned long int ulnt;
// Compile: g++ 01_GenotypeInputFromVCF.cpp -o mergeVCF

// Step 1. Filtering out multi-allelic SNVs and indels from archaic genomes
// Single Archaic genomes: $vcftools --gzvcf /net/akey/vol2/rcmccoy/archaic_vcf/Altai/chr1_mq25_mapab100.vcf.gz --out Altai_chr1 --remove-indels --min-alleles 1 --max-alleles 2 --recode       ------to remove indels and multiallelic sites, only keep homozygous or biallelic sites
// Multiple Archaic genomes: $bcftools merge -m all --threads 2 ./Altai/chr1_mq25.vcf.gz ./Vindi/chr1_mq25.vcf.gz ./Denisova/chr1_mq25.vcf.gz | bcftools view --threads 2 --min-alleles 1 --max-alleles 2 -V indels -g ^miss -O z -o ./InterSectAVD_chr1_mq25.vcf.gz

// Step 2. Merge archaic vcf with vcf from 1000 genomes
// $gunzip InterSectAVD_chr1_mq25.vcf.gz
// $tabix -h ./ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 1:1-250000000 | ./mergeVCF -a ../Archaic/InterSectAVD_chr1_mq25.vcf -l 6260856 -d 3 -n 2504 -o InterSectAVD_1KGP3_chr1.txt
const int max=10000000;

int main(int argc, char *argv[]){
	ulnt Lnum, eLnum, *pos, posi, prel;
	int c,Indnum,Inum,chr,**gt, geno;
	char str[2048],***pol, **poly;
	FILE *inp, *out;
	
	while ((c = getopt(argc, argv, "a:l:d:n:o:")) >= 0) {
		switch (c) {
			case 'a': inp = fopen(optarg, "rt"); break; // Archaic vcf from step 1
			case 'l': Lnum = atol(optarg); break; // The number of sites for archaic vcf
			case 'd': Inum = atoi(optarg); break; // The number of archaic samples
			case 'n': Indnum = atoi(optarg); break; // The number of modern humans
			case 'o': out = fopen(optarg, "wt"); break; // The output file
		}
	}
	
	pos = new ulnt[Lnum]; gt = new int *[Lnum], pol=new char **[Lnum];
	poly = new char *[2]; poly[0] = new char[max]; poly[1] = new char[max];
	for(ulnt l = 0; l < Lnum; l++){pol[l] = new char *[2]; pol[l][0] = new char[2]; pol[l][1] = new char[2];gt[l] = new int[Inum];}
	
	while(fgetc(inp)=='#'){
		if(fgetc(inp)=='#') fgets(poly[0], max, inp);
		else{
			fscanf(inp,"%*s%*s%*s%*s%*s%*s%*s%*s%*s");
			for(int i = 0; i < Inum; i++) fscanf(inp,"%*s");
			break;
		}
	}
	eLnum = 0;
	for(ulnt l = 0; l < Lnum; l++){
		fscanf(inp,"%d%lu%*s%s%s%*s%*s%*s%*s", &chr, &pos[eLnum], poly[0], poly[1]); // fscanf reads each line, and runs the following commands, then comes back to scan next line and runs the following commands, and so on...
		if(strlen(poly[0])==1 && strlen(poly[1])==1){ 
			strcpy(pol[eLnum][0], poly[0]);
			strcpy(pol[eLnum][1], poly[1]);
			c = 0;
			for(int i = 0; i < Inum; i++){
				fscanf(inp,"%s",str);
				gt[eLnum][i] = 9;
				if(strstr(str,"0/0")!=NULL || strstr(str,"0|0")!=NULL) gt[eLnum][i] = 0;
				if(strstr(str,"0/1")!=NULL || strstr(str,"0|1")!=NULL || strstr(str,"1/0")!=NULL || strstr(str,"1|0")!=NULL ) gt[eLnum][i] = 1;
				if(strstr(str,"1/1")!=NULL || strstr(str,"1|1")!=NULL) gt[eLnum][i] = 2;
				if(gt[eLnum][i] != 9) c += 1;
			}
			if(c > 0) eLnum++;
		}
		else{
			for(int i = 0; i < Inum; i++) fscanf(inp,"%*s");
		}
	}
	fclose(inp);
	printf("%d\t%lu\n",chr, eLnum);
	
	while(fgetc(stdin)=='#'){
		if(fgetc(stdin)=='#') fgets(poly[0], max, stdin);
		else{
			fscanf(stdin,"%*s%*s%*s%*s%*s%*s%*s%*s%*s");
			for(int i = 0; i < Indnum; i++) fscanf(stdin,"%*s");
			break;
		}
	}
	
	prel=0;
	while(!feof(stdin)){
		fscanf(stdin,"%*s%lu%*s%s%s%*s%*s%*s%*s", &posi, poly[0], poly[1]);
		for(ulnt l = prel; l < eLnum; l++){
			if(pos[l] < posi){
				// Assuming all the sites in the bed file from modern humans is callable
				c = 0;
				for(int i = 0; i < Inum; i++) if(gt[l][i]!=9 && gt[l][i]!=0) c+=1;
				if(c > 0){
					fprintf(out,"%d\t%lu\t%s\t%s\t", chr, pos[l], pol[l][0], pol[l][1]);
					for(int i = 0; i < Inum; i++) fprintf(out,"%d\t", gt[l][i]);
					for(int i = 0; i < Indnum; i++) fprintf(out,"0\t");
					fprintf(out,"\n");
				}
				prel = l+1;
			}
			if(pos[l] == posi){
				if(strcmp(pol[l][0], poly[0]) == 0 && ((pol[l][1][0]=='.' && strlen(poly[1])==1) || strcmp(pol[l][1],poly[1])==0)){
					fprintf(out,"%d\t%lu\t%s\t%s\t", chr, pos[l], pol[l][0], poly[1]);
					for(int i = 0; i < Inum; i++) fprintf(out,"%d\t", gt[l][i]);
					for(int i = 0; i < Indnum; i++){
						fscanf(stdin,"%s",str);
						geno = 9;
						if(strstr(str,"0/0")!=NULL || strstr(str,"0|0")!=NULL) geno = 0;
						if(strstr(str,"0/1")!=NULL || strstr(str,"0|1")!=NULL || strstr(str,"1/0")!=NULL || strstr(str,"1|0")!=NULL ) geno = 1;
						if(strstr(str,"1/1")!=NULL || strstr(str,"1|1")!=NULL) geno = 2;
						fprintf(out,"%d\t",geno);
					}
					fprintf(out,"\n");
				}
				// Remove multi-allelic sites
				else{
					for(int i = 0; i < Indnum; i++) fscanf(stdin, "%*s");
				}
				prel = l+1;
				break;
			}
			// The sites in the ancient DNA is not callable
			if(pos[l] > posi){
				for(int i = 0; i < Indnum; i++) fscanf(stdin, "%*s");
				break;
			}
		}
	}
	fclose(out);
	
	for(ulnt l = 0; l < Lnum; l++){delete pol[l][0];delete pol[l][1];delete pol[l];delete gt[l];}
	delete poly[0]; delete poly[1]; delete []poly;delete []pol;delete []pos;delete []gt;
	return 0;
}

