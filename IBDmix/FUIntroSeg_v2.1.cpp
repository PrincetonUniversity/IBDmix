#include <unistd.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef unsigned long int ulnt;

const double minesp = 1e-200;
ulnt **ibd;
double *slod;

// pa: reference allele frequency; pb: alternative allele frequency
double cal_lod(int source_gt, int target_gt, double pb, double aerr, double merr){
    double pa = 1-pb;
    double berr = 1 - aerr;
    double denom = 1 - aerr * berr;
    double term1 = berr * (1 - merr) + aerr * merr;
    double term2 = berr * merr + aerr * (1 - merr);
    if(source_gt == 0 && target_gt == 0)
        return log10(term1 / pa / denom);
    if(source_gt == 2 && target_gt == 2)
        return log10(term1 / pb / denom);

    if(source_gt == 0 && target_gt == 1)
        return log10(0.5 / denom * (term2 / pb + term1 / pa));
    if(source_gt == 2 && target_gt == 1)
        return log10(0.5 / denom * (term1 / pb + term2 / pa));

    if(source_gt == 1 && target_gt == 1)
        return -log10(2*pa*pb*(3-2*denom));

    if(source_gt == 1 && target_gt == 0)
        return -log10(pa * (3 - 2 * denom));
    if(source_gt == 1 && target_gt == 2)
        return -log10(pb * (3 - 2 * denom));

    if(source_gt == 0 && target_gt == 2)
        return log10((term2 + minesp) / pb / (denom + minesp));
    if(source_gt == 2 && target_gt == 0)
        return log10((term2 + minesp) / pa / (denom + minesp));

    return 0;
}

// The sum of the IBD LOD scores for variants in the interval is greater than a specificed value,
// and for which the sum cannot be increased by expanding the interval
int call_segment(double *lod, ulnt Lnum, double threlod){
    int Ibdnum = 0;
    double cumsum;
    
    ibd[Ibdnum][0] = 0; ibd[Ibdnum][1] = 0; slod[Ibdnum] = cumsum = 0;
    
    for(ulnt l = 0; l < Lnum; l++){
        cumsum += lod[l];
        if(cumsum >= slod[Ibdnum]){
            slod[Ibdnum] = cumsum;
            ibd[Ibdnum][1] = l;
        }
        if(cumsum < 0){
            if(slod[Ibdnum] >= threlod){
                l = ibd[Ibdnum][1];
                Ibdnum++;
            }
            ibd[Ibdnum][0] = ibd[Ibdnum][1] = (l+1);
            cumsum = slod[Ibdnum] = 0;
        }
    }
    if(slod[Ibdnum] >= threlod) Ibdnum++;
    return Ibdnum;  
}

double calLOD_segment(double *lod, ulnt start, ulnt end){
    double cumsum = 0;
    for(ulnt l = start; l <= end; l++)
        cumsum += lod[l];
    return cumsum;
}

int main_ibd(int argc, char *argv[]){
    ulnt Lnum,*pos,posi[2],prel;
    int c, chr, chro, Anum, Indnum, Inum, eInum, IBDnum, *sel, *indsel, **mgt, **agt, mafcut=1;
    double acerr=0.01, mderr=0.0025, mderrprop=2, *obsf, *err, **lod, threlod=3, tp, maxlod;
    char str[1024],**id;
    FILE *inp[3],*out;

    static const struct option long_opt[] = {
        {"help", no_argument, 0, 'h'},
        {"geno", required_argument, 0, 'g'},
        {"sample", required_argument, 0, 's'},
        {"output", required_argument, 0, 'o'},
        {"numAcnts", required_argument, 0, 'i'},
        {"numMdns", required_argument, 0, 'n'},
        {"numLoci", required_argument, 0, 'l'},
        {"LODcutoff", required_argument, 0, 'd'},
        {"MAFcutoff", required_argument, 0, 'm'},
        {"AHErr", required_argument, 0, 'a'},
        {"MHErrMax", required_argument, 0, 'e'},
        {"MHErrProp", required_argument, 0, 'c'},
        {"maskRegion", required_argument, 0, 'r'},
        {0, 0, 0, 0}
    };
    
    inp[2] = NULL;
    while ((c = getopt_long(argc, argv, "hg:s:o:i:n:l:d:m:a:e:c:r:", long_opt, NULL)) != -1) {
        switch (c) {
            case 'h':
                printf("Usage: %s -g <genotype file> -s <sample file> -o <output file> -n <#Individuals> -l <#Loci> [OPTIONS]\n",argv[0]);
                printf("-h, --help:         print this help and exit\n");
                printf("-g, --geno:         input file for genotype\n");
                printf("-s, --sample:       sample file\n");
                printf("-o, --output:       output file\n");
                printf("-i, --numAcnts:     number of archaic homoids\n");
                printf("-n, --numMdns:      number of modern humans\n");
                printf("-l, --numLoci:      number of total loci\n");
                printf("-d, --LODcutoff:    Log(odds) cutoff for introgressed segment calling (default: 3.0)\n");
                printf("-m, --MAFcutoff:    minor allele count cutoff for inclusive sites (default: 1)\n");
                printf("-a, --AHErr:        allele error rate for archaic DNA (default: 0.01)\n");
                printf("-e, --MHErrMax:     Maximum allele error rate for modern humans (default: 0.0025)\n");
                printf("-c, --MHErrProp:    Ratio between allele error rate and minor allele frequency (default: 2)\n");
                printf("-r, --maskRegion:   bed file for masked regions (optinal)\n");
                return(0);
            case 'g': inp[0] = fopen(optarg, "rt"); break;
            case 's': inp[1] = fopen(optarg, "rt"); break;
            case 'o': out = fopen(optarg,"wt"); break;
            case 'i': Anum = atoi(optarg); break;
            case 'n': Indnum = atoi(optarg); break;
            case 'l': Lnum = atol(optarg); break;
            case 'd': threlod = atof(optarg); break;
            case 'm': mafcut = atoi(optarg); break;
            case 'a': acerr = atof(optarg); break;
            case 'e': mderr = atof(optarg); break;
            case 'c': mderrprop = atof(optarg); break;
            case 'r': inp[2] = fopen(optarg, "rt"); break;
        }
    }
    
    pos = new ulnt[Lnum]; sel = new int[Lnum];
    indsel = new int[Indnum]; id = new char *[Indnum];
    obsf = new double[Lnum]; err = new double[Lnum];
    lod = new double *[Anum]; ibd = new ulnt *[Lnum]; slod = new double[Lnum];
    for(int i = 0; i < Indnum; i++)
        id[i] = new char[100];
    for(ulnt l = 0; l < Lnum; l++)
        ibd[l] = new ulnt[2];
    for(int a = 0; a < Anum; a++)
        lod[a] = new double[Lnum];
    
    // Input data
    Inum = 0;
    for(int i = 0; i < Indnum; i++){
        fscanf(inp[1],"%s%d",id[Inum], &indsel[i]);
        if(indsel[i]==1)
            Inum++;
    }
    fclose(inp[1]);
    
    agt = new int *[Lnum]; mgt = new int *[Lnum];
    for(ulnt l = 0; l < Lnum; l++){
        agt[l] = new int[Anum];
        mgt[l]=new int[Inum];
    }
    
    for(ulnt l = 0; l < Lnum; l++){
        sel[l] = 1;
        fscanf(inp[0],"%d%lu%*s%*s", &chr, &pos[l]);

        for(int a = 0; a < Anum; a++){
            fscanf(inp[0],"%d", &agt[l][a]);
            if(agt[l][a]==9)
                sel[l] = 0;
        }

        Inum = 0;
        for(int i = 0; i < Indnum; i++){
            fscanf(inp[0],"%d",&c);
            if(indsel[i]==1){
                mgt[l][Inum] = c;
                Inum++;
            }
        }
    }
    fclose(inp[0]);
    
    // Filtering based on minor allele count, and set error rate
    for(ulnt l = 0; l < Lnum; l++){
        eInum = 0; obsf[l] = 0;
        for(int i = 0; i < Inum; i++)
            if(mgt[l][i] != 9){
                eInum++;
                obsf[l] += mgt[l][i];
            }
        if(obsf[l]<=mafcut || (2*eInum-obsf[l]) <= mafcut)
            sel[l] = 0;
        if(eInum > 0)
            obsf[l] = 0.5*obsf[l]/eInum;
        else
            obsf[l] = 0;
        
        err[l] = mderr;
        if(err[l] > obsf[l]*mderrprop)
            err[l] = obsf[l]*mderrprop;
        else{
            if(err[l] > (1-obsf[l])*mderrprop)
                err[l] = (1-obsf[l])*mderrprop;
        }
    }
                       
    // Filtering based on mask files
    if(inp[2] != NULL){
        prel = 0;
        while(!feof(inp[2])){
            fscanf(inp[2],"%d%lu%lu\n",&chro, &posi[0], &posi[1]);
            for(ulnt l = prel; l < Lnum; l++)
                if(sel[l] == 1){
                    if(chr == chro && pos[l] > posi[0] && pos[l] <= posi[1]){
                        sel[l] = 0;
                        prel = l;
                    }
                    if(chro > chr || (chr == chro && pos[l] > posi[1]))
                        break;
                }
        }
        fclose(inp[2]);
    }
    
    // Introgression segment detection
    for(int i = 0; i < Inum; i++){
        // calculate log_odds = log10(Po(.|I)/Po(.|nI)) at each site for testing population
        for(ulnt l = 0; l < Lnum; l++)
            for(int a = 0; a < Anum; a++){
                lod[a][l] = 0;
                if(sel[l] == 0 && ((agt[l][a] == 0 && mgt[l][i] == 2) || (agt[l][a] == 2 && mgt[l][i] == 0)))
                    lod[a][l] = cal_lod(agt[l][a], mgt[l][i], obsf[l], acerr, err[l]);
                if(sel[l] == 1)
                    lod[a][l] = cal_lod(agt[l][a], mgt[l][i], obsf[l], acerr, err[l]);
            }

        for(int a = 0; a < Anum; a++){
            // Segment calling
            IBDnum = call_segment(lod[a], Lnum, threlod);

            for(int d = 0; d < IBDnum; d++){
                // introgressed segment: [start, end)
                if(ibd[d][1] < (Lnum - 1))
                    fprintf(out,"%s\t%d\t%lu\t%lu\t",id[i],chr,pos[ibd[d][0]],pos[ibd[d][1]+1]);
                else
                    fprintf(out,"%s\t%d\t%lu\t%lu\t",id[i],chr,pos[ibd[d][0]],pos[Lnum-1]);

                maxlod = slod[d];
                for(int aa = 0; aa < Anum; aa++){
                    if(aa == a)
                        fprintf(out,"%g\t", slod[d]);
                    else{
                        tp = calLOD_segment(lod[aa], ibd[d][0], ibd[d][1]);
                        fprintf(out,"%g\t", tp);
                        if(tp > maxlod){
                            maxlod = tp;
                        }
                    }
                }
                fprintf(out,"%g\n", maxlod);
            }
        }
    }
    fclose(out);
    
    for(ulnt l = 0; l < Lnum; l++){
        delete agt[l]; delete mgt[l]; delete ibd[l];
    }
    for(int i = 0; i < Indnum; i++)
        delete id[i];
    for(int a = 0; a < Anum; a++)
        delete lod[a];
    delete []pos; delete []obsf; delete []err; delete []agt; delete []mgt;
    delete []sel; delete []indsel; delete []id; delete []lod; delete []slod;
    delete []ibd;
    return 0;
}

