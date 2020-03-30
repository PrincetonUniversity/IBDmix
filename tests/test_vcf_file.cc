#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <iostream>
#include <fstream>

#include "IBDmix/vcf_file.h"

TEST(VcfFile, CanHandleArchaic) {
    std::ofstream vcf_file;
    vcf_file.open("test_vcf.txt");
    vcf_file << "##contig=<ID=Y,length=59373566>\n"
        << "##reference=file:///mnt/solexa/Genomes/hg19_1000g/whole_genome.fa\n"
        << "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	AltaiNea\n"
        << "1	10001	.	T	.	51.05	.	.	GT:DP:GQ:PL:A:C:G:T:IR	0/0:249:21.05:0,21,265:1,2:0,1:1,0:116,129:2\n"
        << "1	10002	.	A	.	36.01	.	.	GT:DP:GQ:PL:A:C:G:T:IR	0/0:250:6.01:0,6,69:171,184:0,9:0,0:0,2:0\n"
        << "1	10003	.	AA	.	36.01	.	.	GT:DP:GQ:PL:A:C:G:T:IR	0/0:250:6.02:0,6,72:222,201:1,0:0,0:0,8:0\n"
        << "1	10004	.	C	A.	39.01	.	.	GT:DP:GQ:PL:A:C:G:T:IR	0/0:250:9.02:0,9,97:0,1:268,211:0,0:0,1:0\n"
        << "1	10005	.	C	.	41.99	.	.	GT:DP:GQ:PL:A:C:G:T:IR	./.:250:12:0,12,119:0,0:341,287:0,0:0,0:0\n"
        << "1	10006	.	C	.	36.01	.	.	GT:DP:GQ:PL:A:C:G:T:IR	1/0:249:6.02:0,6,72:0,1:408,347:0,0:2,0:0\n"
        << "1	10007	.	T	.	39.01	.	.	GT:DP:GQ:PL:A:C:G:T:IR	1/1:250:9.02:0,9,97:0,0:0,0:0,0:506,451:0\n";
    vcf_file.close();
    FILE * file = fopen("test_vcf.txt", "rt");
    std::ostringstream output;
    VCF_File vcf(file, output);
    ASSERT_STREQ(output.str().c_str(), "\tAltaiNea");
    ASSERT_EQ(vcf.number_individuals, 1);

    ASSERT_TRUE(vcf.update());
    ASSERT_EQ(vcf.chromosome, 1);
    ASSERT_EQ(vcf.position, 10001);
    ASSERT_EQ(vcf.reference, 'T');
    ASSERT_EQ(vcf.alternative, '.');
    ASSERT_EQ(vcf.genotypes[0], '0');

    ASSERT_TRUE(vcf.update());
    ASSERT_EQ(vcf.chromosome, 1);
    ASSERT_EQ(vcf.position, 10002);
    ASSERT_EQ(vcf.reference, 'A');
    ASSERT_EQ(vcf.alternative, '.');
    ASSERT_EQ(vcf.genotypes[0], '0');

    // skip multiple ref/alt and no genotypes
    ASSERT_TRUE(vcf.update(true));
    ASSERT_EQ(vcf.chromosome, 1);
    ASSERT_EQ(vcf.position, 10006);
    ASSERT_EQ(vcf.reference, 'C');
    ASSERT_EQ(vcf.alternative, '.');
    ASSERT_EQ(vcf.genotypes[0], '1');

    ASSERT_TRUE(vcf.update(true));
    ASSERT_EQ(vcf.chromosome, 1);
    ASSERT_EQ(vcf.position, 10007);
    ASSERT_EQ(vcf.reference, 'T');
    ASSERT_EQ(vcf.alternative, '.');
    ASSERT_EQ(vcf.genotypes[0], '2');

    ASSERT_TRUE(!vcf.update(true));
}

TEST(VcfFile, CanHandleModern) {
    std::ofstream vcf_file;
    vcf_file.open("test_vcf.txt");
    vcf_file << "##log_file=shapeit_10122015_22h01m13s_3f764d75-2fbb-42df-ab75-8c2dfd5731ce.log\n"
        << "##FORMAT=<ID=GT,Number=1,Type=String,Description='Phased Genotype'>\n"
        << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tI1\tI2\tI3\tI4\tI5\n"
        << "1\t846687\t1:846687_C_T\tC\tT\t.\tPASS\t.\tGT\t0|0\t0|0\t0|0\t0|0\t0|0\n"
        << "1\t846688\t1:846688_G_A\tG\tA\t.\tPASS\t.\tGT\t1|0\t1|1\t0|1\t0|0\t.|.\n"
        << "1\t846742\t1:846742_C_T\tC\tTT\t.\tPASS\t.\tGT\t0|0\t0|0\t0|0\t0|0\t0|0\n"
        << "1\t846758\t1:846758_G_A\tG\t.A\t.\tPASS\t.\tGT\t0|0\t0|0\t0|0\t0|0\t0|0\n"
        << "1\t846808\t1:846808_C_T\tC\tT\t.\tPASS\t.\tGT\t0|0\t0|0\t0|1\t0|0\t0|0\n";
    vcf_file.close();
    FILE * file = fopen("test_vcf.txt", "rt");
    std::ostringstream output;
    VCF_File vcf(file, output);
    ASSERT_STREQ(output.str().c_str(), "\tI1\tI2\tI3\tI4\tI5");
    ASSERT_EQ(vcf.number_individuals, 5);

    ASSERT_TRUE(vcf.update());
    ASSERT_EQ(vcf.chromosome, 1);
    ASSERT_EQ(vcf.position, 846687);
    ASSERT_EQ(vcf.reference, 'C');
    ASSERT_EQ(vcf.alternative, 'T');
    ASSERT_EQ(vcf.genotypes[0], '0');
    ASSERT_EQ(vcf.genotypes[2], '0');
    ASSERT_EQ(vcf.genotypes[4], '0');
    ASSERT_EQ(vcf.genotypes[6], '0');
    ASSERT_EQ(vcf.genotypes[8], '0');

    ASSERT_TRUE(vcf.update());
    ASSERT_EQ(vcf.chromosome, 1);
    ASSERT_EQ(vcf.position, 846688);
    ASSERT_EQ(vcf.reference, 'G');
    ASSERT_EQ(vcf.alternative, 'A');
    ASSERT_EQ(vcf.genotypes[0], '1');
    ASSERT_EQ(vcf.genotypes[2], '2');
    ASSERT_EQ(vcf.genotypes[4], '1');
    ASSERT_EQ(vcf.genotypes[6], '0');
    ASSERT_EQ(vcf.genotypes[8], '9');

    // indels
    ASSERT_TRUE(vcf.update());
    ASSERT_EQ(vcf.chromosome, 1);
    ASSERT_EQ(vcf.position, 846742);
    ASSERT_FALSE(vcf.isvalid);

    ASSERT_TRUE(vcf.update());
    ASSERT_EQ(vcf.chromosome, 1);
    ASSERT_EQ(vcf.position, 846758);
    ASSERT_FALSE(vcf.isvalid);

    ASSERT_TRUE(vcf.update());
    ASSERT_EQ(vcf.chromosome, 1);
    ASSERT_EQ(vcf.position, 846808);
    ASSERT_EQ(vcf.reference, 'C');
    ASSERT_EQ(vcf.alternative, 'T');
    ASSERT_EQ(vcf.genotypes[0], '0');
    ASSERT_EQ(vcf.genotypes[2], '0');
    ASSERT_EQ(vcf.genotypes[4], '1');
    ASSERT_EQ(vcf.genotypes[6], '0');
    ASSERT_EQ(vcf.genotypes[8], '0');

    ASSERT_TRUE(!vcf.update());
}
