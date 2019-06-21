#define BOOST_TEST_MODULE mergeVCF
#include <boost/test/included/unit_test.hpp>
#include "../IBDmix/mergeVCF.cpp"
#include <iostream>
#include <fstream>

BOOST_AUTO_TEST_CASE(test_vcf_arch)
{
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
    VCF_File vcf(file);
    BOOST_REQUIRE_EQUAL(vcf.number_individuals, 1);

    BOOST_REQUIRE(vcf.update());
    BOOST_REQUIRE_EQUAL(vcf.chromosome, 1);
    BOOST_REQUIRE_EQUAL(vcf.position, 10001);
    BOOST_REQUIRE_EQUAL(vcf.reference, 'T');
    BOOST_REQUIRE_EQUAL(vcf.alternative, '.');
    BOOST_REQUIRE_EQUAL(vcf.genotypes[0], 0);

    BOOST_REQUIRE(vcf.update());
    BOOST_REQUIRE_EQUAL(vcf.chromosome, 1);
    BOOST_REQUIRE_EQUAL(vcf.position, 10002);
    BOOST_REQUIRE_EQUAL(vcf.reference, 'A');
    BOOST_REQUIRE_EQUAL(vcf.alternative, '.');
    BOOST_REQUIRE_EQUAL(vcf.genotypes[0], 0);

    // skip multiple ref/alt and no genotypes
    BOOST_REQUIRE(vcf.update(true));
    BOOST_REQUIRE_EQUAL(vcf.chromosome, 1);
    BOOST_REQUIRE_EQUAL(vcf.position, 10006);
    BOOST_REQUIRE_EQUAL(vcf.reference, 'C');
    BOOST_REQUIRE_EQUAL(vcf.alternative, '.');
    BOOST_REQUIRE_EQUAL(vcf.genotypes[0], 1);

    BOOST_REQUIRE(vcf.update(true));
    BOOST_REQUIRE_EQUAL(vcf.chromosome, 1);
    BOOST_REQUIRE_EQUAL(vcf.position, 10007);
    BOOST_REQUIRE_EQUAL(vcf.reference, 'T');
    BOOST_REQUIRE_EQUAL(vcf.alternative, '.');
    BOOST_REQUIRE_EQUAL(vcf.genotypes[0], 2);

    BOOST_REQUIRE(!vcf.update(true));
}

BOOST_AUTO_TEST_CASE(test_vcf_mod)
{
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
    VCF_File vcf(file);
    BOOST_REQUIRE_EQUAL(vcf.number_individuals, 5);

    BOOST_REQUIRE(vcf.update());
    BOOST_REQUIRE_EQUAL(vcf.chromosome, 1);
    BOOST_REQUIRE_EQUAL(vcf.position, 846687);
    BOOST_REQUIRE_EQUAL(vcf.reference, 'C');
    BOOST_REQUIRE_EQUAL(vcf.alternative, 'T');
    BOOST_REQUIRE_EQUAL(vcf.genotypes[0], 0);
    BOOST_REQUIRE_EQUAL(vcf.genotypes[1], 0);
    BOOST_REQUIRE_EQUAL(vcf.genotypes[2], 0);
    BOOST_REQUIRE_EQUAL(vcf.genotypes[3], 0);
    BOOST_REQUIRE_EQUAL(vcf.genotypes[4], 0);

    BOOST_REQUIRE(vcf.update());
    BOOST_REQUIRE_EQUAL(vcf.chromosome, 1);
    BOOST_REQUIRE_EQUAL(vcf.position, 846688);
    BOOST_REQUIRE_EQUAL(vcf.reference, 'G');
    BOOST_REQUIRE_EQUAL(vcf.alternative, 'A');
    BOOST_REQUIRE_EQUAL(vcf.genotypes[0], 1);
    BOOST_REQUIRE_EQUAL(vcf.genotypes[1], 2);
    BOOST_REQUIRE_EQUAL(vcf.genotypes[2], 1);
    BOOST_REQUIRE_EQUAL(vcf.genotypes[3], 0);
    BOOST_REQUIRE_EQUAL(vcf.genotypes[4], 9);

    BOOST_REQUIRE(vcf.update());
    BOOST_REQUIRE_EQUAL(vcf.chromosome, 1);
    BOOST_REQUIRE_EQUAL(vcf.position, 846808);
    BOOST_REQUIRE_EQUAL(vcf.reference, 'C');
    BOOST_REQUIRE_EQUAL(vcf.alternative, 'T');
    BOOST_REQUIRE_EQUAL(vcf.genotypes[0], 0);
    BOOST_REQUIRE_EQUAL(vcf.genotypes[1], 0);
    BOOST_REQUIRE_EQUAL(vcf.genotypes[2], 1);
    BOOST_REQUIRE_EQUAL(vcf.genotypes[3], 0);
    BOOST_REQUIRE_EQUAL(vcf.genotypes[4], 0);

    BOOST_REQUIRE(!vcf.update());
}
