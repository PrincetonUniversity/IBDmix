#define BOOST_TEST_MODULE Genotype_Reader
#include <boost/test/included/unit_test.hpp>
#include "../IBDmix/Genotype_Reader.hpp"
#include <iostream>
#include <fstream>

BOOST_AUTO_TEST_CASE(setup){
    std::ofstream genotype;
    genotype.open("test_genotype.txt");
    genotype << "chrom\tpos\tref\talt\tn1\tm1\tm2\tm3\tm4\n"
        << "1\t2\tA\tT\t1\t0\t0\t0\t0\n"
        << "1\t3\tA\tT\t2\t0\t0\t0\t0\n"
        << "1\t4\tA\tT\t1\t0\t1\t1\t1\n"
        << "1\t104\tA\tT\t1\t0\t1\t1\t1\n"
        << "1\t105\tA\tT\t0\t2\t1\t1\t1\n"
        << "2\t125\tA\tT\t0\t2\t2\t1\t1\n"
        << "3\t126\tA\tT\t0\t2\t2\t2\t2\n"
        ;
    genotype.close();
}

BOOST_AUTO_TEST_CASE(test_initializer)
{
    FILE * file = fopen("test_genotype.txt", "rt");
    // use explicit values in case defaults change
    Genotype_Reader reader = Genotype_Reader(
            file, nullptr, 0.01, 0.002, 2, 1e-200);
    BOOST_REQUIRE_EQUAL(4, reader.initialize_arrays(nullptr));
    rewind(file);
    std::ofstream samples;
    samples.open("test_samples.txt");
    samples << "m2\nm4\n";
    samples.close();
    BOOST_REQUIRE_EQUAL(2, reader.initialize_arrays(fopen("test_samples.txt", "rt")));
    rewind(file);
    samples.open("test_samples.txt");
    samples << "m2\nm4";
    samples.close();
    BOOST_REQUIRE_EQUAL(2, reader.initialize_arrays(fopen("test_samples.txt", "rt")));
    rewind(file);
    samples.open("test_samples.txt");
    samples << "m2";
    samples.close();
    BOOST_REQUIRE_EQUAL(1, reader.initialize_arrays(fopen("test_samples.txt", "rt")));
}

BOOST_AUTO_TEST_CASE(test_find_token)
{
    BOOST_REQUIRE_EQUAL(2, find_token("tes", "test\tte\ttes"));
    BOOST_REQUIRE_EQUAL(-1, find_token("not in", "test\tstring\t"));
    BOOST_REQUIRE_EQUAL(0, find_token("test", "test\tstring\t"));
    BOOST_REQUIRE_EQUAL(1, find_token("string", "test\tstring\t"));
    BOOST_REQUIRE_EQUAL(1, find_token("string", "test\tstring"));
    BOOST_REQUIRE_EQUAL(1, find_token("string\tthing", "test\tstring"));
    BOOST_REQUIRE_EQUAL(1, find_token("string\tthing", "test\tstring\tthings"));
    BOOST_REQUIRE_EQUAL(-1, find_token("strinG", "test\tstring"));
    BOOST_REQUIRE_EQUAL(-1, find_token("string", "test\tstrinG"));
    BOOST_REQUIRE_EQUAL(-1, find_token("string1", "test\tstring"));
    BOOST_REQUIRE_EQUAL(-1, find_token("string", "test\tstring1"));
    BOOST_REQUIRE_EQUAL(-1, find_token("String", "test\tstring1"));
    BOOST_REQUIRE_EQUAL(-1, find_token("\0", "test\tstring1"));
}

BOOST_AUTO_TEST_CASE(test_find_archaic)
{
    Genotype_Reader reader = Genotype_Reader(nullptr);
    reader.samples = new char[100];
    reader.buffer = new char[100];
    strcpy(reader.buffer, "n1\ts1\ts2");

    // test both null, delete first token in samples
    strcpy(reader.samples, "n1\ts1\ts2");
    reader.find_archaic(nullptr, nullptr);
    BOOST_REQUIRE_EQUAL(0, reader.archaic_index);
    BOOST_REQUIRE_EQUAL(reader.samples,  "s1\ts2");

    // test with sample line, delete nothing
    strcpy(reader.samples, "n1\ts1\ts2");
    reader.find_archaic(nullptr, "");
    BOOST_REQUIRE_EQUAL(0, reader.archaic_index);
    BOOST_REQUIRE_EQUAL(reader.samples,  "n1\ts1\ts2");

    // test with archaic
    strcpy(reader.samples, "n1\ts1\ts2");
    reader.find_archaic("s1", nullptr);
    BOOST_REQUIRE_EQUAL(1, reader.archaic_index);
    BOOST_REQUIRE_EQUAL(reader.samples,  "n1\ts2");

    strcpy(reader.samples, "n1\ts1\ts2");
    reader.find_archaic("s2", nullptr);
    BOOST_REQUIRE_EQUAL(2, reader.archaic_index);
    BOOST_REQUIRE_EQUAL(reader.samples,  "n1\ts1\t");

    strcpy(reader.samples, "n1\ts1\ts2\t");
    reader.find_archaic("s2", nullptr);
    BOOST_REQUIRE_EQUAL(2, reader.archaic_index);
    BOOST_REQUIRE_EQUAL(reader.samples,  "n1\ts1\t");

    strcpy(reader.samples, "n1\ts1\ts2\t");
    reader.find_archaic("s2", "");
    BOOST_REQUIRE_EQUAL(2, reader.archaic_index);
    BOOST_REQUIRE_EQUAL(reader.samples,  "n1\ts1\ts2\t");
}

BOOST_AUTO_TEST_CASE(test_determine_sample_mapping)
{
    Genotype_Reader reader = Genotype_Reader(nullptr);
    reader.buffer = new char[100];
    strcpy(reader.buffer, "n1\ts1\ts2");
    reader.sample_to_index = new int[10];
    reader.num_samples = 10;

    reader.archaic_index = 0;
    reader.determine_sample_mapping(NULL);
    for(int i = 0; i < reader.num_samples; i++)
        BOOST_REQUIRE_EQUAL(i+1, reader.sample_to_index[i]);

    reader.archaic_index = 10;
    reader.determine_sample_mapping(NULL);
    for(int i = 0; i < reader.num_samples; i++)
        BOOST_REQUIRE_EQUAL(i, reader.sample_to_index[i]);

    reader.archaic_index = 5;
    reader.determine_sample_mapping(NULL);
    for(int i = 0; i < 5; i++)
        BOOST_REQUIRE_EQUAL(i, reader.sample_to_index[i]);
    for(int i = 5; i < reader.num_samples; i++)
        BOOST_REQUIRE_EQUAL(i+1, reader.sample_to_index[i]);

    strcpy(reader.buffer, "n1\ts1\ts2\ts3\ts4\ts5\ts6\ts7");
    reader.num_samples = 3;
    reader.determine_sample_mapping("s6\tn1\ts3");
    int result[] = {6, 0, 3};
    for(int i = 0; i < reader.num_samples; i++)
        BOOST_REQUIRE_EQUAL(result[i], reader.sample_to_index[i]);

    reader.determine_sample_mapping("s1\ts1\ts1");
    int result2[] = {1, 1, 1};
    for(int i = 0; i < reader.num_samples; i++)
        BOOST_REQUIRE_EQUAL(result2[i], reader.sample_to_index[i]);
}

BOOST_AUTO_TEST_CASE(test_yield_sample)
{
    Genotype_Reader reader = Genotype_Reader(nullptr);
    reader.num_samples = 3;
    reader.samples = new char[100];
    strcpy(reader.samples, "n1\ts1\ts2");
    char *ptr = nullptr;
    int count = 0;

    BOOST_REQUIRE(reader.yield_sample(ptr, count++));
    BOOST_REQUIRE_EQUAL(ptr, "n1");

    BOOST_REQUIRE(reader.yield_sample(ptr, count++));
    BOOST_REQUIRE_EQUAL(ptr, "s1");

    BOOST_REQUIRE(reader.yield_sample(ptr, count++));
    BOOST_REQUIRE_EQUAL(ptr, "s2");

    BOOST_REQUIRE(!reader.yield_sample(ptr, count++));
    BOOST_REQUIRE_EQUAL(reader.samples, "n1\ts1\ts2");
}

double original_cal_lod(int source_gt, int target_gt, double pb, double aerr, double merr,
        double minesp){
    double pa = 1-pb;
    if(source_gt == 0 && target_gt == 0)
        return log10(((1-aerr)*(1-merr)+aerr*merr) / pa / (1-aerr*(1-aerr)));
    if(source_gt == 0 && target_gt == 1){
        double temp = 0.5 * ((1-aerr)*merr+aerr*(1-merr)) / pb /
                (1-aerr*(1-aerr)) + 0.5 * ((1-aerr)*(1-merr)+aerr*merr) /
                pa / (1-aerr*(1-aerr));
        return log10(temp);
    }
    if(source_gt == 0 && target_gt == 2)
        return log10(((1-aerr)*merr+aerr*(1-merr) + minesp) /
                pb / (1-aerr*(1-aerr) + minesp));
    if(source_gt == 1 && target_gt == 0)
        return -log10(pa * (1+2*aerr*(1-aerr)));
    if(source_gt == 1 && target_gt == 1)
        return -log10(2*pa*pb*(1+2*aerr*(1-aerr)));
    if(source_gt == 1 && target_gt == 2)
        return -log10(pb * (1+2*aerr*(1-aerr)));
    if(source_gt == 2 && target_gt == 0)
        return log10(((1-aerr)*merr+aerr*(1-merr) + minesp) /
                pa / (1-aerr*(1-aerr) + minesp));
    if(source_gt == 2 && target_gt == 1)
        return log10(0.5*((1-aerr)*(1-merr)+aerr*merr) / pb /
                (1-aerr*(1-aerr)) + 0.5 * ((1-aerr)*merr+aerr*(1-merr)) /
                pa / (1-aerr*(1-aerr)));
    if(source_gt == 2 && target_gt == 2)
        return log10(((1-aerr)*(1-merr)+aerr*merr) / pb / (1-aerr*(1-aerr)));
    return 0;
}

BOOST_AUTO_TEST_CASE(test_calculate_lod)
{
    Genotype_Reader reader = Genotype_Reader(nullptr);

    double pbs[] = {0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.99, 1};
    double aerrs[] = {0.01, 0.02, 0.03};
    double merrs[] = {0.0001, 0.001, 0.002};
    double minesps[] = {1e-200, 1e-190, 1e-100};
    char gts[] = {'0', '1', '2', '9'};
    for(double pb : pbs) for(double aerr : aerrs) for(double merr : merrs)
        for(double minesp : minesps) for(char arch : gts){
            reader.archaic_error = aerr;
            reader.minesp = minesp;
            reader.update_lod_cache(arch, pb, merr);
            for(char mod : gts){
                double read = reader.calculate_lod(mod);
                double orig = original_cal_lod(arch-'0', mod-'0', pb, aerr, merr, minesp);
                if(isinf(read) && isinf(orig))
                    continue;
                BOOST_REQUIRE_CLOSE(read, orig, 0.00001);
            }
            reader.update_lod_cache(arch, pb, merr, false); // not selected
            for(char mod : gts){
                double read = reader.calculate_lod(mod);
                double orig = original_cal_lod(arch-'0', mod-'0', pb, aerr, merr, minesp);
                if((mod == '0' && arch == '2') || (mod == '2' && arch == '0')){
                    if(isinf(read) && isinf(orig))
                        continue;
                    //std::cout << "mod=" << mod-'0' << " arch=" << arch-'0' << 
                    //        " pb=" << pb << " aerr=" << aerr << " merr=" << merr
                    //        << " minesp=" << minesp << "\n"
                    //        << "reader=" << read
                    //        << "\norig=  " << orig
                    //        << "\n";
                    BOOST_REQUIRE_CLOSE(read, orig, 0.00001);
                }
                else{
                    BOOST_REQUIRE_EQUAL(0, read);
                }
            }
        }

}

BOOST_AUTO_TEST_CASE(test_get_frequency)
{
    double result = 0;
    Genotype_Reader reader = Genotype_Reader(nullptr);
    reader.sample_to_index = new int[4];
    for(int i = 0; i < 4; i++)
        reader.sample_to_index[i] = i+1;
    reader.num_samples = 4;
    reader.buffer = new char[100];
    reader.minor_allele_cutoff = 1;

    strcpy(reader.buffer, "0 0 1 0 0");
    // fail since cutoff too low
    BOOST_REQUIRE(!reader.get_frequency(result));
    reader.minor_allele_cutoff = 0;
    BOOST_REQUIRE(reader.get_frequency(result));
    BOOST_REQUIRE_EQUAL(0.125, result);

    reader.minor_allele_cutoff = 1;
    strcpy(reader.buffer, "0 0 2 0 0");
    BOOST_REQUIRE(reader.get_frequency(result));
    BOOST_REQUIRE_EQUAL(0.25, result);

    reader.minor_allele_cutoff = 0;
    strcpy(reader.buffer, "0 0 2 9 9");
    BOOST_REQUIRE(reader.get_frequency(result));
    BOOST_REQUIRE_EQUAL(0.5, result);

    strcpy(reader.buffer, "0 9 2 9 9");
    BOOST_REQUIRE(!reader.get_frequency(result));
    BOOST_REQUIRE_EQUAL(1, result);

    strcpy(reader.buffer, "0 9 9 9 9");
    BOOST_REQUIRE(!reader.get_frequency(result));
    BOOST_REQUIRE_EQUAL(0, result);

    strcpy(reader.buffer, "0 0 0 0 0");
    BOOST_REQUIRE(!reader.get_frequency(result));
    BOOST_REQUIRE_EQUAL(0, result);

}

BOOST_AUTO_TEST_CASE(test_get_modern_error)
{
    Genotype_Reader reader = Genotype_Reader(nullptr);
    reader.modern_error_max = 0.1;
    reader.modern_error_proportion = 2;

    BOOST_REQUIRE_CLOSE(0.1, reader.get_modern_error(0.4), 0.00001);
    BOOST_REQUIRE_CLOSE(0.1, reader.get_modern_error(0.6), 0.00001);
    BOOST_REQUIRE_CLOSE(0.1, reader.get_modern_error(0.1), 0.00001);
    BOOST_REQUIRE_CLOSE(0.1, reader.get_modern_error(0.9), 0.00001);
    BOOST_REQUIRE_CLOSE(0.02, reader.get_modern_error(0.01), 0.00001);
    BOOST_REQUIRE_CLOSE(0.02, reader.get_modern_error(0.99), 0.00001);
}

BOOST_AUTO_TEST_CASE(test_process_line_buffer){
    Genotype_Reader reader = Genotype_Reader(nullptr);
    reader.buffer = new char[100];
    reader.num_samples = 4;
    reader.sample_to_index = new int[4];
    for(int i = 0; i < 4; i++)
        reader.sample_to_index[i] = 4-i;
    reader.archaic_index = 0;
    reader.lod_scores = new double[4];

    strcpy(reader.buffer, "0 0 1 2 9");
    reader.process_line_buffer(true);
    double result[] = {0, -1.617, 0.00432, 0.300};
    for(int i = 0; i < 4; i++)
        BOOST_REQUIRE_CLOSE(result[i], reader.lod_scores[i], 0.1);

    reader.process_line_buffer(false);
    result[0] = 0;
    result[1] = -1.617;
    result[2] = 0;
    result[3] = 0;
    for(int i = 0; i < 4; i++)
        BOOST_REQUIRE_CLOSE(result[i], reader.lod_scores[i], 0.1);

    strcpy(reader.buffer, "9 0 1 2 9");
    reader.process_line_buffer(true);
    result[1] = 0;
    for(int i = 0; i < 4; i++)
        BOOST_REQUIRE_CLOSE(result[i], reader.lod_scores[i], 0.1);

    reader.process_line_buffer(false);
    for(int i = 0; i < 4; i++)
        BOOST_REQUIRE_CLOSE(result[i], reader.lod_scores[i], 0.1);

    strcpy(reader.buffer, "0 0 1 0 9");
    reader.process_line_buffer(true);
    for(int i = 0; i < 4; i++)
        BOOST_REQUIRE_CLOSE(result[i], reader.lod_scores[i], 0.1);

    strcpy(reader.buffer, "1 1 1 2 9");
    reader.process_line_buffer(true);
    result[1] = 0.16758;
    result[2] = 0.34366;
    result[3] = 0.34366;
    for(int i = 0; i < 4; i++)
        BOOST_REQUIRE_CLOSE(result[i], reader.lod_scores[i], 0.1);
        //std::cout << reader.lod_scores[i] << "\n";
}

BOOST_AUTO_TEST_CASE(test_in_mask){
    std::ofstream mask;
    mask.open("test_mask.txt");
    mask << "1 100 120\n"
        << "1 130 140\n"
        << "1 160 161\n"
        << "1 190 200\n"
        << "1 260 281\n"
        << "2 130 140\n"
        << "4 130 140\n"
        << "8 130 140\n"
        ;
    mask.close();
    Genotype_Reader reader = Genotype_Reader(
            nullptr, fopen("test_mask.txt", "rt"));
    BOOST_REQUIRE_EQUAL(reader.mask_chromosome, -1);

    reader.chromosome = 1;
    reader.position = 90;
    BOOST_REQUIRE(!reader.in_mask());
    BOOST_REQUIRE_EQUAL(reader.mask_chromosome, 1);
    BOOST_REQUIRE_EQUAL(reader.mask_start, 100);
    BOOST_REQUIRE_EQUAL(reader.mask_end, 120);

    // check same position
    BOOST_REQUIRE(!reader.in_mask());
    BOOST_REQUIRE(!reader.in_mask());
    BOOST_REQUIRE(!reader.in_mask());

    // should not update mask
    BOOST_REQUIRE_EQUAL(reader.mask_chromosome, 1);
    BOOST_REQUIRE_EQUAL(reader.mask_start, 100);
    BOOST_REQUIRE_EQUAL(reader.mask_end, 120);

    reader.position = 100;
    BOOST_REQUIRE(!reader.in_mask());

    reader.position = 101;
    BOOST_REQUIRE(reader.in_mask());

    reader.position = 102;
    BOOST_REQUIRE(reader.in_mask());

    reader.position = 120;
    BOOST_REQUIRE(reader.in_mask());
    BOOST_REQUIRE_EQUAL(reader.mask_chromosome, 1);
    BOOST_REQUIRE_EQUAL(reader.mask_start, 100);
    BOOST_REQUIRE_EQUAL(reader.mask_end, 120);

    // skip a range
    reader.position = 161;
    BOOST_REQUIRE(reader.in_mask());
    BOOST_REQUIRE_EQUAL(reader.mask_chromosome, 1);
    BOOST_REQUIRE_EQUAL(reader.mask_start, 160);
    BOOST_REQUIRE_EQUAL(reader.mask_end, 161);

    reader.position = 261;
    BOOST_REQUIRE(reader.in_mask());
    BOOST_REQUIRE_EQUAL(reader.mask_chromosome, 1);
    BOOST_REQUIRE_EQUAL(reader.mask_start, 260);
    BOOST_REQUIRE_EQUAL(reader.mask_end, 281);
    
    //skip chromosomes
    reader.chromosome = 2;
    reader.position = 130;
    BOOST_REQUIRE(!reader.in_mask());
    BOOST_REQUIRE_EQUAL(reader.mask_chromosome, 2);
    BOOST_REQUIRE_EQUAL(reader.mask_start, 130);
    BOOST_REQUIRE_EQUAL(reader.mask_end, 140);

    reader.chromosome = 3;
    reader.position = 130;
    BOOST_REQUIRE(!reader.in_mask());
    BOOST_REQUIRE_EQUAL(reader.mask_chromosome, 4);
    BOOST_REQUIRE_EQUAL(reader.mask_start, 130);
    BOOST_REQUIRE_EQUAL(reader.mask_end, 140);

    reader.chromosome = 4;
    reader.position = 131;
    BOOST_REQUIRE(reader.in_mask());
    BOOST_REQUIRE_EQUAL(reader.mask_chromosome, 4);
    BOOST_REQUIRE_EQUAL(reader.mask_start, 130);
    BOOST_REQUIRE_EQUAL(reader.mask_end, 140);

    reader.chromosome = 5;
    reader.position = 131;
    BOOST_REQUIRE(!reader.in_mask());
    BOOST_REQUIRE_EQUAL(reader.mask_chromosome, 8);
    BOOST_REQUIRE_EQUAL(reader.mask_start, 130);
    BOOST_REQUIRE_EQUAL(reader.mask_end, 140);

    reader.chromosome = 6;
    reader.position = 131;
    BOOST_REQUIRE(!reader.in_mask());
    BOOST_REQUIRE_EQUAL(reader.mask_chromosome, 8);
    BOOST_REQUIRE_EQUAL(reader.mask_start, 130);
    BOOST_REQUIRE_EQUAL(reader.mask_end, 140);

    reader.chromosome = 8;
    reader.position = 131;
    BOOST_REQUIRE(reader.in_mask());
    BOOST_REQUIRE_EQUAL(reader.mask_chromosome, 8);
    BOOST_REQUIRE_EQUAL(reader.mask_start, 130);
    BOOST_REQUIRE_EQUAL(reader.mask_end, 140);

    reader.chromosome = 9;
    reader.position = 131;
    BOOST_REQUIRE(!reader.in_mask());
    BOOST_REQUIRE(!reader.in_mask());
    BOOST_REQUIRE(!reader.in_mask());
    BOOST_REQUIRE(!reader.in_mask());
    BOOST_REQUIRE_EQUAL(reader.mask_chromosome, 8);
    BOOST_REQUIRE_EQUAL(reader.mask_start, 130);
    BOOST_REQUIRE_EQUAL(reader.mask_end, 140);

    reader.chromosome = 8;
    // still retain the last read values (not actual use case)
    BOOST_REQUIRE(reader.in_mask());
    fclose(reader.mask);
    // setting mask to null should short the region check
    reader.mask = nullptr;
    BOOST_REQUIRE(!reader.in_mask());
}

BOOST_AUTO_TEST_CASE(test_update_defaults)
{
    Genotype_Reader reader = Genotype_Reader(fopen("test_genotype.txt", "rt"));
    reader.initialize_arrays();
    BOOST_REQUIRE_EQUAL(4, reader.num_samples);

    //"1\t2\tA\tT\t1\t0\t0\t0\t0\n"  fails allele check
    BOOST_REQUIRE(reader.update());
    BOOST_REQUIRE_EQUAL(1, reader.chromosome);
    BOOST_REQUIRE_EQUAL(2, reader.position);
    BOOST_REQUIRE_CLOSE(0, reader.lod_scores[0], 0.0001);
    BOOST_REQUIRE_CLOSE(0, reader.lod_scores[1], 0.0001);
    BOOST_REQUIRE_CLOSE(0, reader.lod_scores[2], 0.0001);
    BOOST_REQUIRE_CLOSE(0, reader.lod_scores[3], 0.0001);

    // "1\t3\tA\tT\t2\t0\t0\t0\t0\n" fails check, 2/0 override
    BOOST_REQUIRE(reader.update());
    BOOST_REQUIRE_EQUAL(1, reader.chromosome);
    BOOST_REQUIRE_EQUAL(3, reader.position);
    BOOST_REQUIRE_CLOSE(-1.99568, reader.lod_scores[0], 0.0001);
    BOOST_REQUIRE_CLOSE(-1.99568, reader.lod_scores[1], 0.0001);
    BOOST_REQUIRE_CLOSE(-1.99568, reader.lod_scores[2], 0.0001);
    BOOST_REQUIRE_CLOSE(-1.99568, reader.lod_scores[3], 0.0001);

    // "1\t4\tA\tT\t1\t0\t1\t1\t1\n"  passes check
    BOOST_REQUIRE(reader.update());
    BOOST_REQUIRE_EQUAL(1, reader.chromosome);
    BOOST_REQUIRE_EQUAL(4, reader.position);
    BOOST_REQUIRE_CLOSE(0.195605, reader.lod_scores[0], 0.0001);
    BOOST_REQUIRE_CLOSE(0.320544, reader.lod_scores[1], 0.0001);
    BOOST_REQUIRE_CLOSE(0.320544, reader.lod_scores[2], 0.0001);
    BOOST_REQUIRE_CLOSE(0.320544, reader.lod_scores[3], 0.0001);

    // "1\t104\tA\tT\t1\t0\t1\t1\t1\n" same as above
    BOOST_REQUIRE(reader.update());
    BOOST_REQUIRE_EQUAL(1, reader.chromosome);
    BOOST_REQUIRE_EQUAL(104, reader.position);
    BOOST_REQUIRE_CLOSE(0.195605, reader.lod_scores[0], 0.0001);
    BOOST_REQUIRE_CLOSE(0.320544, reader.lod_scores[1], 0.0001);
    BOOST_REQUIRE_CLOSE(0.320544, reader.lod_scores[2], 0.0001);
    BOOST_REQUIRE_CLOSE(0.320544, reader.lod_scores[3], 0.0001);

    // "1\t105\tA\tT\t0\t2\t1\t1\t1\n"  passes check
    BOOST_REQUIRE(reader.update());
    BOOST_REQUIRE_EQUAL(1, reader.chromosome);
    BOOST_REQUIRE_EQUAL(105, reader.position);
    BOOST_REQUIRE_CLOSE(-1.713827, reader.lod_scores[0], 0.0001);
    BOOST_REQUIRE_CLOSE(0.127177, reader.lod_scores[1], 0.0001);
    BOOST_REQUIRE_CLOSE(0.127177, reader.lod_scores[2], 0.0001);
    BOOST_REQUIRE_CLOSE(0.127177, reader.lod_scores[3], 0.0001);

    // "2\t125\tA\tT\t0\t2\t2\t1\t1\n" change chromosome
    BOOST_REQUIRE(reader.update());
    BOOST_REQUIRE_EQUAL(2, reader.chromosome);
    BOOST_REQUIRE_EQUAL(125, reader.position);
    BOOST_REQUIRE_CLOSE(-1.79301, reader.lod_scores[0], 0.0001);
    BOOST_REQUIRE_CLOSE(-1.79301, reader.lod_scores[1], 0.0001);
    BOOST_REQUIRE_CLOSE(0.301874, reader.lod_scores[2], 0.0001);
    BOOST_REQUIRE_CLOSE(0.301874, reader.lod_scores[3], 0.0001);

    // "3\t126\tA\tT\t0\t2\t2\t2\t2\n" change chrom, 0/2 override check
    BOOST_REQUIRE(reader.update());
    BOOST_REQUIRE_EQUAL(3, reader.chromosome);
    BOOST_REQUIRE_EQUAL(126, reader.position);
    BOOST_REQUIRE_CLOSE(-1.99568, reader.lod_scores[0], 0.0001);
    BOOST_REQUIRE_CLOSE(-1.99568, reader.lod_scores[1], 0.0001);
    BOOST_REQUIRE_CLOSE(-1.99568, reader.lod_scores[2], 0.0001);
    BOOST_REQUIRE_CLOSE(-1.99568, reader.lod_scores[3], 0.0001);

    // eof
    BOOST_REQUIRE(!reader.update());
    BOOST_REQUIRE(!reader.update());
    BOOST_REQUIRE(!reader.update());
}

BOOST_AUTO_TEST_CASE(test_update_samples)
{
    Genotype_Reader reader = Genotype_Reader(fopen("test_genotype.txt", "rt"));
    std::ofstream samples;
    samples.open("test_samples.txt");
    samples << "m4\nm2\n";
    samples.close();
    reader.initialize_arrays(fopen("test_samples.txt", "rt"), "m1");
    BOOST_REQUIRE_EQUAL(2, reader.num_samples);

    //"1\t2\tA\tT\t1\t0\t0\t0\t0\n"  fails allele check
    BOOST_REQUIRE(reader.update());
    BOOST_REQUIRE_EQUAL(1, reader.chromosome);
    BOOST_REQUIRE_EQUAL(2, reader.position);
    BOOST_REQUIRE_CLOSE(0, reader.lod_scores[0], 0.0001);
    BOOST_REQUIRE_CLOSE(0, reader.lod_scores[1], 0.0001);

    // "1\t3\tA\tT\t2\t0\t0\t0\t0\n" fails check
    BOOST_REQUIRE(reader.update());
    BOOST_REQUIRE_EQUAL(1, reader.chromosome);
    BOOST_REQUIRE_EQUAL(3, reader.position);
    BOOST_REQUIRE_CLOSE(0, reader.lod_scores[0], 0.0001);
    BOOST_REQUIRE_CLOSE(0, reader.lod_scores[1], 0.0001);

    // "1\t4\tA\tT\t1\t0\t1\t1\t1\n"  passes check
    BOOST_REQUIRE(reader.update());
    BOOST_REQUIRE_EQUAL(1, reader.chromosome);
    BOOST_REQUIRE_EQUAL(4, reader.position);
    BOOST_REQUIRE_CLOSE(0.00432094, reader.lod_scores[0], 0.0001);
    BOOST_REQUIRE_CLOSE(0.00432094, reader.lod_scores[1], 0.0001);

    // "1\t104\tA\tT\t1\t0\t1\t1\t1\n" same as above
    BOOST_REQUIRE(reader.update());
    BOOST_REQUIRE_EQUAL(1, reader.chromosome);
    BOOST_REQUIRE_EQUAL(104, reader.position);
    BOOST_REQUIRE_CLOSE(0.00432094, reader.lod_scores[0], 0.0001);
    BOOST_REQUIRE_CLOSE(0.00432094, reader.lod_scores[1], 0.0001);

    // "1\t105\tA\tT\t0\t2\t1\t1\t1\n"  passes check
    BOOST_REQUIRE(reader.update());
    BOOST_REQUIRE_EQUAL(1, reader.chromosome);
    BOOST_REQUIRE_EQUAL(105, reader.position);
    BOOST_REQUIRE_CLOSE(0.00432094, reader.lod_scores[0], 0.0001);
    BOOST_REQUIRE_CLOSE(0.00432094, reader.lod_scores[1], 0.0001);

    // "2\t125\tA\tT\t0\t2\t2\t1\t1\n" change chromosome
    BOOST_REQUIRE(reader.update());
    BOOST_REQUIRE_EQUAL(2, reader.chromosome);
    BOOST_REQUIRE_EQUAL(125, reader.position);
    BOOST_REQUIRE_CLOSE(0, reader.lod_scores[0], 0.0001);
    BOOST_REQUIRE_CLOSE(0, reader.lod_scores[1], 0.0001);

    // "3\t126\tA\tT\t0\t2\t2\t2\t2\n" change chrom, 0/2 override check
    BOOST_REQUIRE(reader.update());
    BOOST_REQUIRE_EQUAL(3, reader.chromosome);
    BOOST_REQUIRE_EQUAL(126, reader.position);
    BOOST_REQUIRE_CLOSE(0, reader.lod_scores[0], 0.0001);
    BOOST_REQUIRE_CLOSE(0, reader.lod_scores[1], 0.0001);

    // eof
    BOOST_REQUIRE(!reader.update());
    BOOST_REQUIRE(!reader.update());
    BOOST_REQUIRE(!reader.update());
}

BOOST_AUTO_TEST_CASE(test_update_mask)
{
    Genotype_Reader reader = Genotype_Reader(
            fopen("test_genotype.txt", "rt"),
            fopen("test_mask.txt", "rt"));
    reader.initialize_arrays();
    BOOST_REQUIRE_EQUAL(4, reader.num_samples);

    //"1\t2\tA\tT\t1\t0\t0\t0\t0\n"  fails allele check
    BOOST_REQUIRE(reader.update());
    BOOST_REQUIRE_EQUAL(1, reader.chromosome);
    BOOST_REQUIRE_EQUAL(2, reader.position);
    BOOST_REQUIRE_CLOSE(0, reader.lod_scores[0], 0.0001);
    BOOST_REQUIRE_CLOSE(0, reader.lod_scores[1], 0.0001);
    BOOST_REQUIRE_CLOSE(0, reader.lod_scores[2], 0.0001);
    BOOST_REQUIRE_CLOSE(0, reader.lod_scores[3], 0.0001);

    // "1\t3\tA\tT\t2\t0\t0\t0\t0\n" fails check, 2/0 override
    BOOST_REQUIRE(reader.update());
    BOOST_REQUIRE_EQUAL(1, reader.chromosome);
    BOOST_REQUIRE_EQUAL(3, reader.position);
    BOOST_REQUIRE_CLOSE(-1.99568, reader.lod_scores[0], 0.0001);
    BOOST_REQUIRE_CLOSE(-1.99568, reader.lod_scores[1], 0.0001);
    BOOST_REQUIRE_CLOSE(-1.99568, reader.lod_scores[2], 0.0001);
    BOOST_REQUIRE_CLOSE(-1.99568, reader.lod_scores[3], 0.0001);

    // "1\t4\tA\tT\t1\t0\t1\t1\t1\n"  passes check
    BOOST_REQUIRE(reader.update());
    BOOST_REQUIRE_EQUAL(1, reader.chromosome);
    BOOST_REQUIRE_EQUAL(4, reader.position);
    BOOST_REQUIRE_CLOSE(0.195605, reader.lod_scores[0], 0.0001);
    BOOST_REQUIRE_CLOSE(0.320544, reader.lod_scores[1], 0.0001);
    BOOST_REQUIRE_CLOSE(0.320544, reader.lod_scores[2], 0.0001);
    BOOST_REQUIRE_CLOSE(0.320544, reader.lod_scores[3], 0.0001);

    // "1\t104\tA\tT\t1\t0\t1\t1\t1\n" in mask
    BOOST_REQUIRE(reader.update());
    BOOST_REQUIRE_EQUAL(1, reader.chromosome);
    BOOST_REQUIRE_EQUAL(104, reader.position);
    BOOST_REQUIRE_CLOSE(0, reader.lod_scores[0], 0.0001);
    BOOST_REQUIRE_CLOSE(0, reader.lod_scores[1], 0.0001);
    BOOST_REQUIRE_CLOSE(0, reader.lod_scores[2], 0.0001);
    BOOST_REQUIRE_CLOSE(0, reader.lod_scores[3], 0.0001);

    // "1\t105\tA\tT\t0\t2\t1\t1\t1\n"  in mask, 0, 2 override
    BOOST_REQUIRE(reader.update());
    BOOST_REQUIRE_EQUAL(1, reader.chromosome);
    BOOST_REQUIRE_EQUAL(105, reader.position);
    BOOST_REQUIRE_CLOSE(-1.713827, reader.lod_scores[0], 0.0001);
    BOOST_REQUIRE_CLOSE(0, reader.lod_scores[1], 0.0001);
    BOOST_REQUIRE_CLOSE(0, reader.lod_scores[2], 0.0001);
    BOOST_REQUIRE_CLOSE(0, reader.lod_scores[3], 0.0001);

    // "2\t125\tA\tT\t0\t2\t2\t1\t1\n" change chromosome
    BOOST_REQUIRE(reader.update());
    BOOST_REQUIRE_EQUAL(2, reader.chromosome);
    BOOST_REQUIRE_EQUAL(125, reader.position);
    BOOST_REQUIRE_CLOSE(-1.79301, reader.lod_scores[0], 0.0001);
    BOOST_REQUIRE_CLOSE(-1.79301, reader.lod_scores[1], 0.0001);
    BOOST_REQUIRE_CLOSE(0.301874, reader.lod_scores[2], 0.0001);
    BOOST_REQUIRE_CLOSE(0.301874, reader.lod_scores[3], 0.0001);

    // "3\t126\tA\tT\t0\t2\t2\t2\t2\n" change chrom, 0/2 override check
    BOOST_REQUIRE(reader.update());
    BOOST_REQUIRE_EQUAL(3, reader.chromosome);
    BOOST_REQUIRE_EQUAL(126, reader.position);
    BOOST_REQUIRE_CLOSE(-1.99568, reader.lod_scores[0], 0.0001);
    BOOST_REQUIRE_CLOSE(-1.99568, reader.lod_scores[1], 0.0001);
    BOOST_REQUIRE_CLOSE(-1.99568, reader.lod_scores[2], 0.0001);
    BOOST_REQUIRE_CLOSE(-1.99568, reader.lod_scores[3], 0.0001);

    // eof
    BOOST_REQUIRE(!reader.update());
    BOOST_REQUIRE(!reader.update());
    BOOST_REQUIRE(!reader.update());
}
