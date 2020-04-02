add_library(vcf_file STATIC vcf_file.cc ${IBDmix_SOURCE_DIR}/include/IBDmix/vcf_file.h)
target_include_directories(vcf_file PUBLIC ../include)

add_executable(generate_gt generate_gt.cc)
target_include_directories(generate_gt PUBLIC ../include)
target_link_libraries(generate_gt vcf_file)

add_library(ibd_stack STATIC IBD_Stack.cc ${IBDmix_SOURCE_DIR}/include/IBDmix/IBD_Stack.h)
target_include_directories(ibd_stack PUBLIC ../include)

add_library(genotype_reader STATIC Genotype_Reader.cc ${IBDmix_SOURCE_DIR}/include/IBDmix/Genotype_Reader.h)
target_include_directories(genotype_reader PUBLIC ../include)

add_library(ibd_segment IBD_Segment.cc
    ${IBDmix_SOURCE_DIR}/include/IBDmix/IBD_Segment.h)
target_include_directories(ibd_segment PUBLIC ../include)
target_link_libraries(ibd_segment 
    genotype_reader ibd_stack)

add_library(ibd_collection IBD_Collection.cc
    ${IBDmix_SOURCE_DIR}/include/IBDmix/IBD_Collection.h)
target_include_directories(ibd_collection PUBLIC ../include)
target_link_libraries(ibd_collection 
    genotype_reader ibd_segment ibd_stack)

add_executable(ibdmix main.cc)
target_include_directories(ibdmix PUBLIC ../include)
target_link_libraries(ibdmix
    ibd_collection genotype_reader ibd_stack)
