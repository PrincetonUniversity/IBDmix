add_library(vcf_file STATIC vcf_file.cc ${IBDmix_SOURCE_DIR}/include/IBDmix/vcf_file.h)
target_include_directories(vcf_file PUBLIC ../include)

add_executable(generate_gt generate_gt.cc)
target_include_directories(generate_gt PUBLIC ../include)
target_link_libraries(generate_gt vcf_file CLI11::CLI11)

add_library(ibd_stack STATIC IBD_Stack.cc ${IBDmix_SOURCE_DIR}/include/IBDmix/IBD_Stack.h)
target_include_directories(ibd_stack PUBLIC ../include)

add_library(mask_reader STATIC Mask_Reader.cc ${IBDmix_SOURCE_DIR}/include/IBDmix/Mask_Reader.h)
target_include_directories(mask_reader PUBLIC ../include)

add_library(lod_calculator STATIC lod_calculator.cc ${IBDmix_SOURCE_DIR}/include/IBDmix/lod_calculator.h)
target_include_directories(lod_calculator PUBLIC ../include)

add_library(sample_mapper STATIC Sample_Mapper.cc ${IBDmix_SOURCE_DIR}/include/IBDmix/Sample_Mapper.h)
target_include_directories(sample_mapper PUBLIC ../include)

add_library(genotype_reader STATIC Genotype_Reader.cc ${IBDmix_SOURCE_DIR}/include/IBDmix/Genotype_Reader.h)
target_include_directories(genotype_reader PUBLIC ../include)
target_link_libraries(genotype_reader
    mask_reader sample_mapper lod_calculator)

add_library(recorders STATIC Segment_Recorders.cc ${IBDmix_SOURCE_DIR}/include/IBDmix/Segment_Recorders.h)
target_include_directories(recorders PUBLIC ../include)
target_link_libraries(recorders
    ibd_stack genotype_reader)

add_library(ibd_segment IBD_Segment.cc
    ${IBDmix_SOURCE_DIR}/include/IBDmix/IBD_Segment.h)
target_include_directories(ibd_segment PUBLIC ../include)
target_link_libraries(ibd_segment
    genotype_reader ibd_stack recorders)

add_library(ibd_collection IBD_Collection.cc
    ${IBDmix_SOURCE_DIR}/include/IBDmix/IBD_Collection.h)
target_include_directories(ibd_collection PUBLIC ../include)
target_link_libraries(ibd_collection
    genotype_reader ibd_segment ibd_stack)

add_executable(ibdmix main.cc)
target_include_directories(ibdmix PUBLIC ../include)
target_link_libraries(ibdmix
    ibd_collection genotype_reader ibd_stack CLI11::CLI11)

add_executable(gt_lods tabulate_lods.cc)
target_include_directories(gt_lods PUBLIC ../include)
target_link_libraries(gt_lods
    genotype_reader CLI11::CLI11)

install(
  TARGETS
    gt_lods
    ibdmix
    generate_gt
  DESTINATION
    bin
  )
