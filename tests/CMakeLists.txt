# add in macros for testing
macro(package_add_test TESTNAME FILES LIBRARIES)
    # create an exectuable in which the tests will be stored
    add_executable(${TESTNAME} ${FILES})
    # link the Google test infrastructure, mocking library,
    target_link_libraries(${TESTNAME} ${LIBRARIES})
    target_link_libraries(${TESTNAME} gtest gmock gtest_main)
    # gtest_discover_tests replaces gtest_add_tests,
    gtest_discover_tests(${TESTNAME}
        WORKING_DIRECTORY ${PROJECT_DIR}
        PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_DIR}"
        )
    set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
    target_include_directories(${TESTNAME} PUBLIC include)
    add_test(NAME, ${TESTNAME} COMMAND ${TESTNAME})
endmacro()

package_add_test(ibd_stack_test test_IBD_Stack.cc ibd_stack)
package_add_test(ibd_segment_test test_IBD_Segment.cc ibd_segment)
package_add_test(mask_reader_test test_Mask_Reader.cc mask_reader)
package_add_test(sample_mapper_test test_Sample_Mapper.cc sample_mapper)
package_add_test(genotype_reader_test test_Genotype_Reader.cc genotype_reader)
package_add_test(vcf_file_test test_vcf_file.cc vcf_file)
