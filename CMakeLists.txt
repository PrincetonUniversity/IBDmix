cmake_minimum_required(VERSION 3.13)

project(IBDmix VERSION 0.2
              DESCRIPTION "Estimate introgression by IBD"
              LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
option(BUILD_TESTS "Build unit tests" OFF)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif()

# Get the git command
find_package(Git QUIET)

if(GIT_FOUND AND EXISTS "${IBDmix_SOURCE_DIR}/.git")
    # Update submodules as needed
    option(IBDMIX_SUBMODULE "Check submodules during build" ON)
    if(IBDMIX_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init
            WORKING_DIRECTORY ${IBDmix_SOURCE_DIR}
            RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()

    execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
        WORKING_DIRECTORY "${IBDmix_SOURCE_DIR}"
        OUTPUT_VARIABLE IBDMIX_GIT_VERSION
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE)

else()
    set(IBDMIX_GIT_VERSION "unknown")
endif()

if( NOT EXISTS "${IBDmix_SOURCE_DIR}/extern/CLI11/CMakeLists.txt")
    message(FATAL_ERROR "The submodules were not downloaded! IBDMIX_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif()

# Adding simple libraries
set(CLI11_TESTING OFF)
add_subdirectory("extern/CLI11")
mark_as_advanced(CLI_CXX_STD CLI_EXAMPLES CLI_SINGLE_FILE CLI_SINGLE_FILE_TESTS CLI_TESTING)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTS)
  message(STATUS "TESTING NOW")
  enable_testing()
  include(GoogleTest)
  add_subdirectory(tests)

  # Download and unpack googletest at configure time
  configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
  if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
  if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
  endif()

  # Prevent overriding the parent project's compiler/linker
  # settings on Windows
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

  # Add googletest directly to our build. This defines
  # the gtest and gtest_main targets.
  add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                   ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                   EXCLUDE_FROM_ALL)
endif()

add_subdirectory(src)
