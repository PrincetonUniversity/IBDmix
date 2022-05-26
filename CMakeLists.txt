cmake_minimum_required(VERSION 3.14)

project(IBDmix VERSION 1.0.0
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
    message(STATUS "Downloading dependencies...")
    if( NOT EXISTS "${IBDmix_SOURCE_DIR}/extern/CLI11/CMakeLists.txt")
        file(DOWNLOAD
            https://github.com/CLIUtils/CLI11/archive/b9a2f320b912aae623c4acddea06cba399b98654.zip
            "${IBDmix_SOURCE_DIR}/extern/CLI11.zip"
            SHOW_PROGRESS
            EXPECTED_HASH MD5=222245b96453401aa366328e398ae1ad)
        message(STATUS "Extracting CLI11")
        execute_process(COMMAND ${CMAKE_COMMAND} -E tar -xf
            "${IBDmix_SOURCE_DIR}/extern/CLI11.zip"
            WORKING_DIRECTORY
            "${IBDmix_SOURCE_DIR}/extern/"
            )
        file(RENAME
            "${IBDmix_SOURCE_DIR}/extern/CLI11-b9a2f320b912aae623c4acddea06cba399b98654/"
            "${IBDmix_SOURCE_DIR}/extern/CLI11/"
            )
        file(REMOVE
            "${IBDmix_SOURCE_DIR}/extern/CLI11.zip"
            )
    endif()
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
  include(FetchContent)
  FetchContent_Declare(
      googletest
      URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
      )
  # For Windows: Prevent overriding the parent project's compiler/linker settings
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)
  include(GoogleTest)
  add_subdirectory(tests)

endif()

add_subdirectory(src)
