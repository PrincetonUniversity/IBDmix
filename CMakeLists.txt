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

include(FetchContent)
FetchContent_Declare(
  cli11
  GIT_REPOSITORY https://github.com/CLIUtils/CLI11
  GIT_TAG        v2.2.0
)

FetchContent_MakeAvailable(cli11)

# Adding simple libraries
# set(CLI11_TESTING OFF)
# mark_as_advanced(CLI_CXX_STD CLI_EXAMPLES CLI_SINGLE_FILE CLI_SINGLE_FILE_TESTS CLI_TESTING)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTS)
  message(STATUS "TESTING NOW")
  enable_testing()
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
