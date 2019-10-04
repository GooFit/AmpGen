################################################################################
# Package: AmpGen
################################################################################

# Support standalone build
if(COMMAND gaudi_subdir)
    gaudi_subdir(AmpGen v1r2)

    find_package(ROOT COMPONENTS RIO Hist Matrix Graf Minuit2 Tree MathMore MathCore Physics)
    find_package(TBB REQUIRED)
    find_package(GSL REQUIRED)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR})

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -Wno-unused-parameter -ldl")
    set(AmpGen_VERSION_MAJOR "1")
    set(AmpGen_VERSION_MINOR "1")
    set(AMPGEN_CXX ${CMAKE_CXX_COMPILER}  CACHE FILEPATH "This should be the path to compiler (use which c++ for macOS)" )
    configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/AmpGen/Version.h.in" "${CMAKE_CURRENT_SOURCE_DIR}/AmpGen/Version.h")
    gaudi_add_library(AmpGen
                      src/*.cpp
                      src/Lineshapes/*.cpp
                      NO_PUBLIC_HEADERS
              INCLUDE_DIRS ROOT 
              LINK_LIBRARIES ROOT TBB GSL)
    target_compile_definitions(AmpGen PUBLIC "AMPGEN_CXX=\"${AMPGEN_CXX}\"")
    foreach(app Generator Fitter ConvertToSourceCode)
      gaudi_add_executable(${app} apps/${app}.cpp
                           LINK_LIBRARIES AmpGen)
    endforeach()
else()
    cmake_minimum_required(VERSION 3.3)
    project(AmpGen LANGUAGES CXX VERSION 1.2)
    include(Standalone.cmake)
endif()
