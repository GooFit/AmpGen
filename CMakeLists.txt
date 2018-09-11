################################################################################
# Package: AmpGen
################################################################################

# Support standalone build
if(COMMAND gaudi_subdir)
    gaudi_subdir(AmpGen v1r1)

    find_package(ROOT COMPONENTS RIO Hist Matrix Graf Minuit2 Tree MathMore MathCore Physics)

    include_directories(${CMAKE_CURRENT_SOURCE_DIR})

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -Wno-unused-parameter -ldl")

    gaudi_add_library(AmpGen
                      src/*.cpp
                      src/Lineshapes/*.cpp
                      NO_PUBLIC_HEADERS
              INCLUDE_DIRS ROOT
              LINK_LIBRARIES ROOT)

    foreach(app Generator Fitter ConvertToSourceCode)
      gaudi_add_executable(${app} apps/${app}.cpp
                           LINK_LIBRARIES AmpGen)
    endforeach()
else()
    cmake_minimum_required(VERSION 3.3)
    project(AmpGen LANGUAGES CXX VERSION 1.1)
    include(Standalone.cmake)
endif()
