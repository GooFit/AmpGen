set(AMPGEN_CXX ${CMAKE_CXX_COMPILER} CACHE FILEPATH "This should be the path to compiler (use which c++ for macOS)" )

file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

file(GLOB_RECURSE AMPGEN_SRC src/*)
file(GLOB_RECURSE AMPGEN_HDR AmpGen/*)

if( NOT "${CMAKE_CXX_STANDARD}" ) 
  set(CMAKE_CXX_STANDARD 17) 
endif() 

SET(USE_OPENMP TRUE    CACHE BOOL   "USE_OPENMP") # flag to use openmp for threading 
SET(USE_SIMD   "AVX2d" CACHE STRING "USE_SIMD")   # AVX instruction set + precision to use 
SET(USE_MVEC   TRUE    CACHE BOOL   "USE_MVEC")   # flag to use vector math library mvec 

set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY         "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY         "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY         "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG   "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_TEST_OUTPUT_DIRECTORY            "${CMAKE_BINARY_DIR}/bin/test")

include(CMakeDependentOption)
include(CMakePrintHelpers)
include(GNUInstallDirs)

option(AMPGEN_DEBUG "AmpGen Debug printout")
option(AMPGEN_TRACE "AmpGen Trace printout")


configure_file ("${PROJECT_SOURCE_DIR}/AmpGen/Version.h.in" "${CMAKE_BINARY_DIR}/AmpGenVersion.h")

add_library(${PROJECT_NAME} SHARED ${AMPGEN_SRC} ${AMPGEN_HDR})
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

if(DEFINED ENV{ROOTSYS})
  list(APPEND CMAKE_MODULE_PATH "$ENV{ROOTSYS}/etc/cmake/")
endif()

find_package(ROOT CONFIG REQUIRED COMPONENTS Matrix MathMore MathCore Gpad Tree Graf)

if( USE_OPENMP )
  find_package(OpenMP)
endif()

cmake_print_variables(CMAKE_SOURCE_DIR)

# Default build type from the Kitware Blog
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
    STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

target_include_directories(AmpGen PUBLIC $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}>
                                         $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}> )

target_include_directories(AmpGen SYSTEM PUBLIC "${ROOT_INCLUDE_DIRS}")

target_link_libraries(AmpGen PUBLIC -lm ${ROOT_LIBRARIES} ${CMAKE_DL_LIBS} )

find_library(libmvec mvec)

if ( USE_MVEC AND libmvec  )
  message( STATUS "Using libmvec for vectorised math operations")
  target_link_libraries(AmpGen PUBLIC mvec)
  target_compile_definitions(AmpGen PUBLIC "USE_MVEC=1") 
else()
  message( STATUS "Warning: libmvec not found, with use scalar math where necessary.")
  target_compile_definitions(AmpGen PUBLIC "USE_MVEC=0") 
endif()



if( ( NOT TARGET ROOT::Minuit2 AND NOT TARGET Minuit2 ) OR "${extern_minuit2}" )
  message( STATUS "Use external Minuit2")
  add_subdirectory("extern/Minuit2")
  set_target_properties(Minuit2     PROPERTIES FOLDER extern)
  target_compile_options(Minuit2 PUBLIC -fPIC -Wno-suggest-override)
  set_target_properties(Minuit2Math PROPERTIES FOLDER extern)
  add_library(ROOT::Minuit2 ALIAS Minuit2)
  target_include_directories( AmpGen PUBLIC "${CMAKE_SOURCE_DIR}/extern/Minuit2/inc/")
else()
  message( STATUS "Use ROOT::Minuit2")
endif()
if ( TARGET Minuit2 AND NOT TARGET ROOT::Minuit2 )
  find_package( ROOT CONFIG REQUIRED COMPONENTS Minuit2)
  add_library(ROOT::Minuit2 ALIAS Minuit2)
endif()

target_link_libraries(AmpGen PUBLIC ROOT::Minuit2 )

if( USE_OPENMP )
  if(OpenMP_FOUND OR OpenMP_CXX_FOUND)
    if(NOT TARGET OpenMP::OpenMP_CXX)
      add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
      set_property(TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
      set_property(TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_LINK_LIBRARIES  ${OpenMP_CXX_FLAGS})
      if(CMAKE_VERSION VERSION_LESS 3.4)
        set_property(TARGET OpenMP::OpenMP_CXX APPEND PROPERTY INTERFACE_LINK_LIBRARIES -pthread)
      else()
        find_package(Threads REQUIRED)
        set_property(TARGET OpenMP::OpenMP_CXX APPEND PROPERTY INTERFACE_LINK_LIBRARIES Threads::Threads)
      endif()
    endif()
    target_link_libraries(AmpGen PUBLIC OpenMP::OpenMP_CXX)
  else()
    message(STATUS "OpenMP not found for CXX, you might have forgotten lb-run ROOT bash or CXX=`which g++` in CERN stack")
  endif()
endif()

set(RAPIDSIM_DATA "")

if( "${USE_RAPIDSIM}" ) 
  include(${CMAKE_ROOT}/Modules/FetchContent.cmake) 
  FetchContent_Declare( RapidSim GIT_REPOSITORY https://github.com/gcowan/RapidSim/ )
  if(NOT RapidSim_POPULATED)
    message("Fetching RapidSim from: https://github.com/gcowan/RapidSim/")
    FetchContent_Populate(RapidSim)
  endif()
  FetchContent_GetProperties(RapidSim)
  set(RAPIDSIM_DATA "${CMAKE_BINARY_DIR}/_deps/rapidsim-src/rootfiles/" )
endif()

if( RAPIDSIM_DATA )
  message("Set RAPIDSIM_DATA = ${RAPIDSIM_DATA}")
endif()

# Default to XROOTD only if on CMT system. Can be overridden with -DAMPGEN_XROOTD=ON
if(DEFINED ENV{CMTCONFIG})
  set(AMPGEN_XROOTD_DEFAULT ON)
else()
  set(AMPGEN_XROOTD_DEFAULT OFF)
endif()

cmake_dependent_option(AMPGEN_XROOTD "Turn on XROOTD discovery" ON "AMPGEN_XROOTD_DEFAULT" OFF)

if(AMPGEN_XROOTD)
  find_library(XROOTD_LIB NAMES libXrdCl.so
    HINTS "/cvmfs/lhcb.cern.ch/lib/lcg/releases/LCG_89/xrootd/4.6.0/$ENV{CMTCONFIG}/lib64")
  target_link_libraries(AmpGen PUBLIC ${XROOTD_LIB})
endif()

target_compile_definitions(AmpGen PRIVATE
  "AMPGENROOT_CMAKE=\"${CMAKE_BINARY_DIR}/bin\""
  "AMPGENROOT=\"${PROJECT_SOURCE_DIR}\""
  "AMPGEN_CXX=\"${AMPGEN_CXX}\""
  "USE_OPENMP=\"${USE_OPENMP}\""
  $<$<BOOL:${AMPGEN_DEBUG}>:DEBUGLEVEL=1>
  $<$<BOOL:${AMPGEN_TRACE}>:TRACELEVEL=1>)


target_compile_options(AmpGen
  INTERFACE
  -Wall -Wextra -Wpedantic -g3
  -Wno-unused-parameter
  -Wno-unknown-pragmas
  $<$<CONFIG:Release>:-O3>)

# if( NOT ${USE_SIMD} MATCHES "" ) 
  if ( ${USE_SIMD} MATCHES "AVX2d" )
  message(STATUS "Enabling AVX2 [double precision]")
  target_compile_definitions(AmpGen PUBLIC "ENABLE_AVX=1" "ENABLE_AVX2d=1") 
  target_compile_options(AmpGen PUBLIC -march=native -ftree-vectorize -mavx2 -ffast-math -DHAVE_AVX2_INSTRUCTIONS)
  elseif ( ${USE_SIMD} MATCHES "AVX2f" )
  message(STATUS "Enabling AVX2 [single precision]")
  target_compile_definitions(AmpGen PUBLIC "ENABLE_AVX=1" "ENABLE_AVX2f=1") 
  target_compile_options(AmpGen PUBLIC -march=native -ftree-vectorize -mavx2 -ffast-math -DHAVE_AVX2_INSTRUCTIONS)
  elseif ( ${USE_SIMD} MATCHES "AVX512d" )
  message(STATUS "Enabling AVX2 [double precision]")
  target_compile_definitions(AmpGen PUBLIC "ENABLE_AVX=1" "ENABLE_AVX512=1") 
  target_compile_options(AmpGen PUBLIC -march=native -ftree-vectorize -mavx512f -ffast-math -DHAVE_AVX512_INSTRUCTIONS)
  endif()
  if("${CMAKE_CXX_COMPILER_ID}" MATCHES "AppleClang" OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" )
   target_compile_options(AmpGen PUBLIC -mfma)
  endif()
  # else()
  #   message("SIMD disabled, resorting to scalar build : ${USE_SIMD}")
  # endif()

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "AppleClang" )
  target_link_libraries(AmpGen PUBLIC stdc++ )
  message(STATUS "Using OSX specific flags: -lm -lstdc++ -lSystem")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lm -lstdc++ -lSystem")
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  target_link_libraries(AmpGen PUBLIC stdc++)
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lm -lstdc++")
else()
  target_compile_options(AmpGen PUBLIC -Wno-suggest-override)
endif()

file(GLOB_RECURSE applications apps/*.cpp examples/*.cpp )

foreach( file ${applications} )
  get_filename_component( Executable ${file} NAME_WE )
  #   cmake_print_variables(Executable)
  add_executable(${Executable} ${file})
  target_compile_options(${Executable} PUBLIC -g3 -Ofast)
  target_link_libraries(${Executable} PUBLIC AmpGen  )
endforeach()

file(GLOB_RECURSE options_files options/*.*)
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_BINARY_DIR}/bin")
foreach(file ${options_files})
  get_filename_component(OptionFile "${file}" NAME)
  # cmake_print_variables(OptionFile)
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink "${file}" "${CMAKE_BINARY_DIR}/bin/${OptionFile}")
endforeach()

enable_testing()
set(Boost_NO_BOOST_CMAKE ON)
add_subdirectory(test)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(AmpGenVersion.cmake VERSION ${PACKAGE_VERSION} COMPATIBILITY AnyNewerVersion)
configure_file(AmpGenConfig.cmake.in AmpGenConfig.cmake) #  @ONLY)

export( TARGETS AmpGen NAMESPACE AmpGen:: FILE AmpGenTargets.cmake )
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE AmpGen)
