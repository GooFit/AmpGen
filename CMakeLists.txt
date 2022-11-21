################################################################################
# Package: AmpGen
################################################################################

cmake_minimum_required(VERSION 3.9...3.15)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project(AmpGen LANGUAGES CXX VERSION 2.3)

set(AMPGEN_CXX ${CMAKE_CXX_COMPILER} CACHE FILEPATH "This should be the path to compiler (use which c++ for macOS)" )

file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
  message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

include(${PROJECT_SOURCE_DIR}/cmake/FindAVX.cmake)
CHECK_FOR_AVX()

set(CMAKE_CXX_STANDARD "17"    CACHE STRING "CMAKE_CXX_STANDARD")

if( AVX2_FOUND )
  set(USE_SIMD           "AVX2d" CACHE STRING "USE_SIMD")       # AVX instruction set + precision to use 
  set(USE_MVEC           TRUE    CACHE BOOL   "USE_MVEC")       # flag to use vector math library mvec 
else()
  set(USE_SIMD           "0"     CACHE STRING "USE_SIMD")       # AVX instruction set + precision to use 
  set(USE_MVEC           FALSE   CACHE BOOL   "USE_MVEC")       # flag to use vector math library mvec 
endif()

set(USE_OPENMP         TRUE    CACHE BOOL   "USE_OPENMP")     # flag to use openmp for threading 
set(ENABLE_FITTING     TRUE    CACHE BOOL   "ENABLE_FITTING") # flag to enable fitting (requires Minuit2) 
set(ENABLE_TESTS       TRUE    CACHE BOOL   "ENABLE_TESTS")   # flag to build unit tests
set(BUILTIN_GSL        TRUE    CACHE BOOL   "BUILTIN_GSL")    # flag to use ROOT's builtin GSL 
set(ENABLE_INSTALL     TRUE    CACHE BOOL   "ENABLE_INSTALL") # flag to enable installs (switched off to make Travis happier) 
set(AMPGENROOT         "${PROJECT_SOURCE_DIR}" CACHE STRING "AMPGENROOT")
set(AMPGENROOT_CMAKE   "${CMAKE_BINARY_DIR}/bin" CACHE STRING "AMPGENROOT_CMAKE")
set(USE_ROOT           TRUE    CACHE BOOL   "Flag to use root (required for fitting, plotting & I/O)")
set(USE_TBB           FALSE    CACHE BOOL   "USE_TBB")     # flag to use TBB (only used internally by ROOT, but linkage seems unreliable)
set(EXEC_EXTENSION     "" CACHE STRING "EXEC_EXTENSION")
set(PACKED_CMAKE_PREFIX_PATH "" CACHE STRING "PACKED_CMAKE_PREFIX_PATH")
set(CMAKE_SHARED_LIBRARY_PREFIX "") # such that we can have libAmpGen as a target AND AmpGen i.e. the application

if( NOT PACKED_CMAKE_PREFIX_PATH STREQUAL "" )
  string( REPLACE "|" " " CMAKE_PREFIX_PATH ${PACKED_CMAKE_PREFIX_PATH})
  separate_arguments(PACKED_CMAKE_PREFIX_PATH)
endif()


file(GLOB_RECURSE AMPGEN_SRC src/*)
file(GLOB_RECURSE AMPGEN_HDR AmpGen/*)

if( NOT ENABLE_FITTING )
  list(REMOVE_ITEM AMPGEN_SRC ${PROJECT_SOURCE_DIR}/src/Minimiser.cpp)
  list(REMOVE_ITEM AMPGEN_SRC ${PROJECT_SOURCE_DIR}/src/FitResult.cpp)
  list(REMOVE_ITEM AMPGEN_HDR ${PROJECT_SOURCE_DIR}/AmpGen/Minimiser.h)
  list(REMOVE_ITEM AMPGEN_HDR ${PROJECT_SOURCE_DIR}/AmpGen/FitResult.h)
endif()

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

add_library(lib${PROJECT_NAME} SHARED ${AMPGEN_SRC} ${AMPGEN_HDR})
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS lib${PROJECT_NAME})

if( USE_ROOT AND NOT DEFINED ROOT_LIBRARIES )
  if(DEFINED ENV{ROOTSYS})
    list(APPEND CMAKE_MODULE_PATH "$ENV{ROOTSYS}/etc/cmake/")
  endif()
  find_package(ROOT COMPONENTS Matrix MathCore Gpad Tree Graf)
endif() 

if ( NOT DEFINED GSL_LIBRARIES )
  find_package(GSL)
endif()

if( USE_OPENMP )
  find_package(OpenMP)
endif()

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

target_include_directories(libAmpGen PUBLIC $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}> $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}> )
target_include_directories(libAmpGen PUBLIC ${ROOT_INCLUDE_DIRS} ${GSL_INCLUDE_DIR} ) 
target_link_libraries(     libAmpGen PUBLIC -lm ${ROOT_LIBRARIES} ${CMAKE_DL_LIBS} ${GSL_LIBRARIES} )

find_library(libmvec mvec)

if ( USE_MVEC AND libmvec  )
  message( DEBUG "Using libmvec for vectorised math operations")
  target_link_libraries(libAmpGen PUBLIC mvec)
else()
  message( DEBUG "libmvec not found, with use scalar math where necessary.")
endif()

if( ENABLE_FITTING AND USE_ROOT )
  if( ( NOT TARGET ROOT::Minuit2 AND NOT TARGET Minuit2 ) OR "${extern_minuit2}" )
    message( STATUS "Use external Minuit2")
    add_subdirectory("extern/Minuit2")
    set_target_properties(Minuit2     PROPERTIES FOLDER extern)
    target_compile_options(Minuit2 PUBLIC -fPIC -Wno-suggest-override)
    set_target_properties(Minuit2Math PROPERTIES FOLDER extern)
    add_library(ROOT::Minuit2 ALIAS Minuit2)
    target_include_directories( libAmpGen PUBLIC "${CMAKE_SOURCE_DIR}/extern/Minuit2/inc/")
  else()
    message( STATUS "Use ROOT::Minuit2")
  endif()
  if ( TARGET Minuit2 AND NOT TARGET ROOT::Minuit2 )
    find_package( ROOT CONFIG REQUIRED COMPONENTS Minuit2)
    add_library(ROOT::Minuit2 ALIAS Minuit2)
  endif()
  target_link_libraries(libAmpGen PUBLIC ROOT::Minuit2 )
endif()

if( USE_OPENMP )
  list( GET OpenMP_CXX_LIBRARIES 0 openmp_lib) 
  message(DEBUG "OpenMP_CXX_LIBRARIES: ${openmp_lib} CXX_FLAGS: ${OpenMP_CXX_FLAGS}")
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
    target_link_libraries(libAmpGen PUBLIC OpenMP::OpenMP_CXX)
    get_filename_component( libpath ${openmp_lib} DIRECTORY)
    message(DEBUG "OpenMP libpath: ${libpath}")
    target_compile_definitions(libAmpGen PRIVATE "AMPGEN_OPENMP_FLAGS=\"-L${libpath} ${OpenMP_CXX_FLAGS}\"")
  else()
    message(STATUS "OpenMP not found for CXX")
  endif()
else()
  target_compile_definitions(libAmpGen PRIVATE "AMPGEN_OPENMP_FLAGS=\"\"")
endif()

if( DEFINED TBB_DIR )
  list(APPEND CMAKE_MODULE_PATH "${TBB_DIR}")
  find_package(TBB REQUIRED)
  target_link_libraries(libAmpGen PUBLIC TBB::tbb)

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
  target_link_libraries(libAmpGen PUBLIC ${XROOTD_LIB})
endif()

target_compile_definitions(libAmpGen PUBLIC
  "USE_MVEC=$<BOOL:${USE_MVEC} AND ${libmvec}>"
  "ENABLE_FITTING=${ENABLE_FITTING}"
)

target_compile_definitions(libAmpGen PRIVATE
  "AMPGENROOT_CMAKE=\"${AMPGENROOT_CMAKE}\""
  "AMPGENROOT=\"${AMPGENROOT}\""
  "AMPGEN_CXX=\"${AMPGEN_CXX}\""
  $<$<BOOL:${AMPGEN_DEBUG}>:DEBUGLEVEL=1>
  $<$<BOOL:${AMPGEN_TRACE}>:TRACELEVEL=1>)

target_compile_options(libAmpGen
  INTERFACE
  -Wall -Wextra -Wpedantic -g3
  -Wno-unused-parameter
  -Wno-unknown-pragmas
  $<$<CONFIG:Release>:-O3>)

if ( ${USE_SIMD} MATCHES "AVX2d" )
  message(STATUS "Enabling AVX2 [double precision]")
  target_compile_definitions(libAmpGen PUBLIC "ENABLE_AVX=1" "INSTRUCTION_SET=2" ) 
  target_compile_options(libAmpGen PUBLIC -ftree-vectorize -mfma -mavx2 -ffast-math -fno-finite-math-only -DHAVE_AVX2_INSTRUCTIONS)
elseif ( ${USE_SIMD} MATCHES "AVX2f" )
  message(STATUS "Enabling AVX2 [single precision]")
  target_compile_definitions(libAmpGen PUBLIC "ENABLE_AVX=1" "INSTRUCTION_SET=1") 
  target_compile_options(libAmpGen PUBLIC -ftree-vectorize -mfma -mavx2 -ffast-math -fno-finite-math-only -DHAVE_AVX2_INSTRUCTIONS)
elseif ( ${USE_SIMD} MATCHES "AVX512d" )
  message(STATUS "Enabling AVX512 [double precision]")
  target_compile_definitions(libAmpGen PUBLIC "ENABLE_AVX=1" "INSTRUCTION_SET=3") 
  message(WARNING "AVX512 implementation not tested: proceed with caution")
  target_compile_options(libAmpGen PUBLIC -ftree-vectorize -mfma -mavx512f -ffast-math -DHAVE_AVX512_INSTRUCTIONS)
elseif( ${USE_SIMD} MATCHES "ARM128d" )
  message(STATUS "Enabling ARM128d [double precision]")
  target_compile_definitions(libAmpGen PUBLIC "ENABLE_AVX=1" "INSTRUCTION_SET=10") 
  message(WARNING "ARM128 (double precision) implementation not tested: proceed with caution")
  target_compile_options(libAmpGen PUBLIC -march=armv8.5-a  )
else()
  target_compile_definitions(libAmpGen PUBLIC "ENABLE_AVX=0" "INSTRUCTION_SET=0")
  target_compile_options(libAmpGen PUBLIC -march=native)
  message(STATUS "SIMD disabled, resorting to scalar build : ${USE_SIMD}")
endif()

if(${CMAKE_CXX_COMPILER_ID} MATCHES "AppleClang" OR 
    ${CMAKE_CXX_COMPILER_ID} MATCHES "Clang" )
  target_compile_options(libAmpGen PUBLIC -mfma)
endif()

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "AppleClang" )
  target_link_libraries(libAmpGen PUBLIC stdc++ )
  message(STATUS "Using OSX specific flags: -lm -lstdc++ -lSystem")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lm -lstdc++ -lSystem")
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  target_link_libraries(libAmpGen PUBLIC stdc++)
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lm -lstdc++")
else()
  target_compile_options(libAmpGen PUBLIC -Wno-suggest-override)
endif()

# build the applications and examples 

if ( ENABLE_FITTING )
  file(GLOB_RECURSE applications apps/*.cpp examples/*.cpp )
else() 
  set(applications apps/AmpGen.cpp)
endif()

foreach( file ${applications} )
  get_filename_component( Executable ${file} NAME_WE )
  add_executable(${Executable}${EXEC_EXTENSION} ${file})
  target_compile_options(${Executable}${EXEC_EXTENSION} PUBLIC -g3 -O3)
  target_link_libraries(${Executable}${EXEC_EXTENSION} PUBLIC libAmpGen  )
  target_include_directories(${Executable}${EXEC_EXTENSION} PRIVATE ${ROOT_INCLUDE_DIRS})
endforeach()

# symlink the options files (so as there is the mass_width etc.) to the build directory, and also the release notes 
file(GLOB_RECURSE options_files options/*.*)
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_BINARY_DIR}/bin")
foreach(file ${options_files})
  get_filename_component(OptionFile "${file}" NAME)
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink "${file}" "${CMAKE_BINARY_DIR}/bin/${OptionFile}")
endforeach()
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink "${CMAKE_SOURCE_DIR}/doc/release.notes" "${CMAKE_BINARY_DIR}/bin/release.notes")

if( ENABLE_TESTS )
  enable_testing()
  set(Boost_NO_BOOST_CMAKE ON)
  add_subdirectory(test)
endif()


# stuff related to export, install and version tagging

include(CMakePackageConfigHelpers)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/AmpGenVersion.cmake VERSION ${PACKAGE_VERSION} COMPATIBILITY AnyNewerVersion)

set(INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_INCLUDEDIR})
set(LIB_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR})
set(DATA_INSTALL_DIR ${CMAKE_INSTALL_DATADIR}/AmpGen)

configure_package_config_file(AmpGenConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/AmpGenConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_DATADIR}/AmpGen/cmake
    PATH_VARS INCLUDE_INSTALL_DIR LIB_INSTALL_DIR DATA_INSTALL_DIR
    )

if( ENABLE_INSTALL ) 
  install( FILES ${CMAKE_CURRENT_BINARY_DIR}/AmpGenConfig.cmake
                 ${CMAKE_CURRENT_BINARY_DIR}/AmpGenVersion.cmake
                 DESTINATION ${CMAKE_INSTALL_DATADIR}/AmpGen/cmake)
  
  install(TARGETS libAmpGen
    EXPORT AmpGenTargets
    LIBRARY DESTINATION lib
    PUBLIC_HEADER DESTINATION AmpGen
  )
  
  install(TARGETS AmpGen${EXEC_EXTENSION}
    DESTINATION bin
  )
  
  install(EXPORT "AmpGenTargets"
          NAMESPACE "AmpGen::"
          DESTINATION ${CMAKE_INSTALL_DATADIR}/AmpGen/cmake
  )
  
  export( TARGETS libAmpGen NAMESPACE AmpGen:: FILE AmpGenTargets.cmake )
  set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
  export(PACKAGE AmpGen)
endif()
