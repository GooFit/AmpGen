set(AMPGEN_CXX ${CMAKE_CXX_COMPILER}  CACHE FILEPATH "This should be the path to compiler (use which c++ for macOS)" )

file(GLOB_RECURSE AMPGEN_SRC src/*)
file(GLOB_RECURSE AMPGEN_HDR AmpGen/*)

if( NOT "${CMAKE_CXX_STANDARD}" ) 
  set(CMAKE_CXX_STANDARD 14) 
endif() 

SET(USE_OPENMP TRUE CACHE BOOL "USE_OPENMP")

message(STATUS "USE_OPENMP = ${USE_OPENMP}")

set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_TEST_OUTPUT_DIRECTORY    "${CMAKE_BINARY_DIR}/bin/test")

include(CMakeDependentOption)
include(CMakePrintHelpers)

option(AMPGEN_DEBUG "AmpGen Debug printout")
option(AMPGEN_TRACE "AmpGen Trace printout")

configure_file ("${PROJECT_SOURCE_DIR}/AmpGen/Version.h.in" "${CMAKE_BINARY_DIR}/AmpGenVersion.h")

add_library(AmpGen SHARED ${AMPGEN_SRC} ${AMPGEN_HDR})

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

message( STATUS "ROOT_INCLUDE_DIRS = ${ROOT_INCLUDE_DIRS}")

target_include_directories(AmpGen PUBLIC "${CMAKE_SOURCE_DIR}")
target_include_directories(AmpGen PUBLIC "${CMAKE_BINARY_DIR}")
target_include_directories(AmpGen SYSTEM PUBLIC "${ROOT_INCLUDE_DIRS}")

target_link_libraries(AmpGen PUBLIC ${ROOT_LIBRARIES} ${CMAKE_DL_LIBS})

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

target_compile_definitions(AmpGen
  PUBLIC
  "AMPGENROOT_CMAKE=\"${CMAKE_BINARY_DIR}/bin\""
  "AMPGEN_CXX=\"${AMPGEN_CXX}\""
  $<$<BOOL:${AMPGEN_DEBUG}>:DEBUGLEVEL=1>
  $<$<BOOL:${AMPGEN_TRACE}>:TRACELEVEL=1>)

target_compile_options(AmpGen
  PUBLIC
  -Wall -Wextra -Wpedantic -g3
  -Wno-unused-parameter
  -Wno-unknown-pragmas
  -Wnon-virtual-dtor
  -Wno-overloaded-virtual
  -march=native
  $<$<CONFIG:Release>:-Ofast>)

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lm -lstdc++")
else()
  target_compile_options(AmpGen PUBLIC -Wno-suggest-override)
endif()

file(GLOB_RECURSE applications apps/*.cpp )
file(GLOB_RECURSE examples examples/*.cpp )

foreach( file ${applications} )
  get_filename_component( Executable ${file} NAME_WE )
  cmake_print_variables(Executable)
  add_executable(${Executable} ${file})
  target_compile_options(${Executable} PUBLIC -g3 -Ofast)
  target_link_libraries(${Executable} PUBLIC AmpGen)
endforeach()

foreach( file ${examples} )
  get_filename_component( Executable ${file} NAME_WE )
  cmake_print_variables(Executable)
  add_executable(${Executable} ${file})
  target_link_libraries(${Executable} PUBLIC AmpGen)
endforeach()

file(GLOB_RECURSE options_files options/*.*)
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_BINARY_DIR}/bin")
foreach(file ${options_files})
  get_filename_component(OptionFile "${file}" NAME)
  cmake_print_variables(OptionFile)
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink "${file}" "${CMAKE_BINARY_DIR}/bin/${OptionFile}")
endforeach()

enable_testing()
set(Boost_NO_BOOST_CMAKE ON)

find_package(Boost 1.67.0 COMPONENTS unit_test_framework)
if ( Boost_FOUND )
  include_directories (${Boost_INCLUDE_DIRS})
  file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} test/*.cpp)
  foreach(testSrc ${TEST_SRCS})
    get_filename_component(testName ${testSrc} NAME_WE)
    add_executable(${testName} ${testSrc})
    set_target_properties(${testName} 
      PROPERTIES 
      RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_TEST_OUTPUT_DIRECTORY}" 
      RUNTIME_OUTPUT_DIRECTORY_DEBUG   "${CMAKE_TEST_OUTPUT_DIRECTORY}" 
      RUNTIME_OUTPUT_DIRECTORY         "${CMAKE_TEST_OUTPUT_DIRECTORY}" 
      EXECUTABLE_OUTPUT_DIRECTORY      "${CMAKE_TEST_OUTPUT_DIRECTORY}" 
    )
    target_link_libraries(${testName} PUBLIC ${Boost_LIBRARIES} AmpGen)
    add_test(NAME ${testName} WORKING_DIRECTORY ${CMAKE_TEST_OUTPUT_DIRECTORY} COMMAND ${CMAKE_TEST_OUTPUT_DIRECTORY}/${testName} ) 
  endforeach(testSrc)
else()
  message( WARNING "Warning: Boost (version >= 1.67.0) required to build unit tests\n")
endif()

