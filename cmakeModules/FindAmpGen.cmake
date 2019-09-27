#//-------------------------------------------------------------------------
#//
#// Description:
#//      cmake module for finding AmpGen
#//      AmpGen installation location is defined by environment variable $AMPGENROOT
#//
#//      following variables are defined here:
#//      AMPGEN_FOUND       - flag if everything went smoothly
#//      AMPGEN_SOURCE_DIR  - AMPGEN source directory
#//      AMPGEN_INCLUDE_DIR - AMPGEN header directory
#//
#//      Example usage:
#//          find_package(AMPGEN REQUIRED)
#//
#//-------------------------------------------------------------------------


set(AMPGEN_FOUND FALSE)
set(AMPGEN_ERROR_REASON "")

if($ENV{AMPGENROOT} STREQUAL "")
  set(AMPGEN_ERROR_REASON "${AMPGEN_ERROR_REASON} Environment variable AMPGENROOT is not set.")
else()
  set(AMPGEN_FOUND TRUE)

  set(AMPGEN_INCLUDE_DIR "$ENV{AMPGENROOT}/AmpGen")
  if(NOT EXISTS "${AMPGEN_INCLUDE_DIR}")
    set(AMPGEN_FOUND FALSE)
    set(AMPGEN_ERROR_REASON "${AMPGEN_ERROR_REASON} Directory '${AMPGEN_INCLUDE_DIR}' does not exist.")
  endif()

  set(AMPGEN_SOURCE_DIR "$ENV{AMPGENROOT}/src")
  if(NOT EXISTS "${AMPGEN_INCLUDE_DIR}")
    set(AMPGEN_FOUND FALSE)
    set(AMPGEN_ERROR_REASON "${AMPGEN_ERROR_REASON} Directory '${AMPGEN_INCLUDE_DIR}' does not exist.")
  endif()

  if(NOT ROOT_FOUND)
    message(WARNING "ROOT not found. Trying to find it from here...")
    find_package(ROOT CONFIG REQUIRED COMPONENTS Matrix MathMore MathCore Gpad Tree Graf)
  endif()

  if(NOT ROOT_FOUND)
    set(AMPGEN_FOUND FALSE)
    set(AMPGEN_ERROR_REASON "${AMPGEN_ERROR_REASON} ROOT was not found.")
  endif()
endif()

SET(USE_OPENMP TRUE CACHE BOOL "USE_OPENMP")
message(STATUS "USE_OPENMP = ${USE_OPENMP}")
if(USE_OPENMP)
  if(NOT OpenMP_FOUND AND NOT OpenMP_CXX_FOUND)
    message(STATUS "OpenMP not found. Trying to find it from here")
    find_package(OpenMP)
    if(NOT OpenMP_FOUND AND NOT OpenMP_CXX_FOUND)
      set(AMPGEN_FOUND FALSE)
      set(AMPGEN_ERROR_REASON "${AMPGEN_ERROR_REASON} OpenMP was not found.")
    endif()
  endif()
endif()

# make variables changeable
mark_as_advanced(AMPGEN_SOURCE_DIR AMPGEN_INCLUDE_DIR)

# Everything OK. Report result and build
if(AMPGEN_FOUND)
  project(AmpGen LANGUAGES CXX VERSION 1.2)
  message(STATUS "Found AMPGEN in '$ENV{AMPGENROOT}'.")
  message(STATUS "Using AMPGEN source directory '${AMPGEN_SOURCE_DIR}'.")
  message(STATUS "Using AMPGEN include directory '${AMPGEN_INCLUDE_DIR}'.")
  message(STATUS "Will create dynamic library in '${AMPGEN_LIBRARY_DIR}' on demand.")

  # this seems to be important as it switches off the --std=gnu++17 flag which breaks everything
  set(CMAKE_CXX_EXTENSIONS OFF)

  file(GLOB_RECURSE AMPGEN_SRC ${AMPGEN_SOURCE_DIR}/*)
  file(GLOB_RECURSE AMPGEN_HDR ${AMPGEN_INCLUDE_DIR}/*)

  include(CMakePrintHelpers)

  option(AMPGEN_DEBUG "AmpGen Debug printout")
  option(AMPGEN_TRACE "AmpGen Trace printout")

  configure_file("${AMPGEN_INCLUDE_DIR}/Version.h.in" "${CMAKE_BINARY_DIR}/AmpGenVersion.h")

  add_library(AmpGen SHARED ${AMPGEN_SRC} ${AMPGEN_HDR})

  target_include_directories(AmpGen PUBLIC "${CMAKE_BINARY_DIR}")
  target_include_directories(AmpGen PUBLIC "$ENV{AMPGENROOT}")
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

  file(GLOB_RECURSE applications $ENV{AMPGENROOT}/apps/*.cpp )
  file(GLOB_RECURSE examples $ENV{AMPGENROOT}/examples/*.cpp )
  file(GLOB_RECURSE options_files $ENV{AMPGENROOT}/options/*.*)

  foreach( file ${applications} )
    get_filename_component( Executable ${file} NAME_WE )
    cmake_print_variables(Executable)
    add_executable(${Executable} ${file})
    target_compile_options(${Executable} PUBLIC -g3 -Ofast)
    target_link_libraries(${Executable} PUBLIC AmpGen -lMathMore)
  endforeach()

  foreach( file ${examples} )
    get_filename_component( Executable ${file} NAME_WE )
    cmake_print_variables(Executable)
    add_executable(${Executable} ${file})
    target_link_libraries(${Executable} PUBLIC AmpGen -lMathMore)
  endforeach()

  execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_BINARY_DIR}/bin")
  foreach(file ${options_files})
    get_filename_component(OptionFile "${file}" NAME)
    cmake_print_variables(OptionFile)
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink "${file}" "${CMAKE_BINARY_DIR}/bin/${OptionFile}")
  endforeach()

  enable_testing()
  set(Boost_NO_BOOST_CMAKE ON)
  if ( NOT Boost_FOUND )
    find_package(Boost 1.67.0 COMPONENTS unit_test_framework)
  else()
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
  endif()

else()
  if(AMPGEN_FIND_REQUIRED)
    message(FATAL_ERROR "Unable to find requested AMPGEN installation:${AMPGEN_ERROR_REASON}")
  else()
    if(NOT AMPGEN_FIND_QUIETLY)
      message(STATUS "AMPGEN was not found:${AMPGEN_ERROR_REASON}")
    endif()
  endif()
endif()
