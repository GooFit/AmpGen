find_package(Boost 1.67.0 COMPONENTS unit_test_framework)

if ( Boost_FOUND )
  include_directories (${Boost_INCLUDE_DIRS})
  file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
  foreach(testSrc ${TEST_SRCS})
    get_filename_component(testName ${testSrc} NAME_WE)
    add_executable(${testName} ${testSrc})
    target_compile_options(${testName} PRIVATE -Wno-deprecated-builtins -Wno-deprecated-declarations -Wno-unused-but-set-variable )
    set_target_properties(${testName} 
      PROPERTIES 
      RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_TEST_OUTPUT_DIRECTORY}" 
      RUNTIME_OUTPUT_DIRECTORY_DEBUG   "${CMAKE_TEST_OUTPUT_DIRECTORY}" 
      RUNTIME_OUTPUT_DIRECTORY         "${CMAKE_TEST_OUTPUT_DIRECTORY}" 
      EXECUTABLE_OUTPUT_DIRECTORY      "${CMAKE_TEST_OUTPUT_DIRECTORY}" 
    )
    target_link_libraries(${testName} PRIVATE ${Boost_LIBRARIES} libAmpGen)
    if( UNIX AND NOT APPLE )
      target_link_libraries( ${testName} PRIVATE rt)
    endif()
    add_test(NAME ${testName} WORKING_DIRECTORY ${CMAKE_TEST_OUTPUT_DIRECTORY} COMMAND ${CMAKE_TEST_OUTPUT_DIRECTORY}/${testName} ) 
  endforeach(testSrc)
else()
  message( WARNING "Warning: Boost (version >= 1.67.0) required to build unit tests\n")
endif()
