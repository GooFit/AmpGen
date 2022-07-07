# Check if SSE instructions are available on the machine where
# the project is compiled.
include(TestCXXAcceptsFlag)

macro(CHECK_FOR_AVX)
  if(CMAKE_SYSTEM_NAME MATCHES "Linux")
    if(CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "x86_64" OR CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "i.86")
      exec_program(cat ARGS "/proc/cpuinfo" OUTPUT_VARIABLE CPUINFO)

      string(REGEX REPLACE "^.*(avx).*$" "\\1" _SSE_THERE ${CPUINFO})
      string(COMPARE EQUAL "avx" "${_SSE_THERE}" _AVX_TRUE)
      CHECK_CXX_ACCEPTS_FLAG("-mavx" _AVX_OK)

      string(REGEX REPLACE "^.*(avx2).*$" "\\1" _SSE_THERE ${CPUINFO})
      string(COMPARE EQUAL "avx2" "${_SSE_THERE}" _AVX2_TRUE)
      CHECK_CXX_ACCEPTS_FLAG("-mavx2" _AVX2_OK)
    endif()
  elseif(CMAKE_SYSTEM_NAME MATCHES "FreeBSD")
    if(CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "amd64" OR CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "i.86")
      exec_program(cat ARGS "/var/run/dmesg.boot | grep Features" OUTPUT_VARIABLE CPUINFO)

      string(REGEX REPLACE "^.*(AVX).*$" "\\1" _SSE_THERE ${CPUINFO})
      string(COMPARE EQUAL "AVX" "${_SSE_THERE}" _AVX_TRUE)
      CHECK_CXX_ACCEPTS_FLAG("-mavx" _AVX_OK)

      string(REGEX REPLACE "^.*(AVX2).*$" "\\1" _SSE_THERE ${CPUINFO})
      string(COMPARE EQUAL "AVX2" "${_SSE_THERE}" _AVX2_TRUE)
      CHECK_CXX_ACCEPTS_FLAG("-mavx2" _AVX2_OK)
    endif()
  elseif(CMAKE_SYSTEM_NAME MATCHES "Darwin")
    if(NOT CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "arm")
      exec_program("/usr/sbin/sysctl -n machdep.cpu.features machdep.cpu.leaf7_features" OUTPUT_VARIABLE CPUINFO)

      string(REGEX REPLACE "^.*(AVX).*$" "\\1" _SSE_THERE ${CPUINFO})
      string(COMPARE EQUAL "AVX" "${_SSE_THERE}" _AVX_TRUE)
      CHECK_CXX_ACCEPTS_FLAG("-mavx" _AVX_OK)

      string(REGEX REPLACE "^.*(AVX2).*$" "\\1" _SSE_THERE ${CPUINFO})
      string(COMPARE EQUAL "AVX2" "${_SSE_THERE}" _AVX2_TRUE)
      CHECK_CXX_ACCEPTS_FLAG("-mavx2" _AVX2_OK)
    endif()
  elseif(CMAKE_SYSTEM_NAME MATCHES "Windows")
    # TODO
    if(ARCH STREQUAL win32 OR ARCH STREQUAL x64)
      set(_SSE_TRUE true)
      set(_SSE_OK   true)
      set(_SSE2_TRUE true)
      set(_SSE2_OK   true)
    endif()
  endif()

  if ( _AVX2_TRUE )
    set( AVX2_FOUND TRUE )
  endif()
  mark_as_advanced(AVX2_FOUND)

  unset(_AVX_TRUE)
  unset(_AVX_OK)
  unset(_AVX_OK CACHE)
  unset(_AVX2_TRUE)
  unset(_AVX2_OK)
  unset(_AVX2_OK CACHE)

endmacro(CHECK_FOR_AVX)
