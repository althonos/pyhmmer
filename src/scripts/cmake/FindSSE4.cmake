#.rst:
# FindSSE4
# --------
#
# Finds SSE4 support
#
# This module can be used to detect SSE4 support in a C compiler.  If
# the compiler supports SSE4, the flags required to compile with
# SSE4 support are returned in variables for the different languages.
# The variables may be empty if the compiler does not need a special
# flag to support SSE4.
#
# The following variables are set:
#
# ::
#
#    SSE4_C_FLAGS - flags to add to the C compiler for SSE4 support
#    SSE4_FOUND - true if SSE4 is detected
#
#=============================================================================

set(_SSE4_REQUIRED_VARS)
set(CMAKE_REQUIRED_QUIET_SAVE ${CMAKE_REQUIRED_QUIET})
set(CMAKE_REQUIRED_QUIET ${SSE4_FIND_QUIETLY})

# sample SSE4 source code to test
set(SSE4_C_TEST_SOURCE
"
#if defined(_MSC_VER)
#include <intrin.h>
#else
#include <smmintrin.h>
#endif
int foo() {
    __m128  a = _mm_set1_ps(1.0);
            a = _mm_blend_ps(a, _mm_setzero_ps(), 0b1111);
    float   x = _mm_extract_ps(a, 1);
    return (x == 1) ? 0 : 1;
}
int main(void) { return foo(); }
")

# if these are set then do not try to find them again,
# by avoiding any try_compiles for the flags
if((DEFINED SSE4_C_FLAGS) OR (DEFINED HAVE_SSE4))
else()
  if(WIN32)
    set(SSE4_C_FLAG_CANDIDATES
      #Empty, if compiler automatically accepts SSE4
      " "
      "/arch:SSE4")
  else()
    set(SSE4_C_FLAG_CANDIDATES
      #Empty, if compiler automatically accepts SSE4
      " "
      #clang
      "-msse4.1"
    )
  endif()

  include(CheckCSourceCompiles)

  foreach(FLAG IN LISTS SSE4_C_FLAG_CANDIDATES)
    set(SAFE_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
    set(CMAKE_REQUIRED_FLAGS "${FLAG}")
    unset(HAVE_SSE4 CACHE)
    if(NOT CMAKE_REQUIRED_QUIET)
      message(STATUS "Try SSE4 C flag = [${FLAG}]")
    endif()
    check_c_source_compiles("${SSE4_C_TEST_SOURCE}" HAVE_SSE4)
    set(CMAKE_REQUIRED_FLAGS "${SAFE_CMAKE_REQUIRED_FLAGS}")
    if(HAVE_SSE4)
      set(SSE4_C_FLAGS_INTERNAL "${FLAG}")
      break()
    endif()
  endforeach()

  unset(SSE4_C_FLAG_CANDIDATES)

  set(SSE4_C_FLAGS "${SSE4_C_FLAGS_INTERNAL}"
    CACHE STRING "C compiler flags for SSE4 intrinsics")
endif()

list(APPEND _SSE4_REQUIRED_VARS SSE4_C_FLAGS)

set(CMAKE_REQUIRED_QUIET ${CMAKE_REQUIRED_QUIET_SAVE})

if(_SSE4_REQUIRED_VARS)
  include(FindPackageHandleStandardArgs)

  find_package_handle_standard_args(SSE4
                                    REQUIRED_VARS ${_SSE4_REQUIRED_VARS})

  mark_as_advanced(${_SSE4_REQUIRED_VARS})

  unset(_SSE4_REQUIRED_VARS)
else()
  message(SEND_ERROR "FindSSE4 requires C or CXX language to be enabled")
endif()
