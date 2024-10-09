#.rst:
# FindVMX
# --------
#
# Finds VMX support
#
# This module can be used to detect VMX support in a C compiler.  If
# the compiler supports VMX, the flags required to compile with
# VMX support are returned in variables for the different languages.
# The variables may be empty if the compiler does not need a special
# flag to support VMX.
#
# The following variables are set:
#
# ::
#
#    VMX_C_FLAGS - flags to add to the C compiler for VMX support
#    VMX_FOUND - true if VMX is detected
#
#=============================================================================

set(_VMX_REQUIRED_VARS)
set(CMAKE_REQUIRED_QUIET_SAVE ${CMAKE_REQUIRED_QUIET})
set(CMAKE_REQUIRED_QUIET ${VMX_FIND_QUIETLY})

# sample VMX source code to test
set(VMX_C_TEST_SOURCE
"
#include <altivec.h>
int foo() {
    vector float a = vec_splats(1.0);
    float f;
    vec_ste(a, 0, &f);
    return (f == 1) ? 0 : 1;
}
int main(void) { return foo(); }
")

# if these are set then do not try to find them again,
# by avoiding any try_compiles for the flags
if((DEFINED VMX_C_FLAGS) OR (DEFINED HAVE_VMX))
else()
  if(WIN32)
    set(VMX_C_FLAG_CANDIDATES
      #Empty, if compiler automatically accepts VMX
      " "
      "/arch:VMX")
  else()
    set(VMX_C_FLAG_CANDIDATES
      #Empty, if compiler automatically accepts VMX
      " "
    )
  endif()

  include(CheckCSourceCompiles)

  foreach(FLAG IN LISTS VMX_C_FLAG_CANDIDATES)
    set(SAFE_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
    set(CMAKE_REQUIRED_FLAGS "${FLAG}")
    unset(HAVE_VMX CACHE)
    if(NOT CMAKE_REQUIRED_QUIET)
      message(STATUS "Try VMX C flag = [${FLAG}]")
    endif()
    check_c_source_compiles("${VMX_C_TEST_SOURCE}" HAVE_VMX)
    set(CMAKE_REQUIRED_FLAGS "${SAFE_CMAKE_REQUIRED_FLAGS}")
    if(HAVE_VMX)
      set(VMX_C_FLAGS_INTERNAL "${FLAG}")
      break()
    endif()
  endforeach()

  unset(VMX_C_FLAG_CANDIDATES)

  set(VMX_C_FLAGS "${VMX_C_FLAGS_INTERNAL}"
    CACHE STRING "C compiler flags for VMX intrinsics")
endif()

list(APPEND _VMX_REQUIRED_VARS VMX_C_FLAGS)

set(CMAKE_REQUIRED_QUIET ${CMAKE_REQUIRED_QUIET_SAVE})

if(_VMX_REQUIRED_VARS)
  include(FindPackageHandleStandardArgs)

  find_package_handle_standard_args(VMX
                                    REQUIRED_VARS ${_VMX_REQUIRED_VARS})

  mark_as_advanced(${_VMX_REQUIRED_VARS})

  unset(_VMX_REQUIRED_VARS)
else()
  message(SEND_ERROR "FindVMX requires C or CXX language to be enabled")
endif()
