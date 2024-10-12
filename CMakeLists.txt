cmake_minimum_required(VERSION 3.20)
project(${SKBUILD_PROJECT_NAME} VERSION ${SKBUILD_PROJECT_VERSION} LANGUAGES C)

set(BUILD_SHARED_LIBS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set_property(GLOBAL PROPERTY PYTHON_EXTENSIONS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

# --- Detect SIMD --------------------------------------------------------------

include("src/scripts/cmake/FindNEON.cmake")
include("src/scripts/cmake/FindSSE2.cmake")
include("src/scripts/cmake/FindSSE4.cmake")
include("src/scripts/cmake/FindVMX.cmake")
include("src/scripts/cmake/Destatic.cmake")
include("src/scripts/cmake/CythonExtension.cmake")

if(HAVE_SSE4)
    set(eslENABLE_SSE4 true)
    string(STRIP "${SSE4_C_FLAGS}" IMPL_FLAGS)
    if(IMPL_FLAGS)
        add_compile_options(${IMPL_FLAGS})
    endif()
endif()

if(HAVE_SSE2)
    set(eslENABLE_SSE true)
    set(HMMER_IMPL "SSE")
    string(STRIP "${SSE2_C_FLAGS}" IMPL_FLAGS)
elseif(HAVE_NEON)
    set(eslENABLE_NEON true)
    set(HMMER_IMPL "NEON")
    string(STRIP "${NEON_C_FLAGS}" IMPL_FLAGS)
elseif(HAVE_VMX)
    set(eslENABLE_VMX true)
    set(HMMER_IMPL "VMX")
    string(STRIP "${VMX_C_FLAGS}" IMPL_FLAGS)
else()
    message(FATAL_ERROR "Unsupported platform, none of SSE2, NEON or AltiVec available.")
endif()

if(IMPL_FLAGS)
    add_compile_options(${IMPL_FLAGS})
endif()

# --- Setup include folders ----------------------------------------------------

set(CYTHON_HEADERS_DIR ${CMAKE_CURRENT_LIST_DIR}/include)

# --- Compile source code ------------------------------------------------------

add_subdirectory("src")
