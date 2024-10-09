cmake_minimum_required(VERSION 3.20)
project(${SKBUILD_PROJECT_NAME} VERSION ${SKBUILD_PROJECT_VERSION} LANGUAGES C)

set(BUILD_SHARED_LIBS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set_property(GLOBAL PROPERTY PYTHON_EXTENSIONS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

# --- Detect SIMD --------------------------------------------------------------

include("src/scripts/cmake/FindNEON.cmake")
include("src/scripts/cmake/FindSSE2.cmake")
include("src/scripts/cmake/FindSSE4.cmake")
# include("src/scripts/cmake/FindVMX.cmake")
include("src/scripts/cmake/Destatic.cmake")
include("src/scripts/cmake/CythonExtension.cmake")

if(HAVE_SSE4)
    set(eslENABLE_SSE4 true)
    if(SSE4_C_FLAGS)
        add_compile_options(${SSE4_C_FLAGS})
    endif()
endif()

if(HAVE_SSE2)
    set(eslENABLE_SSE true)
    set(HMMER_IMPL "SSE")
    if(SSE2_C_FLAGS)
        add_compile_options(${SSE2_C_FLAGS})
    endif()
elseif(HAVE_NEON)
    set(eslENABLE_NEON true)
    set(HMMER_IMPL "NEON")
    if(NEON_C_FLAGS)
        add_compile_options(${NEON_C_FLAGS})
    endif()
else()
    error("Unsupported platform, neither SSE2 nor NEON available.")
endif()

# --- Setup include folders ----------------------------------------------------

# include_directories(BEFORE PUBLIC ${CMAKE_CURRENT_LIST_DIR}/src/pyhmmer)
# include_directories(AFTER  PUBLIC vendor/easel/include)
# include_directories(AFTER  PUBLIC vendor/hmmer)

set(CYTHON_HEADERS_DIR ${CMAKE_CURRENT_LIST_DIR}/include)

# --- Compile Python extension -------------------------------------------------

add_subdirectory("src")
