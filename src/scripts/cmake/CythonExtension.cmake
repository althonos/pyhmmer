find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
get_property(PYTHON_EXTENSIONS_SOURCE_DIR GLOBAL PROPERTY PYTHON_EXTENSIONS_SOURCE_DIR)

# --- Detect PyInterpreterState_GetID ------------------------------------------

include(CheckSymbolExists)

set(SAFE_CMAKE_REQUIRED_INCLUDES "${CMAKE_REQUIRED_INCLUDES}")
set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${Python_INCLUDE_DIRS})
check_symbol_exists(PyInterpreterState_GetID "stdint.h;stdlib.h;Python.h" HAVE_PYINTERPRETERSTATE_GETID)
set(CMAKE_REQUIRED_INCLUDES "${SAFE_CMAKE_REQUIRED_INCLUDES}")

# --- Detect implementation ----------------------------------------------------

string(TOLOWER "${Python_INTERPRETER_ID}" SYS_IMPLEMENTATION_NAME)

# --- Prepare Cython directives and constants ----------------------------------

set(CYTHON_DIRECTIVES
    -X cdivision=True
    -X nonecheck=False
    -E SSE2_BUILD_SUPPORT=$<IF:$<BOOL:${HAVE_SSE2}>,True,False>
    -E AVX2_BUILD_SUPPORT=$<IF:$<BOOL:${HAVE_AVX2}>,True,False>
    -E NEON_BUILD_SUPPORT=$<IF:$<BOOL:${HAVE_NEON}>,True,False>
    -E MMX_BUILD_SUPPORT=False
    -E AVX512_BUILD_SUPPORT=False
    -E SYS_IMPLEMENTATION_NAME=${SYS_IMPLEMENTATION_NAME}
    -E SYS_VERSION_INFO_MAJOR=${Python_VERSION_MAJOR}
    -E SYS_VERSION_INFO_MINOR=${Python_VERSION_MINOR}
    -E TARGET_CPU=${CMAKE_SYSTEM_PROCESSOR}
    -E TARGET_SYSTEM=${CMAKE_SYSTEM_NAME}
    -E SYS_BYTEORDER=$<IF:$<STREQUAL:${CMAKE_C_BYTE_ORDER},BIG_ENDIAN>,big,little>
    -E PYPY=$<IF:$<STREQUAL:${Python_INTERPRETER_ID},PyPy>,True,False>
    -E PROJECT_VERSION=${CMAKE_PROJECT_VERSION}
)

if(CMAKE_BUILD_TYPE STREQUAL Debug)
  set(CYTHON_DIRECTIVES
    ${CYTHON_DIRECTIVES}
    -X cdivision_warnings=True
    -X warn.undeclared=True
    -X warn.unreachable=True
    -X warn.maybe_uninitialized=True
    -X warn.unused=True
    -X warn.unused_arg=True
    -X warn.unused_result=True
    -X warn.multiple_declarators=True
  )
  if(NOT Python_INTERPRETER_ID STREQUAL PyPy)
    set(CYTHON_DIRECTIVES
      ${CYTHON_DIRECTIVES}
      -X linetrace=true
    )
  endif()
else()
  set(CYTHON_DIRECTIVES
    ${CYTHON_DIRECTIVES}
    -X boundscheck=False
    -X wraparound=False
  )
endif()

# --- Declare Cython extension -------------------------------------------------

macro(cython_extension _name)
  set(multiValueArgs LINKS DIRECTIVES)
  cmake_parse_arguments(CYTHON_EXTENSION "" "" "${multiValueArgs}" ${ARGN} )

  # Make sure that the source directory is known
  if(NOT DEFINED PYTHON_EXTENSIONS_SOURCE_DIR)
    message(FATAL_ERROR "The PYTHON_EXTENSIONS_SOURCE_DIR variable has not been set.")
  endif()

  # Generate C file from Cython file
  add_custom_command(
    OUTPUT 
      ${_name}.c
    COMMENT
      "Cythonizing ${_name}.pyx to ${_name}.c"
    COMMAND
      Python::Interpreter -m cython
        "${CMAKE_CURRENT_SOURCE_DIR}/${_name}.pyx"
        --output-file ${_name}.c
        --depfile
        -I "${CYTHON_HEADERS_DIR}"
        ${CYTHON_DIRECTIVES}
        ${CYTHON_EXTENSION_DIRECTIVES}
    MAIN_DEPENDENCY
      ${_name}.pyx
    DEPFILE
      ${_name}.c.dep
    VERBATIM)

  # Build fully-qualified module name as the target name
  string(REGEX REPLACE "^${PYTHON_EXTENSIONS_SOURCE_DIR}/?" "" _dest_folder ${CMAKE_CURRENT_SOURCE_DIR})
  string(REPLACE "/" "." _target ${_dest_folder}.${_name})

  # Add Python module
  set(EXTENSION_SOURCES ${_name}.pyx ${_name}.c)
  if(EXISTS ${_name}.pxd)
    set(EXTENSION_SOURCES ${EXTENSION_SOURCES} ${_name}.pxd)
  endif()
  python_add_library(${_target} MODULE WITH_SOABI ${EXTENSION_SOURCES})
  set_target_properties(${_target} PROPERTIES OUTPUT_NAME ${_name} )
  target_include_directories(${_target} AFTER PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})  
  target_link_libraries(${_target} PUBLIC ${CYTHON_EXTENSION_LINKS})
  if(HAVE_PYINTERPRETERSTATE_GETID)
    target_compile_definitions(${_target} PUBLIC HAVE_PYINTERPRETERSTATE_GETID)
  endif()

  if(CMAKE_BUILD_TYPE STREQUAL Debug)
    if(NOT Python_INTERPRETER_ID STREQUAL PyPy)
      target_compile_definitions(${_target} PUBLIC CYTHON_TRACE_NOGIL=1)
    endif()
  else()
    target_compile_definitions(${_target} PUBLIC CYTHON_WITHOUT_ASSERTIONS=1)
  endif()

  # Preserve the relative project structure in the install directory
  string(REGEX REPLACE "^${PYTHON_EXTENSIONS_SOURCE_DIR}/?" "" _dest_folder ${CMAKE_CURRENT_SOURCE_DIR})
  install(TARGETS ${_target} DESTINATION ${_dest_folder} )
  message(DEBUG "Install folder for extension ${_name}: ${_dest_folder}")

  # Add the targets to the list of Cython extensions
  get_property(_ext GLOBAL PROPERTY PYRODIGAL_CYTHON_EXTENSIONS)
  list(APPEND _ext ${_target})
  set_property(GLOBAL PROPERTY PYRODIGAL_CYTHON_EXTENSIONS ${_ext})
endmacro()
