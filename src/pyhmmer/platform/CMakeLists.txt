include(CythonExtension)

add_library(platform INTERFACE)
cmake_path(GET CMAKE_CURRENT_LIST_DIR PARENT_PATH _pyhmmer_folder)
target_include_directories(platform INTERFACE "${_pyhmmer_folder}")
message(STATUS ${_pyhmmer_folder})

if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    cython_extension(win32  LINKS libeasel platform)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    cython_extension(linux  LINKS libeasel platform)
elseif((CMAKE_SYSTEM_NAME STREQUAL "Darwin") OR (CMAKE_SYSTEM_NAME MATCHES ".*BSD"))
    cython_extension(bsd    LINKS libeasel platform)
endif()

