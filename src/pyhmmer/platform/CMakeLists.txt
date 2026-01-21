include(CythonExtension)

if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    cython_extension(win32  LINKS libeasel)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    cython_extension(linux  LINKS libeasel)
elseif((CMAKE_SYSTEM_NAME STREQUAL "Darwin") OR (CMAKE_SYSTEM_NAME MATCHES ".*BSD"))
    cython_extension(bsd    LINKS libeasel)
endif()

if(DEFINED PYHMMER_INSTALL_LIBS_DIR)
    cmake_path(APPEND SKBUILD_PLATLIB_DIR "${PYHMMER_INSTALL_LIBS_DIR}" "include" "pyhmmer" "platform" OUTPUT_VARIABLE INCLUDE_DEST_FOLDER)
    install(FILES "linux.h" "bsd.h" "util.h" DESTINATION "${INCLUDE_DEST_FOLDER}")
endif()
