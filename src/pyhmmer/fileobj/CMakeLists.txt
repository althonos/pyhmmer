include(CythonExtension)

if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    cython_extension(win32  LINKS libeasel)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    cython_extension(cookie LINKS libeasel)
elseif((CMAKE_SYSTEM_NAME STREQUAL "Darwin") OR (CMAKE_SYSTEM_NAME MATCHES "*BSD"))
    cython_extension(fun    LINKS libeasel)
endif()

add_library(fileobj INTERFACE)
target_include_directories(fileobj INTERFACE ${CMAKE_CURRENT_LIST_DIR})