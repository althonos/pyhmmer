include(CythonExtension)

# cython_extension(errors LINKS libeasel)
# cython_extension(easel  LINKS libeasel)
# cython_extension(plan7  LINKS libeasel libhmmer DIRECTIVES -E HMMER_IMPL=${HMMER_IMPL})
# cython_extension(daemon LINKS libeasel libhmmer DIRECTIVES -E HMMER_IMPL=${HMMER_IMPL})

if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    cython_extension(win32 LINKS libeasel)
endif()
