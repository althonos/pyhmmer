include(CythonExtension)

add_subdirectory(platform)

cython_extension(errors   LINKS libeasel)
cython_extension(easel    LINKS libeasel platform)
cython_extension(plan7    LINKS libeasel platform libhmmer DIRECTIVES -E HMMER_IMPL=${HMMER_IMPL})
cython_extension(daemon   LINKS libeasel platform libhmmer DIRECTIVES -E HMMER_IMPL=${HMMER_IMPL})
