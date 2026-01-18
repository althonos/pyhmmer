include(CythonExtension)

add_subdirectory(fileobj)

cython_extension(errors   LINKS libeasel)
cython_extension(easel    LINKS libeasel fileobj)
cython_extension(plan7    LINKS libeasel libhmmer fileobj DIRECTIVES -E HMMER_IMPL=${HMMER_IMPL})
cython_extension(daemon   LINKS libeasel libhmmer fileobj DIRECTIVES -E HMMER_IMPL=${HMMER_IMPL})
