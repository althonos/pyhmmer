
file(COPY ../vendor/easel DESTINATION .)
file(COPY ./easel DESTINATION .)
destatic(../vendor/easel/esl_sqio_ascii.c ${CMAKE_CURRENT_BINARY_DIR}/easel/esl_sqio_ascii.full.c)

file(COPY ../vendor/hmmer DESTINATION .)
file(COPY ./hmmer DESTINATION .)
destatic(../vendor/hmmer/src/p7_hmmfile.c ${CMAKE_CURRENT_BINARY_DIR}/hmmer/src/p7_hmmfile.full.c)

add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/easel ${CMAKE_CURRENT_BINARY_DIR}/easel/build)
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/hmmer ${CMAKE_CURRENT_BINARY_DIR}/hmmer/build)
add_subdirectory(pyhmmer)
