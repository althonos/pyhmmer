# --- Handle config.h generation *a la* ./configure ----------------------------

include(CheckIncludeFile)
include(CheckFunctionExists)

CHECK_INCLUDE_FILE("stdio.h" HAVE_STDIO_H)
CHECK_INCLUDE_FILE("stdlib.h" HAVE_STLIB_H)
CHECK_INCLUDE_FILE("string.h" HAVE_STRING_H)
CHECK_INCLUDE_FILE("inttypes.h" HAVE_INTTYPES_H)
CHECK_INCLUDE_FILE("stdint.h" HAVE_STDINT_H)
CHECK_INCLUDE_FILE("strings.h" HAVE_STRINGS_H)
CHECK_INCLUDE_FILE("unistd.h" HAVE_UNISTD_H)
CHECK_INCLUDE_FILE("endian.h" HAVE_ENDIAN_H)
CHECK_INCLUDE_FILE("netinet/in.h" HAVE_NETINET_IN_H)
CHECK_INCLUDE_FILE("sys/mman.h" HAVE_SYS_MMAN_H)
CHECK_INCLUDE_FILE("sys/param.h" HAVE_SYS_PARAM_H)
CHECK_INCLUDE_FILE("sys/stat.h" HAVE_SYS_STAT_H)
CHECK_INCLUDE_FILE("sys/sysctl.h" HAVE_SYS_SYSCTL_H)
CHECK_INCLUDE_FILE("sys/types.h" HAVE_SYS_TYPES_H)
CHECK_INCLUDE_FILE("syslog.h" HAVE_SYSLOG_H)

CHECK_FUNCTION_EXISTS(aligned_alloc HAVE_ALIGNED_ALLOC)
CHECK_FUNCTION_EXISTS(erfc HAVE_ERFC)
CHECK_FUNCTION_EXISTS(getpid HAVE_GETPID)
CHECK_FUNCTION_EXISTS(_mm_malloc HAVE__MM_MALLOC)
CHECK_FUNCTION_EXISTS(popen HAVE_POPEN)
CHECK_FUNCTION_EXISTS(posix_memalign HAVE_POSIX_MEMALIGN)
CHECK_FUNCTION_EXISTS(strcasecmp HAVE_STRCASECMP)
CHECK_FUNCTION_EXISTS(strsep HAVE_STRSEP)
CHECK_FUNCTION_EXISTS(sysconf HAVE_SYSCONF)
CHECK_FUNCTION_EXISTS(sysctl HAVE_SYSCTL)
CHECK_FUNCTION_EXISTS(times HAVE_TIMES)
CHECK_FUNCTION_EXISTS(fseeko HAVE_FSEEKO)

set(EASEL_DATE "Aug 2023")
set(EASEL_COPYRIGHT "Copyright (C) 2023 Howard Hughes Medical Institute.")
set(EASEL_LICENSE  "Freely distributed under the BSD open source license.")
set(EASEL_VERSION "0.49")
set(EASEL_URL "http://bioeasel.org/")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/esl_config.h.cmake esl_config.h)
file(READ ${CMAKE_CURRENT_BINARY_DIR}/esl_config.h ESL_CONFIG_H)

message(DEBUG "-- Show configuration ")
message(DEBUG "HAVE_STLIB_H: ${HAVE_STLIB_H}")
message(DEBUG "HAVE_STRINGS_H: ${HAVE_STRINGS_H}")
message(DEBUG "HAVE_UNISTD_H: ${HAVE_UNISTD_H}")
message(DEBUG "HAVE_STDINT_H: ${HAVE_STDINT_H}")
message(DEBUG "-- Dumping esl_config.h --")
message(DEBUG "${ESL_CONFIG_H}")
message(DEBUG "-- End esl_config.h --")

# --- Manually declare required sources only -----------------------------------

set(EASEL_SOURCE_DIR
    "${PROJECT_SOURCE_DIR}/vendor/easel"
)

set(EASEL_SOURCES
    easel.c
    esl_alloc.c
    esl_alphabet.c
    esl_arr2.c
    esl_arr3.c
    esl_avx512.c
    esl_avx.c
    esl_bitfield.c
    esl_buffer.c
    esl_cluster.c
    esl_composition.c
    esl_cpu.c
    esl_dirichlet.c
    esl_distance.c
    esl_dmatrix.c
    # esl_dsqdata.c  # <-- don't compile esl_dsqdata.c because it is not used, and depends heavily on pthreads
    esl_exponential.c
    esl_fileparser.c
    esl_gamma.c
    esl_gencode.c
    esl_getopts.c
    esl_gev.c
    esl_graph.c
    esl_gumbel.c
    esl_heap.c
    esl_histogram.c
    esl_hmm.c
    esl_huffman.c
    esl_hyperexp.c
    esl_iset.c
    esl_json.c
    esl_keyhash.c
    esl_lognormal.c
    esl_matrixops.c
    esl_mem.c
    esl_minimizer.c
    esl_mixdchlet.c
    esl_mixgev.c
    esl_mpi.c
    esl_msa.c
    esl_msacluster.c
    esl_msafile2.c
    esl_msafile_a2m.c
    esl_msafile_afa.c
    esl_msafile.c
    esl_msafile_clustal.c
    esl_msafile_phylip.c
    esl_msafile_psiblast.c
    esl_msafile_selex.c
    esl_msafile_stockholm.c
    esl_msashuffle.c
    esl_msaweight.c
    esl_neon.c
    esl_normal.c
    esl_paml.c
    esl_quicksort.c
    esl_rand64.c
    esl_random.c
    esl_randomseq.c
    esl_ratematrix.c
    esl_recorder.c
    esl_red_black.c
    esl_regexp.c
    esl_rootfinder.c
    esl_scorematrix.c
    esl_sq.c
    esl_sqio_ascii.c
    esl_sqio.c
    esl_sqio_ncbi.c
    esl_sse.c
    esl_ssi.c
    esl_stack.c
    esl_stats.c
    esl_stopwatch.c
    esl_stretchexp.c
    esl_subcmd.c
    esl_swat.c
    esl_threads.c
    esl_tree.c
    esl_varint.c
    esl_vectorops.c
    esl_vmx.c
    esl_weibull.c
    esl_workqueue.c
    esl_wuss.c
    interface_gsl.c
    interface_lapack.c
)

set(EASEL_HEADERS
    easel.h
    esl_alloc.h
    esl_alphabet.h
    esl_arr2.h
    esl_arr3.h
    esl_avx512.h
    esl_avx.h
    esl_bitfield.h
    esl_buffer.h
    esl_cluster.h
    esl_composition.h
    esl_cpu.h
    esl_dirichlet.h
    esl_distance.h
    esl_dmatrix.h
    esl_dsqdata.h
    esl_exponential.h
    esl_fileparser.h
    esl_gamma.h
    esl_gencode.h
    esl_getopts.h
    esl_gev.h
    esl_graph.h
    esl_gumbel.h
    esl_heap.h
    esl_histogram.h
    esl_hmm.h
    esl_huffman.h
    esl_hyperexp.h
    esl_iset.h
    esl_json.h
    esl_keyhash.h
    esl_lognormal.h
    esl_matrixops.h
    esl_mem.h
    esl_minimizer.h
    esl_mixdchlet.h
    esl_mixgev.h
    esl_mpi.h
    esl_msacluster.h
    esl_msafile2.h
    esl_msafile_a2m.h
    esl_msafile_afa.h
    esl_msafile_clustal.h
    esl_msafile.h
    esl_msafile_phylip.h
    esl_msafile_psiblast.h
    esl_msafile_selex.h
    esl_msafile_stockholm.h
    esl_msa.h
    esl_msashuffle.h
    esl_msaweight.h
    esl_neon.h
    esl_normal.h
    esl_paml.h
    esl_quicksort.h
    esl_rand64.h
    esl_random.h
    esl_randomseq.h
    esl_ratematrix.h
    esl_recorder.h
    esl_red_black.h
    esl_regexp.h
    esl_rootfinder.h
    esl_scorematrix.h
    esl_sq.h
    esl_sqio_ascii.h
    esl_sqio.h
    esl_sqio_ncbi.h
    esl_sse.h
    esl_ssi.h
    esl_stack.h
    esl_stats.h
    esl_stopwatch.h
    esl_stretchexp.h
    esl_subcmd.h
    esl_threads.h
    esl_tree.h
    esl_varint.h
    esl_vectorops.h
    esl_vmx.h
    esl_weibull.h
    esl_workqueue.h
    esl_wuss.h
    interface_gsl.h
    interface_lapack.h
)

set(EASEL_PATCHED_SOURCES
    ${CMAKE_CURRENT_BINARY_DIR}/esl_config.h
)

# --- Copy/patch source files if needed ----------------------------------------

foreach(_file IN ITEMS ${EASEL_SOURCES} ${EASEL_HEADERS})
    if(EXISTS ${PROJECT_SOURCE_DIR}/patches/${_file}.patch)
        if(NOT EXISTS ${EASEL_SOURCE_DIR}/${_file})
            add_custom_command(
                OUTPUT
                    ${_file}
                COMMENT
                    "Creating ${_file}"
                COMMAND
                    cmake -E touch ${CMAKE_CURRENT_BINARY_DIR}/empty.txt
                COMMAND
                    ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/src/scripts/apply_patch.py
                        --input ${CMAKE_CURRENT_BINARY_DIR}/empty.txt
                        --patch ${PROJECT_SOURCE_DIR}/patches/${_file}.patch
                        --output ${CMAKE_CURRENT_BINARY_DIR}/${_file}
                DEPENDS
                    ${PROJECT_SOURCE_DIR}/patches/${_file}.patch
            )
        else()
            add_custom_command(
                OUTPUT
                    ${_file}
                COMMENT
                    "Patching ${_file}"
                COMMAND
                    ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/src/scripts/apply_patch.py
                        --input ${EASEL_SOURCE_DIR}/${_file}
                        --patch ${PROJECT_SOURCE_DIR}/patches/${_file}.patch
                        --output ${CMAKE_CURRENT_BINARY_DIR}/${_file}
                DEPENDS
                    ${EASEL_SOURCE_DIR}/${_file}
                    ${PROJECT_SOURCE_DIR}/patches/${_file}.patch
            )
        endif()
    else()
        add_custom_command(
            OUTPUT
                ${_file}
            COMMENT
                "Copying ${_file}"
            COMMAND
                cmake -E copy ${EASEL_SOURCE_DIR}/${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file}
            DEPENDS
                ${EASEL_SOURCE_DIR}/${_file}
        )
    endif()
    set(EASEL_PATCHED_SOURCES ${EASEL_PATCHED_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/${_file})
endforeach()

# --- Build libeasel with copied/patched sources -------------------------------

add_library(libeasel ${EASEL_PATCHED_SOURCES})
target_include_directories(libeasel PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(libeasel PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# On Windows, compile with 64-bit file offsets to allow reading 
# SSI files created on 64-bit platforms.
if(WIN32)
    target_compile_definitions(libeasel PRIVATE -D_LARGE_FILES=1 -D_FILE_OFFSET_BITS=64)
endif()

# --- Install shared libraries & headers if required ---------------------------

if(DEFINED PYHMMER_INSTALL_LIBS_DIR)
    install(TARGETS libeasel DESTINATION "${SKBUILD_PLATLIB_DIR}/${PYHMMER_INSTALL_LIBS_DIR}")
    if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
        set_target_properties(libeasel PROPERTIES INSTALL_RPATH "@loader_path/")
    else()
        set_target_properties(libeasel PROPERTIES INSTALL_RPATH "\$ORIGIN/")
    endif()
    cmake_path(APPEND SKBUILD_PLATLIB_DIR "${PYHMMER_INSTALL_LIBS_DIR}" "include" "libeasel" OUTPUT_VARIABLE INCLUDE_DEST_FOLDER)
    foreach(_file IN ITEMS ${EASEL_PATCHED_SOURCES})
        cmake_path(GET _file EXTENSION LAST_ONLY _ext)
        if(_ext STREQUAL ".h")
            cmake_path(RELATIVE_PATH _file BASE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" OUTPUT_VARIABLE _local)
            cmake_path(GET _local PARENT_PATH _dest)
            install(FILES "${_file}" DESTINATION "${INCLUDE_DEST_FOLDER}/${_dest}")
        endif()
    endforeach()
  endif()
