include(CheckIncludeFile)
include(CheckFunctionExists)

CHECK_INCLUDE_FILE("stdio.h" HAVE_STDIO_H)
CHECK_INCLUDE_FILE("stdlib.h" HAVE_STLIB_H)
CHECK_INCLUDE_FILE("string.h" HAVE_STRING_H)
CHECK_INCLUDE_FILE("inttypes.h" HAVE_INTTYPES_H)
CHECK_INCLUDE_FILE("stdint.h" HAVE_STDINT_H)
CHECK_INCLUDE_FILE("strings.h" HAVE_STRINGS_H)
CHECK_INCLUDE_FILE("sys/stat.h" HAVE_SYS_STAT_H)
CHECK_INCLUDE_FILE("sys/types.h" HAVE_SYS_TYPES_H)
CHECK_INCLUDE_FILE("unistd.h" HAVE_UNISTD_H)
CHECK_INCLUDE_FILE("endian.h" HAVE_ENDIAN_H)
CHECK_INCLUDE_FILE("netinet/in.h" HAVE_NETINET_IN_H)
CHECK_INCLUDE_FILE("sys/param.h" HAVE_SYS_PARAM_H)
CHECK_INCLUDE_FILE("sys/sysctl.h" HAVE_SYS_SYSCTL_H)

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

set(EASEL_SOURCES
    ${CMAKE_CURRENT_BINARY_DIR}/esl_config.h
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
    esl_dsqdata.c
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
    esl_sqio_ascii.full.c
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

add_library(libeasel ${EASEL_SOURCES})
target_include_directories(libeasel PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(libeasel PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
