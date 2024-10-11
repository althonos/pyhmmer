/* esl_config.h.in  [input to configure]
 * 
 * System-dependent configuration of Easel, by CMake.
 * 
 * This file should be included in all Easel .c files before
 * anything else, because it may set #define's that control
 * behaviour of system includes and system libraries. An example
 * is large file support.
 * 
 */
#ifndef eslCONFIG_INCLUDED
#define eslCONFIG_INCLUDED

/* Version info.
 */
#cmakedefine EASEL_VERSION "@EASEL_VERSION@"
#cmakedefine EASEL_DATE "@EASEL_DATE@"
#cmakedefine EASEL_COPYRIGHT "@EASEL_COPYRIGHT@"
#cmakedefine EASEL_LICENSE "@EASEL_LICENSE@"
#cmakedefine EASEL_URL "@EASEL_URL@"

/* Control of debugging instrumentation */
#cmakedefine eslDEBUGLEVEL   "@eslDEBUGLEVEL@"  // debugging/assertion verbosity level: (0=none;3=most verbose) 
#cmakedefine eslENABLE_ASAN  "@eslENABLE_ASAN@" // some unit tests may need to know if AddressSanitizer is in use
#cmakedefine eslENABLE_TSAN  "@eslENABLE_TSAN@" //   ... ditto, for ThreadSanitizer

/* Optional parallel implementation support */
#cmakedefine eslENABLE_SSE
#cmakedefine eslENABLE_SSE4
#cmakedefine eslENABLE_AVX
#cmakedefine eslENABLE_AVX512
#cmakedefine eslENABLE_NEON
#cmakedefine eslENABLE_VMX

#cmakedefine eslHAVE_NEON_AARCH64

#cmakedefine eslENABLE_CUDA              // Should we build CUDA GPU acceleration?

#cmakedefine HAVE_FLUSH_ZERO_MODE        // On x86 platforms: we can turn off denormalized floating point math,
#cmakedefine HAVE_DENORMALS_ZERO_MODE    //   which often incurs performance penalty. See simdvec.md in HMMER.

#cmakedefine HAVE_MPI
#cmakedefine HAVE_PTHREAD

/* Programs */
#cmakedefine HAVE_GZIP

/* Libraries */
#cmakedefine HAVE_LIBGSL

/* Headers */
#ifndef HAVE_ENDIAN_H
#cmakedefine HAVE_ENDIAN_H @HAVE_ENDIAN_H@
#endif
#ifndef HAVE_INTTYPES_H
#cmakedefine HAVE_INTTYPES_H @HAVE_INTTYPES_H@
#endif
#ifndef HAVE_STDINT_H
#cmakedefine HAVE_STDINT_H @HAVE_STDINT_H@
#endif
#ifndef HAVE_UNISTD_H
#cmakedefine HAVE_UNISTD_H @HAVE_UNISTD_H@
#endif
#ifndef HAVE_SYS_TYPES_H
#cmakedefine HAVE_SYS_TYPES_H @HAVE_SYS_TYPES_H@
#endif
#ifndef HAVE_STRINGS_H
#cmakedefine HAVE_STRINGS_H @HAVE_STRINGS_H@
#endif
#ifndef HAVE_NETINET_IN_H
#cmakedefine HAVE_NETINET_IN_H @HAVE_NETINET_IN_H@	/* On FreeBSD, you need netinet/in.h for struct sockaddr_in */
#endif
#ifndef HAVE_SYS_PARAM_H
#cmakedefine HAVE_SYS_PARAM_H @HAVE_SYS_PARAM_H@
#endif
#ifndef HAVE_SYS_SYSCTL_H
#cmakedefine HAVE_SYS_SYSCTL_H @HAVE_SYS_SYSCTL_H@
#endif

/* Types */
#cmakedefine WORDS_BIGENDIAN
// #cmakedefine int8_t
// #cmakedefine int16_t
// #cmakedefine int32_t
// #cmakedefine int64_t
// #cmakedefine uint8_t
// #cmakedefine uint16_t
// #cmakedefine uint32_t
// #cmakedefine uint64_t
// #cmakedefine off_t

/* Compiler characteristics */
#cmakedefine HAVE_FUNC_ATTRIBUTE_NORETURN // Compiler supports __attribute__((__noreturn__)), helps w/ clang static analysis.
#cmakedefine HAVE_FUNC_ATTRIBUTE_FORMAT   // Compiler supports __attribute__((format(a,b,c))), typechecking printf-like functions

/* Functions */
#ifndef HAVE_ALIGNED_ALLOC
#cmakedefine HAVE_ALIGNED_ALLOC   @HAVE_ALIGNED_ALLOC@  // esl_alloc
#endif
#ifndef HAVE_ERFC
#cmakedefine HAVE_ERFC            @HAVE_ERFC@           // esl_stats
#endif
#ifndef HAVE_GETCWD
#cmakedefine HAVE_GETCWD          @HAVE_GETCWD@         // esl_getcwd
#endif
#ifndef HAVE_GETPID
#cmakedefine HAVE_GETPID          @HAVE_GETPID@         // esl_random
#endif
#ifndef HAVE__MM_MALLOC
#cmakedefine HAVE__MM_MALLOC      @HAVE__MM_MALLOC@     // esl_alloc
#endif
#ifndef HAVE_POPEN
#cmakedefine HAVE_POPEN           @HAVE_POPEN@          // various file parsers that check for piped input
#endif
#ifndef HAVE_POSIX_MEMALIGN
#cmakedefine HAVE_POSIX_MEMALIGN  @HAVE_POSIX_MEMALIGN@ // esl_alloc
#endif
#ifndef HAVE_STRCASECMP
#cmakedefine HAVE_STRCASECMP      @HAVE_STRCASECMP@     // easel::esl_strcasecmp()
#endif
#ifndef HAVE_STRSEP
#cmakedefine HAVE_STRSEP          @HAVE_STRSEP@         // easel::esl_strsep()
#endif
#ifndef HAVE_SYSCONF
#cmakedefine HAVE_SYSCONF         @HAVE_SYSCONF@        // esl_threads, asking system for cpu number
#endif
#ifndef HAVE_SYSCTL
#cmakedefine HAVE_SYSCTL          @HAVE_SYSCTL@         // esl_threads, ""
#endif
#ifndef HAVE_TIMES
#cmakedefine HAVE_TIMES           @HAVE_TIMES@          // esl_stopwatch
#endif
#ifndef HAVE_FSEEKO
#cmakedefine HAVE_FSEEKO
#endif

/* System services */
#cmakedefine _FILE_OFFSET_BITS    // Large file support; possibly archaic now?
#cmakedefine _LARGE_FILES         //  ""
#cmakedefine _LARGEFILE_SOURCE    //  ""

 
/* Function behavior */
#cmakedefine eslSTOPWATCH_HIGHRES

#endif /*eslCONFIG_INCLUDED*/

