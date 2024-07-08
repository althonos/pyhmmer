# --- Generate configuration ---

easel_conf = configuration_data({
  'EASEL_DATE'      : '"Aug 2023"',
  'EASEL_COPYRIGHT' : '"Copyright (C) 2023 Howard Hughes Medical Institute"',
  'EASEL_LICENSE'   : '"Freely distributed under the BSD open source license."',
  'EASEL_VERSION'   : '"0.49"',
  'EASEL_URL'       : '"http://bioeasel.org/"',
})

_easel_includes = [
  'endian.h',
  'inttypes.h',
  'stdint.h',
  'unistd.h',
  'sys/types.h',
  'strings.h',
  'netinet/in.h',
  'sys/param.h',
  'sys/sysctl.h',
]
foreach h : _easel_includes
  easel_conf.set('HAVE_' + h.underscorify().to_upper(), cc.has_header(h))
endforeach

_easel_functions = [
  'aligned_alloc',
  'erfc',
  'getpid',
  '_mm_malloc',
  'popen',
  'posix_memalign',
  'strcasecmp',
  'strsep',
  'sysconf',
  'sysctl',
  'times',
  'fseeko'
] 
foreach f : _easel_functions
  easel_conf.set('HAVE_' + h.underscorify().to_upper(), cc.has_function(f))
endforeach

_easel_attributes = [
  'noreturn',
  'format'
]
foreach a : _easel_attributes
  easel_conf.set('HAVE_FUNC_ATTRIBUTE' + a.underscorify().to_upper(), cc.has_function_attribute(a))
endforeach

if target_machine.cpu_family() == 'x86_64'
  easel_conf.set('eslENABLE_SSE', true)
endif

if target_machine.cpu_family() == 'arm'
  easel_conf.set('eslENABLE_NEON', true)
endif

if target_machine.cpu_family() == 'aarch64'
  easel_conf.set('eslENABLE_NEON', true)
  easel_conf.set('eslHAVE_NEON_AARCH64', true)
endif

if target_machine.endian() == 'big'
  easel_conf.set('WORDS_BIGENDIAN', true)
endif

if mpi.found()
  easel_conf.set('HAVE_MPI', true)
endif

_ftz_code = (
  '''
  #include <xmmintrin.h>,
  int main(int argc, char** argv) {
    _MM_SET_FLUSH_ZERO_MODE (_MM_FLUSH_ZERO_ON);
  }
  '''
)
if cc.compiles(_ftz_code, name : 'flush-to-zero')
  easel_conf.set('HAVE_FLUSH_ZERO_MODE', true)
endif

_daz_code = (
  '''
  #include <pmmintrin.h>,
  int main(int argc, char** argv) {
    _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
  }
  '''
)
if cc.compiles(_daz_code, name : 'denormals-are-zero')
  easel_conf.set('HAVE_DENORMALS_ZERO_MODE', true)
endif

# --- Write header ---

easel_config_h = configure_file(
  output : 'esl_config.h', 
  configuration : easel_conf
)

# --- Compile library ---

easel_sources = [
  'easel.c',
  'esl_alloc.c',
  'esl_alphabet.c',
  'esl_arr2.c',
  'esl_arr3.c',
  'esl_avx.c',
  'esl_avx512.c',
  'esl_bitfield.c',
  'esl_buffer.c',
  'esl_cluster.c',
  'esl_composition.c',
  'esl_cpu.c',
  'esl_dirichlet.c',
  'esl_distance.c',
  'esl_dmatrix.c',
  'esl_dsqdata.c',
  'esl_exponential.c',
  'esl_fileparser.c',
  'esl_gamma.c',
  'esl_gencode.c',
  'esl_getopts.c',
  'esl_gev.c',
  'esl_graph.c',
  'esl_gumbel.c',
  'esl_heap.c',
  'esl_histogram.c',
  'esl_hmm.c',
  'esl_huffman.c',
  'esl_hyperexp.c',
  'esl_iset.c',
  'esl_json.c',
  'esl_keyhash.c',
  'esl_lognormal.c',
  'esl_matrixops.c',
  'esl_mem.c',
  'esl_minimizer.c',
  'esl_mixdchlet.c',
  'esl_mixgev.c',
  'esl_mpi.c',
  'esl_msa.c',
  'esl_msacluster.c',
  'esl_msafile.c',
  'esl_msafile2.c',
  'esl_msafile_a2m.c',
  'esl_msafile_afa.c',
  'esl_msafile_clustal.c',
  'esl_msafile_phylip.c',
  'esl_msafile_psiblast.c',
  'esl_msafile_selex.c',
  'esl_msafile_stockholm.c',
  'esl_msashuffle.c',
  'esl_msaweight.c',
  'esl_neon.c',
  'esl_normal.c',
  'esl_paml.c',
  'esl_quicksort.c',
  'esl_rand64.c',
  'esl_random.c',
  'esl_randomseq.c',
  'esl_ratematrix.c',
  'esl_recorder.c',
  'esl_red_black.c',
  'esl_regexp.c',
  'esl_rootfinder.c',
  'esl_scorematrix.c',
  'esl_sq.c',
  'esl_sqio.c',
  destatic.process('esl_sqio_ascii.c'),
  'esl_sqio_ncbi.c',
  'esl_sse.c',
  'esl_ssi.c',
  'esl_stack.c',
  'esl_stats.c',
  'esl_stopwatch.c',
  'esl_stretchexp.c',
  'esl_subcmd.c',
  # 'esl_swat',
  'esl_threads.c',
  'esl_tree.c',
  'esl_varint.c',
  'esl_vectorops.c',
  'esl_vmx.c',
  'esl_weibull.c',
  'esl_workqueue.c',
  'esl_wuss.c',
  'interface_gsl.c',
  'interface_lapack.c',
]

easel_lib = static_library(
  'easel',
  easel_sources,
  dependencies : [math, threads, mpi],
  objects : easel_config_h
)

easel_dep = declare_dependency(
  link_with : easel_lib,
  sources : easel_config_h,
  include_directories : '.',
)