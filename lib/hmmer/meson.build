hmmer_conf = configuration_data({
  'HMMER_VERSION'     : '"3.4"',
  'HMMER_DATE'        : '"Aug 2023"',
  'HMMER_LICENSE'     : '"Freely distributed under the BSD open source license."',
  'HMMER_URL'         : '"http://hmmer.org/"',
  'HMMER_COPYRIGHT'   : '"Copyright (C) 2023 Howard Hughes Medical Institute."',
  'p7_RAMLIMIT'       : '32',
  'p7_NCPU'           : '"2"',
  'p7_ALILENGTH'      : '50',
  'p7_ETARGET_AMINO'  : '0.59',
  'p7_ETARGET_DNA'    : '0.62',
  'p7_ETARGET_OTHER'  : '1.0',
  'p7_SEQDBENV'       : '"BLASTDB"',
  'p7_HMMDBENV'       : '"PFAMDB"',
  'p7_MAXABET'        : '20',
  'p7_MAXCODE'        : '29',
  'p7_MAX_SC_TXTLEN'  : '11',
  'p7_MAXDCHLET'      : '20',
})

# --- Select compute kernel ---

if target_machine.cpu_family() == 'x86' or target_machine.cpu_family() == 'x86_64'
    subdir('impl_sse')
elif target_machine.cpu_family() == 'arm' or target_machine.cpu_family() == 'aarch64'
    hmmer_impl = 'NEON'
    hmmer_conf.set('eslENABLE_NEON', true)
    impl_sources = []
elif target_machine.cpu_family() == 'ppc64' and target_machine.endian() == 'big'
    hmmer_impl = 'VMX'
    hmmer_conf.set('eslENABLE_VMX', true)
    impl_sources = []
else
    error('Failed to select vector implementation for target machine ' + target_machine.cpu_family())
endif

# --- Write header ---

hmmer_p7_config_h = configure_file(
  output : 'p7_config.h', 
  configuration : hmmer_conf
)

# --- Compile library ---

hmmer_sources = [
  'build.c', 
  'cachedb.c', 
  'cachedb_shard.c', 
  'emit.c', 
  'errors.c', 
  'evalues.c', 
  'eweight.c', 
  'generic_decoding.c', 
  'generic_fwdback.c', 
  'generic_fwdback_chk.c', 
  'generic_fwdback_banded.c', 
  'generic_null2.c', 
  'generic_msv.c', 
  'generic_optacc.c', 
  'generic_stotrace.c', 
  'generic_viterbi.c', 
  'generic_vtrace.c', 
  'h2_io.c', 
  'heatmap.c', 
  'hmmlogo.c', 
  'hmmdmstr.c', 
  'hmmdmstr_shard.c', 
  'hmmd_search_status.c', 
  'hmmdwrkr.c', 
  'hmmdwrkr_shard.c', 
  'hmmdutils.c', 
  'hmmer.c', 
  'logsum.c', 
  'modelconfig.c', 
  'modelstats.c', 
  'mpisupport.c', 
  'seqmodel.c', 
  'tracealign.c', 
  'p7_alidisplay.c', 
  'p7_bg.c', 
  'p7_builder.c', 
  'p7_domain.c', 
  'p7_domaindef.c', 
  'p7_gbands.c', 
  'p7_gmx.c', 
  'p7_gmxb.c', 
  'p7_gmxchk.c', 
  'p7_hit.c', 
  'p7_hmm.c', 
  'p7_hmmcache.c', 
  'p7_hmmd_search_stats.c', 
  destatic.process('p7_hmmfile.c'), 
  'p7_hmmwindow.c', 
  'p7_pipeline.c', 
  'p7_prior.c', 
  'p7_profile.c', 
  'p7_spensemble.c', 
  destatic.process('p7_tophits.c'), 
  'p7_trace.c', 
  'p7_scoredata.c', 
  'hmmpgmd2msa.c', 
  'fm_alphabet.c', 
  'fm_general.c', 
  'fm_sse.c', 
  'fm_ssv.c'
]

hmmer_lib = static_library(
  'hmmer',
  [hmmer_sources, impl_sources],
  dependencies : [math, threads, mpi, easel_dep],
  sources : hmmer_p7_config_h,
)

hmmer_dep = declare_dependency(
  link_with : hmmer_lib,
  sources : hmmer_p7_config_h,
  include_directories : '.',
)