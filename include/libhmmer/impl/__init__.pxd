from . cimport io, p7_omx, p7_oprofile

if HMMER_IMPL == "VMX":
    from libhmmer.impl_vmx cimport (
        p7O_EXTRA_SB,
        impl_Init,
        p7_MSVFilter,
    )
elif HMMER_IMPL == "NEON":
    from libhmmer.impl_neon cimport (
        p7O_EXTRA_SB,
        impl_Init,
        p7_SSVFilter,
        p7_MSVFilter,
    )
elif HMMER_IMPL == "SSE":
    from libhmmer.impl_sse cimport (
        p7O_EXTRA_SB,
        impl_Init,
        p7_SSVFilter,
        p7_MSVFilter,
    )
   