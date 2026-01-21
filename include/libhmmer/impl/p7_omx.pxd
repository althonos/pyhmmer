if HMMER_IMPL == "VMX":
    from libhmmer.impl_vmx.p7_omx cimport (
        P7_OMX, 
        p7_omx_Create,
        p7_omx_GrowTo,
        p7_omx_Reuse,
        p7_omx_Destroy,
    )
elif HMMER_IMPL == "NEON":
    from libhmmer.impl_neon.p7_omx cimport (
        P7_OMX, 
        p7_omx_Create,
        p7_omx_GrowTo,
        p7_omx_Reuse,
        p7_omx_Destroy,
    )
elif HMMER_IMPL == "SSE":
    from libhmmer.impl_sse.p7_omx cimport (
        P7_OMX, 
        p7_omx_Create,
        p7_omx_GrowTo,
        p7_omx_Reuse,
        p7_omx_Destroy,
    )
   
