if HMMER_IMPL == "VMX":
    from libhmmer.impl_vmx.io cimport (
        p7_oprofile_Write, 
        p7_oprofile_ReadMSV, 
        p7_oprofile_ReadInfoMSV,
        p7_oprofile_ReadRest,
        p7_oprofile_Position,
    )
elif HMMER_IMPL == "NEON":
    from libhmmer.impl_neon.io cimport (
        p7_oprofile_Write, 
        p7_oprofile_ReadMSV, 
        p7_oprofile_ReadInfoMSV,
        p7_oprofile_ReadRest,
        p7_oprofile_Position,
    )
elif HMMER_IMPL == "SSE":
    from libhmmer.impl_sse.io cimport (
        p7_oprofile_Write, 
        p7_oprofile_ReadMSV, 
        p7_oprofile_ReadInfoMSV,
        p7_oprofile_ReadRest,
        p7_oprofile_Position,
    )
   