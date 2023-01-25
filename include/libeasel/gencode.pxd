from libc.stdio cimport FILE
from libc.stdint cimport int8_t

from libeasel cimport ESL_DSQ
from libeasel.alphabet cimport ESL_ALPHABET
from libeasel.sq cimport ESL_SQ, ESL_SQ_BLOCK


cdef extern from "esl_gencode.h" nogil:

    ctypedef struct ESL_GENCODE:
        int     transl_table
        char[128]  desc

        ESL_DSQ[64] basic
        int8_t[64]  is_initiator

        ESL_ALPHABET* nt_abc
        ESL_ALPHABET* aa_abc

    ctypedef esl_gencode_workstate_s ESL_GENCODE_WORKSTATE
    ctypedef struct esl_gencode_workstate_s:
        # stateful info (which may get updated with each new seq, strand, and/or window):
        ESL_SQ[3] *psq
        int8_t[3]  in_orf
        int     apos
        int     frame
        int     codon
        int     inval
        bint    is_revcomp
        int     orfcount

        ESL_SQ_BLOCK* orf_block

        # one-time configuration information (from options)
        bint     do_watson;
        bint     do_crick;
        bint     using_initiators
        int     minlen
        FILE   *outfp
        int     outformat


    # Create/Destroy workstate
    # extern void esl_gencode_WorkstateDestroy(ESL_GENCODE_WORKSTATE *wrk);
    # extern ESL_GENCODE_WORKSTATE * esl_gencode_WorkstateCreate(ESL_GETOPTS *go, ESL_GENCODE *gcode);

    # the ESL_GENCODE genetic code object
    ESL_GENCODE *esl_gencode_Create(const ESL_ALPHABET *nt_abc, const ESL_ALPHABET *aa_abc);
    void         esl_gencode_Destroy            (ESL_GENCODE *gcode);
    int          esl_gencode_Set                (ESL_GENCODE *gcode,  int ncbi_transl_table);
    int          esl_gencode_SetInitiatorAny    (ESL_GENCODE *gcode);
    int          esl_gencode_SetInitiatorOnlyAUG(ESL_GENCODE *gcode);

    # reading and writing genetic codes in NCBI format
    # int          esl_gencode_Read(ESL_FILEPARSER *efp, const ESL_ALPHABET *nucleic_abc, const ESL_ALPHABET *amino_abc, ESL_GENCODE **ret_gcode);
    # int          esl_gencode_Write(FILE *ofp, const ESL_GENCODE *gcode, int add_comment);

    # DNA->protein digital translation, allowing ambiguity chars
    int   esl_gencode_GetTranslation(const ESL_GENCODE *gcode, ESL_DSQ *dsqp);
    int   esl_gencode_IsInitiator   (const ESL_GENCODE *gcode, ESL_DSQ *dsqp);

    # Debugging/development utilities
    char *esl_gencode_DecodeDigicodon(const ESL_GENCODE *gcode, int digicodon, char *codon);
    int   esl_gencode_DumpAltCodeTable(FILE *ofp);
    int   esl_gencode_Compare(const ESL_GENCODE *gc1, const ESL_GENCODE *gc2, int metadata_too);

    # Functions for processing ORFs
    int esl_gencode_ProcessOrf(ESL_GENCODE_WORKSTATE *wrk, ESL_SQ *sq);
    void esl_gencode_ProcessStart(ESL_GENCODE *gcode, ESL_GENCODE_WORKSTATE *wrk, ESL_SQ *sq);
    int esl_gencode_ProcessPiece(ESL_GENCODE *gcode, ESL_GENCODE_WORKSTATE *wrk, ESL_SQ *sq);
    int esl_gencode_ProcessEnd(ESL_GENCODE_WORKSTATE *wrk, ESL_SQ *sq);
