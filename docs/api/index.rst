API Reference
==============

.. toctree::
   :hidden:

   hmmer <hmmer>
   easel <easel>
   plan7 <plan7>
   daemon <daemon>
   errors <errors>


.. currentmodule:: pyhmmer

.. automodule:: pyhmmer


.. only:: html

    HMMER
    -----

    .. autosummary::
        :nosignatures:

        pyhmmer.hmmer.hmmsearch
        pyhmmer.hmmer.hmmscan
        pyhmmer.hmmer.phmmer
        pyhmmer.hmmer.nhmmer
        pyhmmer.hmmer.hmmpress
        pyhmmer.hmmer.hmmalign


    Easel
    -----

    Data Structures
    ^^^^^^^^^^^^^^^

    .. autosummary::
       :nosignatures:

       pyhmmer.easel.Bitfield
       pyhmmer.easel.KeyHash


    Sequences
    ^^^^^^^^^

    .. autosummary::
       :nosignatures:

       pyhmmer.easel.Sequence
       pyhmmer.easel.TextSequence
       pyhmmer.easel.DigitalSequence


    Sequence Blocks
    ^^^^^^^^^^^^^^^

    .. autosummary::
       :nosignatures:

       pyhmmer.easel.SequenceBlock
       pyhmmer.easel.TextSequenceBlock
       pyhmmer.easel.DigitalSequenceBlock


    Sequence File
    ^^^^^^^^^^^^^

    .. autosummary::
       :nosignatures:

       pyhmmer.easel.SequenceFile


    Alignments
    ^^^^^^^^^^

    .. autosummary::
       :nosignatures:

       pyhmmer.easel.MSA
       pyhmmer.easel.TextMSA
       pyhmmer.easel.DigitalMSA


    Alignment File
    ^^^^^^^^^^^^^^

    .. autosummary::
       :nosignatures:

       pyhmmer.easel.MSAFile


    Linear Algebra
    ^^^^^^^^^^^^^^

    .. autosummary::
       :nosignatures:

       pyhmmer.easel.Vector
       pyhmmer.easel.VectorF
       pyhmmer.easel.VectorU8
       pyhmmer.easel.Matrix
       pyhmmer.easel.MatrixF
       pyhmmer.easel.MatrixU8


    Miscellaneous
    ^^^^^^^^^^^^^

    .. autosummary::
       :nosignatures:

       pyhmmer.easel.Alphabet
       pyhmmer.easel.Randomness
       pyhmmer.easel.SSIReader
       pyhmmer.easel.SSIWriter



    Plan7
    -----

    Profile HMMs
    ^^^^^^^^^^^^

    .. autosummary::
        :nosignatures:

        pyhmmer.plan7.HMM
        pyhmmer.plan7.Profile
        pyhmmer.plan7.OptimizedProfile
        pyhmmer.plan7.OptimizedProfileBlock


    HMM Files
    ^^^^^^^^^

    .. autosummary::
        :nosignatures:

        pyhmmer.plan7.HMMFile
        pyhmmer.plan7.HMMPressedFile


    Pipelines
    ^^^^^^^^^

    .. autosummary::
        :nosignatures:

        pyhmmer.plan7.Pipeline
        pyhmmer.plan7.LongTargetsPipeline
        pyhmmer.plan7.Builder
        pyhmmer.plan7.Background


    Results
    ^^^^^^^

    .. autosummary::
        :nosignatures:

        pyhmmer.plan7.TopHits
        pyhmmer.plan7.Hit
        pyhmmer.plan7.Domains
        pyhmmer.plan7.Domain
        pyhmmer.plan7.Alignment


    Traces
    ^^^^^^^

    .. autosummary::
        :nosignatures:

        pyhmmer.plan7.TraceAligner
        pyhmmer.plan7.Traces
        pyhmmer.plan7.Trace


    Iterative Searches
    ^^^^^^^^^^^^^^^^^^

    .. autosummary::
       :nosignatures:

       pyhmmer.plan7.IterativeSearch
       pyhmmer.plan7.IterationResult


    Miscellaneous
    ^^^^^^^^^^^^^

    .. autosummary::
        :nosignatures:

        pyhmmer.plan7.Cutoffs
        pyhmmer.plan7.EvalueParameters
        pyhmmer.plan7.Offsets
        pyhmmer.plan7.Transitions


    Daemon
    ------

    .. autosummary::
        :nosignatures:

        pyhmmer.daemon.Client
        pyhmmer.daemon.IterativeSearch


    Errors
    ------

    .. autosummary::
       :nosignatures:

       pyhmmer.errors.AllocationError
       pyhmmer.errors.UnexpectedError
       pyhmmer.errors.EaselError
       pyhmmer.errors.AlphabetMismatch
       pyhmmer.errors.ServerError
       pyhmmer.errors.MissingCutoffs
       pyhmmer.errors.InvalidParameter
       pyhmmer.errors.InvalidHMM
