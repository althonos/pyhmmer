API Reference
==============

.. toctree::
   :hidden:

   hmmer <hmmer>
   easel <easel>
   plan7 <plan7>
   errors <errors>


.. currentmodule:: pyhmmer

.. automodule:: pyhmmer


.. only:: html

    HMMER
    -----

    .. autosummary::
        :nosignatures:

        pyhmmer.hmmer.hmmsearch
        pyhmmer.hmmer.phmmer
        pyhmmer.hmmer.nhmmer
        pyhmmer.hmmer.hmmpress


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
       pyhmmer.easel.SequenceFile


    Alignments
    ^^^^^^^^^^

    .. autosummary::
       :nosignatures:

       pyhmmer.easel.MSA
       pyhmmer.easel.TextMSA
       pyhmmer.easel.DigitalMSA
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

    Hidden Markov Model
    ^^^^^^^^^^^^^^^^^^^

    .. autosummary::
        :nosignatures:

        pyhmmer.plan7.HMM
        pyhmmer.plan7.HMMFile


    Profile
    ^^^^^^^

    .. autosummary::
        :nosignatures:

        pyhmmer.plan7.Profile
        pyhmmer.plan7.OptimizedProfile
        pyhmmer.plan7.Background


    Pipelines
    ^^^^^^^^^

    .. autosummary::
        :nosignatures:

        pyhmmer.plan7.Pipeline
        pyhmmer.plan7.Builder


    Results
    ^^^^^^^

    .. autosummary::
        :nosignatures:

        pyhmmer.plan7.TopHits
        pyhmmer.plan7.Hit
        pyhmmer.plan7.Domains
        pyhmmer.plan7.Domain
        pyhmmer.plan7.Alignment


    Miscellaneous
    ^^^^^^^^^^^^^

    .. autosummary::
        :nosignatures:

        pyhmmer.plan7.Cutoffs
        pyhmmer.plan7.EvalueParameters
        pyhmmer.plan7.Offsets


    Errors
    ------

    .. autosummary::
       :nosignatures:

       pyhmmer.errors.AllocationError
       pyhmmer.errors.UnexpectedError
       pyhmmer.errors.EaselError
