import typing

from ..easel import DigitalSequence, DigitalSequenceBlock, MSA
from ..plan7 import HMM, TraceAligner

# --- hmmalign ---------------------------------------------------------------

def hmmalign(
    hmm: HMM,
    sequences: typing.Iterable[DigitalSequence],
    *,
    digitize: bool = False,
    trim: bool = False,
    all_consensus_cols: bool = True,
) -> MSA:
    """Align several sequences to a reference HMM, and return the MSA.

    Arguments:
        hmm (`~pyhmmer.plan7.HMM`): The reference HMM to use for the
            alignment.
        sequences (iterable of `~pyhmmer.easel.DigitalSequence`): The
            sequences to align to the HMM. If you plan on using the
            same sequences several times, consider storing them into
            a `~pyhmmer.easel.DigitalSequenceBlock` directly.
        trim (`bool`): Trim off any residues that get assigned to
            flanking :math:`N` and :math:`C` states (in profile traces)
            or :math:`I_0` and :math:`I_m` (in core traces).
        digitize (`bool`): If set to `True`, returns a `DigitalMSA`
            instead of a `TextMSA`.
        all_consensus_cols (`bool`): Force a column to be created for
            every consensus column in the model, even if it means having
            all gap character in a column.

    Returns:
        `~pyhmmer.easel.MSA`: A multiple sequence alignment containing
        the aligned sequences, either a `TextMSA` or a `DigitalMSA`
        depending on the value of the ``digitize`` argument.

    See Also:
        The `~pyhmmer.plan7.TraceAligner` class, which lets you inspect the
        intermediate tracebacks obtained for each alignment before building
        a MSA.

    .. versionadded:: 0.4.7

    """
    aligner = TraceAligner()
    if not isinstance(sequences, DigitalSequenceBlock):
        sequences = DigitalSequenceBlock(hmm.alphabet, sequences)
    traces = aligner.compute_traces(hmm, sequences)
    return aligner.align_traces(
        hmm,
        sequences,
        traces,
        trim=trim,
        digitize=digitize,
        all_consensus_cols=all_consensus_cols,
    )

