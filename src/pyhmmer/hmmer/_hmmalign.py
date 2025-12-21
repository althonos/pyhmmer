import os
import typing
import functools
import itertools
from multiprocessing.pool import ThreadPool

import psutil

from ..easel import DigitalSequence, DigitalSequenceBlock, MSA
from ..plan7 import HMM, TraceAligner, Traces

if typing.TYPE_CHECKING:
    from ._base import BACKEND

# --- Utils ------------------------------------------------------------------

def _batch_sequences(sequences: DigitalSequenceBlock, cpus: int) -> typing.Iterable[DigitalSequenceBlock]:
    batchsize = (len(sequences) + cpus - 1) // cpus
    indices = [0]
    for i in range(cpus):
        indices.append(indices[i] + batchsize)
    indices[-1] = len(sequences)
    return (sequences[indices[i]:indices[i+1]] for i in range(cpus))

# --- hmmalign ---------------------------------------------------------------

def hmmalign(
    hmm: HMM,
    sequences: typing.Iterable[DigitalSequence],
    *,
    cpus: int = 0,
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

    Keyword Arguments:
        cpus (`int`): The number of threads to run in parallel. Pass ``1``
            to run everything in the main thread, ``0`` to automatically
            select a suitable number (using `psutil.cpu_count`), or any
            positive number otherwise.
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

    .. versionadded:: 0.11.4
        The ``cpus`` argument for parallel processing with multithreading.

    """
    cpus = cpus if cpus > 0 else psutil.cpu_count(logical=False) or os.cpu_count() or 1

    aligner = TraceAligner()
    if not isinstance(sequences, DigitalSequenceBlock):
        sequences = DigitalSequenceBlock(hmm.alphabet, sequences)

    if len(sequences) < cpus and cpus > 1:
        cpus = len(sequences)

    if cpus == 1:
        traces = aligner.compute_traces(hmm, sequences)
    else:
        batches = _batch_sequences(sequences, cpus)
        with ThreadPool(cpus) as pool:
            traces = Traces()
            _compute_traces = functools.partial(aligner.compute_traces, hmm)
            for batch_traces in pool.map(_compute_traces, batches):
                traces.extend(batch_traces)

    return aligner.align_traces(
        hmm,
        sequences,
        traces,
        trim=trim,
        digitize=digitize,
        all_consensus_cols=all_consensus_cols,
    )

