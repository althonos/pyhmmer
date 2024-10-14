# coding: utf-8
"""Reimplementation of HMMER binaries with the PyHMMER API.

Note:
    Functions of this module handle parallelization using threads to run
    searches in parallel for the different queries. If less queries are
    given, the number of threads will be reduced to avoid spawning idle
    threads.

"""
# import abc
# import contextlib
# import collections
# import copy
# import ctypes
# import itertools
# import io
# import multiprocessing
# import os
# import operator
# import queue
# import sys
# import threading
# import time
# import typing
# from typing import Any

# import psutil

# from ..easel import (
#     Alphabet,
#     DigitalSequence,
#     DigitalMSA,
#     MSA,
#     MSAFile,
#     TextSequence,
#     SequenceFile,
#     SSIWriter,
#     DigitalSequenceBlock,
# )
# from ..plan7 import (
#     Builder,
#     Background,
#     Pipeline,
#     LongTargetsPipeline,
#     TopHits,
#     IterationResult,
#     HMM,
#     HMMFile,
#     HMMPressedFile,
#     Profile,
#     TraceAligner,
#     OptimizedProfile,
#     OptimizedProfileBlock,
# )
# from ..utils import peekable, singledispatchmethod

from ._hmmalign import hmmalign
from ._hmmpress import hmmpress
from ._hmmscan import hmmscan
from ._hmmsearch import hmmsearch
from ._jackhmmer import jackhmmer
from ._nhmmer import nhmmer
from ._phmmer import phmmer










# add a very limited CLI so that this module can be invoked in a shell:
#     $ python -m pyhmmer.hmmsearch <hmmfile> <seqdb>
if __name__ == "__main__":
    import argparse
    import sys

    # don't load target databases in memory if they would take more than
    # 80% of the remainining available memory (because more memory will be
    # needed afterwards to allocate the `TopHits` for each query)
    MAX_MEMORY_LOAD = 0.80

    def _hmmsearch(args: argparse.Namespace) -> int:
        # check the size of the target database and the amount of available memory
        available_memory = psutil.virtual_memory().available
        database_size = os.stat(args.seqdb).st_size

        with SequenceFile(args.seqdb, digital=True) as sequences:
            # pre-load the database if it is small enough to fit in memory
            if database_size < available_memory * MAX_MEMORY_LOAD:
                sequences = sequences.read_block()  # type: ignore
            # load the query HMMs iteratively
            with HMMFile(args.hmmfile) as hmms:
                if hmms.is_pressed():
                    hmms = hmms.optimized_profiles()  # type: ignore
                hits_list = hmmsearch(hmms, sequences, cpus=args.jobs)
                for hits in hits_list:
                    for hit in hits:
                        if hit.reported:
                            print(
                                hit.name.decode(),
                                (hit.accession or b"-").decode(),
                                (hits.query.name or b"-").decode(),
                                (hits.query.accession or b"-").decode(),
                                hit.evalue,
                                hit.score,
                                hit.bias,
                                sep="\t",
                            )

        return 0

    def _phmmer(args: argparse.Namespace) -> int:
        # check the size of the target database and the amount of available memory
        available_memory = psutil.virtual_memory().available
        database_size = os.stat(args.seqdb).st_size

        alphabet = Alphabet.amino()
        with SequenceFile(args.seqdb, digital=True, alphabet=alphabet) as sequences:
            # pre-load the database if it is small enough to fit in memory
            if database_size < available_memory * MAX_MEMORY_LOAD:
                sequences = sequences.read_block()  # type: ignore
            # load the query sequences iteratively
            with SequenceFile(args.seqfile, digital=True, alphabet=alphabet) as queries:
                hits_list = phmmer(queries, sequences, cpus=args.jobs)
                for hits in hits_list:
                    for hit in hits:
                        if hit.reported:
                            print(
                                hit.name.decode(),
                                "-",
                                hit.best_domain.alignment.hmm_accession.decode(),
                                hit.best_domain.alignment.hmm_name.decode(),
                                hit.evalue,
                                hit.score,
                                hit.bias,
                                sep="\t",
                            )

        return 0

    @contextlib.contextmanager
    def open_query_file(
        queryfile: typing.Union["os.PathLike[str]", typing.BinaryIO],
        alphabet: Alphabet,
    ) -> typing.Iterator[typing.Union["SequenceFile[DigitalSequence]", HMMFile]]:
        """Open either a sequence file or an HMM file."""
        try:
            yield SequenceFile(queryfile, digital=True, alphabet=alphabet)
        except ValueError:
            yield HMMFile(queryfile)

    def _jackhmmer(args: argparse.Namespace) -> int:
        # check the size of the target database and the amount of available memory
        available_memory = psutil.virtual_memory().available
        database_size = os.stat(args.seqdb).st_size

        alphabet = Alphabet.amino()
        with SequenceFile(args.seqdb, digital=True, alphabet=alphabet) as sequences:
            # pre-load the database if it is small enough to fit in memory
            if database_size < available_memory * MAX_MEMORY_LOAD:
                sequences = sequences.read_block()  # type: ignore
            # load the query sequences or HMMs iteratively
            with open_query_file(args.queryfile, alphabet) as queries:
                results = jackhmmer(
                    typing.cast(typing.Iterable[HMM], queries),
                    typing.cast(typing.Iterable[DigitalSequence], sequences),
                    checkpoints=False,
                    cpus=typing.cast(int, args.jobs)
                )
                for result in results:
                    for hit in result.hits:
                        if hit.reported:
                            print(
                                hit.name.decode(),
                                "-",
                                hit.best_domain.alignment.hmm_accession.decode(),
                                hit.best_domain.alignment.hmm_name.decode(),
                                hit.evalue,
                                hit.score,
                                hit.bias,
                                sep="\t",
                            )

        return 0

    def _nhmmer(args: argparse.Namespace) -> int:
        # at the moment `LongTargetsPipeline` only support block targets, not files
        with SequenceFile(args.seqdb, digital=True) as seqfile:
            with SequenceFile(args.seqfile, digital=True) as queryfile:
                hits_list = nhmmer(queryfile, seqfile, cpus=args.jobs)
                for hits in hits_list:
                    for hit in hits:
                        if hit.reported:
                            print(
                                hit.name.decode(),
                                "-",
                                hit.best_domain.alignment.hmm_accession.decode(),
                                hit.best_domain.alignment.hmm_name.decode(),
                                hit.evalue,
                                hit.score,
                                hit.bias,
                                sep="\t",
                            )

        return 0

    def _hmmscan(args: argparse.Namespace) -> int:
        # check the size of the target database and the amount of available memory
        available_memory = psutil.virtual_memory().available
        database_size = os.stat(args.hmmdb).st_size

        with SequenceFile(args.seqfile, digital=True) as seqfile:
            with HMMFile(args.hmmdb) as hmms:
                # pre-load profiles is they can fit into memory
                targets = hmms.optimized_profiles() if hmms.is_pressed() else hmms
                if (
                    hmms.is_pressed()
                    and database_size < available_memory * MAX_MEMORY_LOAD
                ):
                    targets = OptimizedProfileBlock(seqfile.alphabet, targets)  # type: ignore
                # load the query sequences iteratively
                hits_list = hmmscan(seqfile, targets, cpus=args.jobs)  # type: ignore
                for hits in hits_list:
                    for hit in hits:
                        if hit.reported:
                            print(
                                hit.name.decode(),
                                (hit.accession or b"-").decode(),
                                (hits.query.name or b"-").decode(),
                                (hits.query.accession or b"-").decode(),
                                hit.evalue,
                                hit.score,
                                hit.bias,
                                sep="\t",
                            )
        return 0

    def _hmmpress(args: argparse.Namespace) -> int:
        for ext in ["h3m", "h3i", "h3f", "h3p"]:
            path = "{}.{}".format(args.hmmfile, ext)
            if os.path.exists(path):
                if args.force:
                    os.remove(path)
                else:
                    print(f"file {path} already exists")
                    return 1

        with HMMFile(args.hmmfile) as hmms:
            hmmpress(hmms, args.hmmfile)

        return 0

    def _hmmalign(args: argparse.Namespace) -> int:
        try:
            with SequenceFile(args.seqfile, args.informat, digital=True) as seqfile:
                sequences: typing.List[DigitalSequence] = list(seqfile)
        except EOFError as err:
            print(err, file=sys.stderr)
            return 1

        with HMMFile(args.hmmfile) as hmms:
            hmm = next(hmms)
            if next(hmms, None) is not None:
                print("HMM file contains more than one HMM, exiting", file=sys.stderr)
                return 1

        msa = hmmalign(hmm, sequences, trim=args.trim)
        if args.output == "-":
            with io.BytesIO() as out:
                msa.write(out, args.outformat)
                print(out.getvalue().decode("ascii"), end="")
        else:
            with open(args.output, "wb") as out:
                msa.write(out, args.outformat)

        return 0

    parser = argparse.ArgumentParser()
    parser.add_argument("-j", "--jobs", required=False, default=0, type=int)
    subparsers = parser.add_subparsers(
        dest="cmd", help="HMMER command to run", required=True
    )

    parser_hmmsearch = subparsers.add_parser("hmmsearch")
    parser_hmmsearch.set_defaults(call=_hmmsearch)
    parser_hmmsearch.add_argument("hmmfile")
    parser_hmmsearch.add_argument("seqdb")

    parser_phmmer = subparsers.add_parser("phmmer")
    parser_phmmer.set_defaults(call=_phmmer)
    parser_phmmer.add_argument("seqfile")
    parser_phmmer.add_argument("seqdb")

    parser_phmmer = subparsers.add_parser("jackhmmer")
    parser_phmmer.set_defaults(call=_jackhmmer)
    parser_phmmer.add_argument("queryfile")  # can be sequences or HMM
    parser_phmmer.add_argument("seqdb")

    parser_nhmmer = subparsers.add_parser("nhmmer")
    parser_nhmmer.set_defaults(call=_nhmmer)
    parser_nhmmer.add_argument("seqfile")
    parser_nhmmer.add_argument("seqdb")

    parser_hmmsearch = subparsers.add_parser("hmmscan")
    parser_hmmsearch.set_defaults(call=_hmmscan)
    parser_hmmsearch.add_argument("hmmdb")
    parser_hmmsearch.add_argument("seqfile")

    parser_hmmpress = subparsers.add_parser("hmmpress")
    parser_hmmpress.set_defaults(call=_hmmpress)
    parser_hmmpress.add_argument("hmmfile")
    parser_hmmpress.add_argument("-f", "--force", action="store_true")

    parser_hmmalign = subparsers.add_parser("hmmalign")
    parser_hmmalign.set_defaults(call=_hmmalign)
    parser_hmmalign.add_argument("hmmfile", metavar="<hmmfile>")
    parser_hmmalign.add_argument(
        "seqfile",
        metavar="<seqfile>",
    )
    parser_hmmalign.add_argument(
        "-o",
        "--output",
        action="store",
        default="-",
        metavar="<f>",
        help="output alignment to file <f>, not stdout",
    )
    parser_hmmalign.add_argument(
        "--trim",
        action="store_true",
        help="trim terminal tails of nonaligned residues from alignment",
    )
    parser_hmmalign.add_argument(
        "--informat",
        action="store",
        metavar="<s>",
        help="assert <seqfile> is in format <s> (no autodetection)",
        choices=SequenceFile._FORMATS.keys(),
    )
    parser_hmmalign.add_argument(
        "--outformat",
        action="store",
        metavar="<s>",
        help="output alignment in format <s>",
        default="stockholm",
        choices=MSAFile._FORMATS.keys(),
    )

    args = parser.parse_args()
    sys.exit(args.call(args))
