# coding: utf-8
# cython: language_level=3, linetrace=True
"""Reimplementation of the HMMPGMD client.
"""

# --- C imports --------------------------------------------------------------

from cpython.bytearray cimport PyByteArray_AS_STRING

cimport libeasel
cimport libhmmer.hmmpgmd
cimport libhmmer.p7_hit
from libc.stdlib cimport free, realloc
from libc.string cimport memset, memcpy
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t
from libhmmer.hmmpgmd cimport HMMD_SEARCH_STATUS_SERIAL_SIZE, HMMD_SEARCH_STATUS, HMMD_SEARCH_STATS
from libhmmer.p7_pipeline cimport p7_pipemodes_e, P7_PIPELINE
from libhmmer.p7_hit cimport P7_HIT

from pyhmmer.easel cimport Sequence, Alphabet
from pyhmmer.errors import UnexpectedError, AllocationError
from pyhmmer.plan7 cimport TopHits, Pipeline, HMM


# --- Python imports ---------------------------------------------------------

import io
import socket
import warnings


# --- Cython classes ---------------------------------------------------------

cdef class Client:

    def __init__(
        self,
        str address="127.0.0.1",
        uint16_t port=51371,
    ):
        self.address = address
        self.port = port
        self.socket = socket.socket()

    def __enter__(self):
        self.connect()
        return self

    def __exit__(self, exc_value, exc_type, traceback):
        self.close()

    def connect(self):
        self.socket.connect((self.address, self.port))

    def close(self):
        self.socket.close()

    cdef bytearray _recvall(self, size_t message_size):

        cdef bytearray buffer = bytearray(message_size)
        cdef object    view   = memoryview(buffer)
        cdef size_t received  = 0
        cdef size_t recv_size = 0

        while received < message_size:
            recv_size = self.socket.recv_into(view)
            if recv_size == 0:
                raise EOFError(f"Expected message of size {message_size}, received {received}")
            received += recv_size
            view = view[recv_size:]
        return buffer

    cdef TopHits _client(
        self,
        object query,
        uint64_t db,
        Pipeline pli,
        p7_pipemodes_e mode,
    ):
        cdef int                status
        cdef HMMD_SEARCH_STATS  search_stats
        cdef HMMD_SEARCH_STATUS search_status

        cdef object             send_buffer
        cdef bytearray          response
        cdef const char*        response_data

        cdef uint32_t           hits_start
        cdef uint32_t           buf_offset    = 0
        cdef TopHits            hits          = TopHits()
        cdef str                options       = "".join(pli.arguments())

        # clean memory for data structures allocated on the stack
        memset(&search_stats, 0, sizeof(HMMD_SEARCH_STATS))
        memset(&search_status, 0, sizeof(HMMD_SEARCH_STATUS))
        search_stats.hit_offsets = NULL

        try:
            # serialize the query over the socket
            self.socket.sendall(f"@--seqdb {db} {options}\n".encode("ascii"))
            query.write(self.socket.makefile("wb"))
            self.socket.sendall(b"//")

            # get the search status back
            response = self._recvall(HMMD_SEARCH_STATUS_SERIAL_SIZE)
            status = libhmmer.hmmpgmd.hmmd_search_status_Deserialize(
                <const uint8_t*> PyByteArray_AS_STRING(response),
                &buf_offset,
                &search_status
            )
            if status != libeasel.eslOK:
                raise UnexpectedError(status, "hmmd_search_status_Deserialize")

            # check if error happened
            if search_status.status != libeasel.eslOK:
                error = self.socket.recv(search_status.msg_size)
                raise RuntimeError(error.decode("utf-8", "replace"))

            # get the response
            response = self._recvall(search_status.msg_size)
            response_data = PyByteArray_AS_STRING(response)

            with nogil:
                # deserialize search_stats
                buf_offset = 0
                status = libhmmer.hmmpgmd.p7_hmmd_search_stats_Deserialize(
                    <const uint8_t*> response_data,
                    &buf_offset,
                    &search_stats
                )
                if status != libeasel.eslOK:
                    raise UnexpectedError(status, "p7_hmmd_search_search_stats_Deserialize")

                # copy input parameters from the pipeline
                memcpy(&hits._pli, pli._pli, sizeof(P7_PIPELINE))

                # copy the search search_stats
                hits._pli.mode                = mode
                hits._pli.nmodels             = search_stats.nmodels
                hits._pli.nseqs               = search_stats.nseqs
                hits._pli.n_past_msv          = search_stats.n_past_msv
                hits._pli.n_past_vit          = search_stats.n_past_vit
                hits._pli.n_past_fwd          = search_stats.n_past_fwd
                hits._pli.Z                   = search_stats.Z
                hits._pli.domZ                = search_stats.domZ
                hits._pli.Z_setby             = search_stats.Z_setby
                hits._pli.domZ_setby          = search_stats.domZ_setby
                hits._th.nreported            = search_stats.nreported
                hits._th.nincluded            = search_stats.nincluded
                hits._th.is_sorted_by_seqidx  = False
                hits._th.is_sorted_by_sortkey = True
                hits._th.N                    = search_stats.nhits

                # reallocate hit arrays
                if search_stats.nhits > 0:
                    hits._th.unsrt = <P7_HIT*> realloc(hits._th.unsrt, search_stats.nhits * sizeof(P7_HIT))
                    if hits._th.unsrt == NULL:
                        raise AllocationError("P7_HIT", sizeof(P7_HIT), search_stats.nhits)
                    hits._th.hit = <P7_HIT**> realloc(hits._th.hit, search_stats.nhits * sizeof(P7_HIT*))
                    if hits._th.hit == NULL:
                        raise AllocationError("P7_HIT*", sizeof(P7_HIT*), search_stats.nhits)

                # deserialize hits
                hits_start = buf_offset
                for i in range(search_stats.nhits):
                    # clean pointers in hit data to force reallocation
                    hits._th.unsrt[i].name = NULL
                    hits._th.unsrt[i].acc  = NULL
                    hits._th.unsrt[i].desc = NULL
                    hits._th.unsrt[i].dcl  = NULL
                    # check buffers match the ones in the search stats
                    if buf_offset - hits_start != search_stats.hit_offsets[i]:
                        with gil:
                            warnings.warn(f"Hit offset {i} did not match expected (expected {search_stats.hit_offsets[i]}, found {buf_offset - hits_start})")
                    # deserialize and record the hit
                    status = libhmmer.p7_hit.p7_hit_Deserialize(
                        <const uint8_t*> response_data,
                        &buf_offset,
                        &hits._th.unsrt[i]
                    )
                    if status != libeasel.eslOK:
                        raise UnexpectedError(status, "p7_hit_Deserialize")
                    hits._th.hit[i] = &hits._th.unsrt[i]

        finally:
            free(search_stats.hit_offsets)

        return hits

    def search_seq(self, Sequence query, uint64_t db = 1, **options):
        cdef Alphabet abc = getattr(query, "alphabet", Alphabet.amino())
        cdef Pipeline pli = Pipeline(abc, **options)
        return self._client(query, db, pli, p7_pipemodes_e.p7_SEARCH_SEQS)

    def search_hmm(self, HMM query, uint64_t db = 1, **options):
        cdef Pipeline pli = Pipeline(query.alphabet, **options)
        return self._client(query, db, pli, p7_pipemodes_e.p7_SEARCH_SEQS)
