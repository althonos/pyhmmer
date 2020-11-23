# coding: utf-8
# cython: language_level=3, linetrace=True
"""Common errors and status codes for the `easel` and `hmmer` modules.
"""

cimport libeasel

statuscode = {
    libeasel.eslOK: "eslOK",
    libeasel.eslFAIL: "eslFAIL",
    libeasel.eslEOL: "eslEOL",
    libeasel.eslEOF: "eslEOF",
    libeasel.eslEOD: "eslEOD",
    libeasel.eslEMEM: "eslEMEM",
    libeasel.eslENOTFOUND: "eslENOTFOUND",
    libeasel.eslEFORMAT: "eslEFORMAT",
    libeasel.eslEAMBIGUOUS: "eslEAMBIGUOUS",
    libeasel.eslEINCOMPAT: "eslEINCOMPAT",
    libeasel.eslEINVAL: "eslEINVAL",
    libeasel.eslESYS: "eslESYS",
    libeasel.eslECORRUPT: "eslECORRUPT",
    libeasel.eslEINCONCEIVABLE: "eslEINCONCEIVABLE",
    libeasel.eslESYNTAX: "eslESYNTAX",
    libeasel.eslERANGE: "eslERANGE",
    libeasel.eslEDUP: "eslEDUP",
    libeasel.eslENORESULT: "eslENORESULT",
    libeasel.eslETYPE: "eslETYPE",
    libeasel.eslEOVERWRITE: "eslEOVERWRITE",
    libeasel.eslENOSPACE: "eslENOSPACE",
    libeasel.eslEUNIMPLEMENTED: "eslEUNIMPLEMENTED",
    libeasel.eslENOFORMAT: "eslENOFORMAT",
    libeasel.eslENOALPHABET: "eslENOALPHABET",
    libeasel.eslEWRITE: "eslEWRITE",
    libeasel.eslEINACCURATE: "eslEINACCURATE",
}


class UnexpectedError(RuntimeError):
    """An unexpected error that happened in the C code.

    As a user of this library, you should never see this exception being
    raised. If you do, please open an issue with steps to reproduce on the
    `bug tracker <https://github.com/althonos/pyhmmer/issues>`_, so that
    proper error handling can be added to the relevant part of the bindings.

    """

    def __init__(self, int code, str function):
        self.code = code
        self.function = function

    def __repr__(self):
        return "{}({!r}, {!r})".format(type(self).__name__, self.code, self.function)

    def __str__(self):
        return "Unexpected error occurred in {!r}: {} (status code {}).".format(
            self.function,
            statuscode.get(self.code, "<unknown>"),
            self.code,
        )


class AllocationError(MemoryError):
    """A memory error that is caused by an unsuccessful allocation.
    """

    def __init__(self, str ctype):
        self.ctype = ctype

    def __repr__(self):
        return "{}({!r})".format(type(self).__name__, self.ctype)

    def __str__(self):
        return "Could not allocate {}.".format(self.ctype)


class EaselError(RuntimeError):
    """An error that was raised from the Easel code.
    """

    def __init__(self, int code, str message):
        self.code = code
        self.message = message

    def __repr__(self):
        return "{}({!r}, {!r})".format(type(self).__name__, self.code, self.message)

    def __str__(self):
        return "Error raised from C code: {}, {} (status code {}).".format(
            statuscode.get(self.code, "<unknown>"),
            self.message,
            self.code
        )
