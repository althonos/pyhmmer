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

    As an library user, you should not see this exception being raised. If you
    do, please open an issue on the bug tracker with steps to reproduce, so
    that proper error handling can be added to the relevant part of the code.

    """

    def __init__(self, code, function):
        self.code = code
        self.function = function

    def __repr__(self):
        return "{}({!r}, {!r})".format(type(self).__name__, self.code, self.function)

    def __str__(self):
        return "Unexpected error occurred in {!r}: {} (status code {}).".format(
            statuscode.get(self.code, "<unknown>"),
            self.code,
            self.function,
        )


class AllocationError(MemoryError):
    """A memory error that is caused by an unsuccessful allocation.
    """

    def __init__(self, ctype):
        self.ctype = ctype

    def __repr__(self):
        return "{}({!r})".format(type(self).__name__, self.ctype)

    def __str__(self):
        return "Could not allocate {}.".format(self.ctype)
