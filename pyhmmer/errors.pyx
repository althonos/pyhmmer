# coding: utf-8
# cython: language_level=3
"""Common errors and status codes for the `easel` and `hmmer` modules.
"""

cimport libeasel
cimport libeasel.alphabet

include "_getid.pxi"

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
        return "Unexpected error occurred in {!r}: {} (status code {})".format(
            self.function,
            statuscode.get(self.code, "<unknown>"),
            self.code,
        )


class AllocationError(MemoryError):
    """A memory error that is caused by an unsuccessful allocation.
    """

    def __init__(self, str ctype, size_t itemsize, size_t count=1):
        self.ctype = ctype
        self.itemsize = itemsize
        self.count = count

    def __repr__(self):
        cdef str ty = type(self).__name__
        if self.count == 1:
            return f"{ty}({self.ctype!r}, {self.itemsize})"
        return f"{ty}({self.ctype!r}, {self.itemsize}, {self.count})"

    def __str__(self):
        if self.count == 1:
            return f"Could not allocate {self.itemsize} bytes for type {self.ctype}"
        return f"Could not allocate {self.itemsize*self.count} bytes for an array of {self.count} {self.ctype}"


class EaselError(RuntimeError):
    """An error that was raised from the Easel code.
    """

    def __init__(self, int code, str message):
        self.code = code
        self.message = message

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}({self.code!r}, {self.message!r})"

    def __str__(self):
        return "Error raised from C code: {}, {} (status code {})".format(
            self.message,
            statuscode.get(self.code, "<unknown>"),
            self.code
        )


class AlphabetMismatch(ValueError):
    """A value error caused by an alphabet mismatch.

    This error is raised on occasions where several arguments to a
    function have mismatching alphabets. For instance, passing a
    `~pyhmmer.easel.DigitalSequence` to a `~pyhmmer.plan7.Builder` with
    different alphabets::

        >>> seq = easel.DigitalSequence(easel.Alphabet.dna())
        >>> builder = plan7.Builder(easel.Alphabet.amino())
        >>> builder.build(seq, plan7.Background(seq.alphabet))
        Traceback (most recent call last):
          ...
        pyhmmer.errors.AlphabetMismatch: Expected ..., found ...

    .. versionadded:: 0.3.0

    """

    def __init__(self, expected, actual):
        super().__init__(self)
        self.expected = expected
        self.actual = actual

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}({self.expected!r}, {self.actual!r})"

    def __str__(self):
        return f"Expected {self.expected.type} alphabet, found {self.actual.type} alphabet"

    def __eq__(self, other):
        if not isinstance(other, AlphabetMismatch):
            return NotImplemented
        return self.actual == other.actual and self.expected == other.expected


class ServerError(RuntimeError):
    """A runtime error that happened in a ``hmmpgmd`` server.
    """

    def __init__(self, int code, str message):
        self.code = code
        self.message = message

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}({self.code!r}, {self.message!r})"

    def __str__(self):
        return "Error in server: {}, {} (status code {})".format(
            self.message,
            statuscode.get(self.code, "<unknown>"),
            self.code,
        )


class MissingCutoffs(ValueError):
    """The model is missing bitscore cutoffs required by the pipeline.

    .. versionadded:: 0.8.2

    .. versionadded:: 0.10.8
       The ``model_name`` and ``bit_cutoffs`` displayed in the error message.

    """

    def __init__(self, str model_name = None, str bit_cutoffs = None):
        self.model_name = model_name
        self.bit_cutoffs = bit_cutoffs

    def __repr__(self):
        cdef str ty = type(self).__name__
        args = []
        if self.model_name is not None:
            args.append(f"model_name={self.model_name!r}")
        if self.bit_cutoffs is not None:
            args.append(f"bit_cutoffs={self.bit_cutoffs!r}")
            return f"{ty}({self.model_name!r})"
        return f"{ty}({ ', '.join(args) })"

    def __str__(self):
        model = f"Model {self.model_name!r}" if self.model_name is not None else "Model"
        bit_cutoffs = f"{self.bit_cutoffs!r} bitscore cutoff" if self.bit_cutoffs is not None else f"bitscore cutoff"
        return f"{model} is missing {bit_cutoffs} required by pipeline"


class InvalidParameter(ValueError):
    """An invalid parameter value was given.

    .. versionadded:: 0.8.2

    """

    def __init__(self, str name, object value, *, list choices=None, str hint=None):
        self.name = name
        self.value = value
        self.choices = choices
        self.hint = hint

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}({self.name!r}, {self.value!r}, {self.choices!r})"

    def __str__(self):
        options = ""
        hint = self.hint

        if hint is None and self.choices is not None:
            choices = [
                x.__name__
                if isinstance(x, type)
                else repr(x)
                for x in self.choices
            ]
            hint = " or ".join([
                ", ".join(choices[:len(self.choices)-1]),
                choices[len(self.choices)-1]
            ])
        if hint is not None:
            options = f" (expected {hint})"

        return f"Invalid {self.name!r} parameter value: {self.value!r}{options}"


class InvalidHMM(ValueError):
    """A value error caused by a HMM that fails validation.

    .. versionadded:: 0.10.3

    """

    def __init__(self, object hmm, str message):
        super().__init__(hmm, message)
        self.message = message
        self.hmm = hmm

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}({self.code!r}, {self.message!r})"

    def __str__(self):
        return "HMM {!r} is invalid: {}".format(
            self.hmm.name.decode('utf-8', 'replace'),
            self.message
        )

