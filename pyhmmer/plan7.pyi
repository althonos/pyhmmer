# coding: utf-8
import collections.abc
import os
import types
import typing

from .easel import Alphabet, Sequence

class Alignment(collections.abc.Sized):
    domain: Domain
    def __len__(self) -> int: ...
    @property
    def hmm_from(self) -> int: ...
    @property
    def hmm_name(self) -> bytes: ...
    @property
    def hmm_sequence(self) -> str: ...
    @property
    def hmm_to(self) -> int: ...
    @property
    def target_from(self) -> int: ...
    @property
    def target_name(self) -> bytes: ...
    @property
    def target_sequence(self) -> str: ...
    @property
    def target_to(self) -> int: ...
    @property
    def identity_sequence(self) -> str: ...

class Background(object):
    def __init__(self, alphabet: Alphabet, uniform: bool = False) -> None: ...
    def __copy__(self) -> Background: ...
    @property
    def L(self) -> int: ...
    @L.setter
    def L(self, L: int) -> None: ...
    def copy(self) -> Background: ...

class Domain(object):
    alignment: Alignment
    hit: Hit
    @property
    def env_from(self) -> int: ...
    @property
    def env_to(self) -> int: ...
    @property
    def ali_from(self) -> int: ...
    @property
    def ali_to(self) -> int: ...
    @property
    def hmm_from(self) -> int: ...
    @property
    def hmm_to(self) -> int: ...
    @property
    def score(self) -> float: ...
    @property
    def bias(self) -> float: ...
    @property
    def correction(self) -> float: ...
    @property
    def envelope_score(self) -> float: ...
    @property
    def c_evalue(self) -> float: ...
    @property
    def i_evalue(self) -> float: ...

class Domains(typing.Sequence[Domain]):
    hit: Hit
    def __len__(self) -> int: ...
    @typing.overload
    def __getitem__(self, index: int) -> Domain: ...
    @typing.overload
    def __getitem__(self, index: slice) -> typing.Sequence[Domain]: ...

class Hit(object):
    hits: TopHits
    @property
    def name(self) -> bytes: ...
    @property
    def accession(self) -> bytes: ...
    @property
    def description(self) -> bytes: ...
    @property
    def score(self) -> float: ...
    @property
    def pre_score(self) -> float: ...
    @property
    def bias(self) -> float: ...
    @property
    def evalue(self) -> float: ...
    @property
    def domains(self) -> Domains: ...
    def is_included(self) -> bool: ...
    def is_reported(self) -> bool: ...
    def is_new(self) -> bool: ...
    def is_dropped(self) -> bool: ...
    def is_duplicate(self) -> bool: ...

class HMM(object):
    alphabet: Alphabet
    def __init__(self, M: int, alphabet: Alphabet) -> None: ...
    @property
    def M(self) -> int: ...
    @property
    def name(self) -> bytes: ...
    @name.setter
    def name(self, names: bytes) -> None: ...
    @property
    def accession(self) -> bytes: ...
    @accession.setter
    def accession(self, accession: bytes) -> None: ...
    @property
    def description(self) -> bytes: ...
    @description.setter
    def description(self, description: bytes) -> None: ...
    def write(self, fh: typing.BinaryIO, binary: bool = False) -> None: ...
    def zero(self) -> None: ...

class HMMFile(typing.ContextManager[HMMFile], typing.Iterator[HMM]):
    def __init__(
        self,
        file: typing.Union[typing.AnyStr, os.PathLike[typing.AnyStr], typing.BinaryIO],
        db: bool = True
    ) -> None: ...
    def __enter__(self) -> HMMFile: ...
    def __exit__(
        self,
        exc_type: typing.Optional[typing.Type[BaseException]],
        exc_value: typing.Optional[BaseException],
        traceback: typing.Optional[types.TracebackType],
    ) -> bool: ...
    def __iter__(self) -> HMMFile: ...
    def __next__(self) -> HMM: ...
    def close(self) -> None: ...

class OptimizedProfile(object):
    def __init__(self, m: int, abc: Alphabet) -> None: ...
    def __copy__(self) -> OptimizedProfile: ...
    def is_local(self) -> bool: ...
    def copy(self) -> OptimizedProfile: ...
    @property
    def offsets(self) -> _Offsets: ...
    def write(self, fh_filter: typing.BinaryIO, fh_profile: typing.BinaryIO) -> None: ...

class _Offsets(object):
    def __copy__(self) -> _Offsets: ...
    @property
    def model(self) -> typing.Optional[int]: ...
    @model.setter
    def model(self, model: typing.Optional[int]) -> None: ...
    @property
    def filter(self) -> typing.Optional[int]: ...
    @filter.setter
    def filter(self, filter: typing.Optional[int]) -> None: ...
    @property
    def profile(self) -> typing.Optional[int]: ...
    @profile.setter
    def profile(self, profile: typing.Optional[int]) -> None: ...

class Pipeline(object):
    alphabet: Alphabet
    background: Background
    profile: typing.Optional[Profile]
    def __init__(
        self,
        alphabet: Alphabet,
        background: typing.Optional[Background] = None,
        *,
        bias_filter: bool = True,
        report_e: float = 10.0,
        null2: bool = True,
        seed: typing.Optional[int] = None,
        Z: typing.Optional[float] = None,
        domZ: typing.Optional[float] = None,
    ) -> None: ...
    @property
    def Z(self) -> typing.Optional[float]: ...
    @Z.setter
    def Z(self, Z: typing.Optional[float]) -> None: ...
    @property
    def domZ(self) -> typing.Optional[float]: ...
    @domZ.setter
    def domZ(self, domZ: typing.Optional[float]) -> None: ...
    def clear(self) -> None: ...
    def search(
        self,
        hmm: HMM,
        sequences: typing.Iterable[Sequence],
        hits: typing.Optional[TopHits] = None,
    ) -> TopHits: ...

class Profile(object):
    alphabet: Alphabet
    def __init__(self, M: int, alphabet: Alphabet) -> None: ...
    def __copy__(self) -> Profile: ...
    def clear(self) -> None: ...
    def configure(
        self,
        hmm: HMM,
        background: Background,
        L: int,
        multihit: bool = True,
        local: bool = True,
    ) -> None: ...
    def copy(self) -> Profile: ...
    def is_local(self) -> bool: ...
    def is_multihit(self) -> bool: ...
    def optimized(self) -> OptimizedProfile: ...

class TopHits(typing.Sequence[Hit]):
    Z: float
    domZ: float
    long_targets: bool
    def __init__(self) -> None: ...
    def __bool__(self) -> bool: ...
    def __len__(self) -> int: ...
    @typing.overload
    def __getitem__(self, index: int) -> Hit: ...
    @typing.overload
    def __getitem__(self, index: slice) -> typing.Sequence[Hit]: ...
    def __iadd__(self, other: TopHits) -> TopHits: ...
    def clear(self) -> None: ...
    def sort(self, by: str = "key") -> None: ...
    def is_sorted(self, by: str = "key") -> bool: ...
