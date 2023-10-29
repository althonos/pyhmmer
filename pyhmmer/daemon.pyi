import socket
import typing
import types

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal  # type: ignore

import pyhmmer.plan7
from pyhmmer.easel import Sequence, MSA
from pyhmmer.plan7 import TopHits, HMM, Builder

BIT_CUTOFFS = Literal["gathering", "trusted", "noise"]

class Client:
    address: str
    port: int
    socket: socket.socket
    def __init__(self, address: str = "127.0.0.1", port: int = 51371) -> None: ...
    def __enter__(self) -> Client: ...
    def __exit__(
        self,
        exc_type: typing.Optional[typing.Type[BaseException]],
        exc_value: typing.Optional[BaseException],
        traceback: typing.Optional[types.TracebackType],
    ) -> bool: ...
    def __repr__(self) -> str: ...
    def connect(self) -> None: ...
    def close(self) -> None: ...
    def search_seq(
        self,
        query: Sequence,
        db: int = 1,
        ranges: typing.Optional[typing.List[typing.Tuple[int, int]]] = None,
        *,
        bias_filter: bool = True,
        null2: bool = True,
        seed: int = 42,
        Z: typing.Optional[float] = None,
        domZ: typing.Optional[float] = None,
        F1: float = 0.02,
        F2: float = 1e-3,
        F3: float = 1e-5,
        E: float = 10.0,
        T: typing.Optional[float] = None,
        domE: float = 10.0,
        domT: typing.Optional[float] = None,
        incE: float = 0.01,
        incT: typing.Optional[float] = None,
        incdomE: float = 0.01,
        incdomT: typing.Optional[float] = None,
        bit_cutoffs: typing.Optional[BIT_CUTOFFS] = None,
    ) -> TopHits: ...
    def search_hmm(
        self,
        query: HMM,
        db: int = 1,
        ranges: typing.Optional[typing.List[typing.Tuple[int, int]]] = None,
        *,
        bias_filter: bool = True,
        null2: bool = True,
        seed: int = 42,
        Z: typing.Optional[float] = None,
        domZ: typing.Optional[float] = None,
        F1: float = 0.02,
        F2: float = 1e-3,
        F3: float = 1e-5,
        E: float = 10.0,
        T: typing.Optional[float] = None,
        domE: float = 10.0,
        domT: typing.Optional[float] = None,
        incE: float = 0.01,
        incT: typing.Optional[float] = None,
        incdomE: float = 0.01,
        incdomT: typing.Optional[float] = None,
        bit_cutoffs: typing.Optional[BIT_CUTOFFS] = None,
    ) -> TopHits: ...
    def scan_seq(
        self,
        query: Sequence,
        db: int = 1,
        *,
        bias_filter: bool = True,
        null2: bool = True,
        seed: int = 42,
        Z: typing.Optional[float] = None,
        domZ: typing.Optional[float] = None,
        F1: float = 0.02,
        F2: float = 1e-3,
        F3: float = 1e-5,
        E: float = 10.0,
        T: typing.Optional[float] = None,
        domE: float = 10.0,
        domT: typing.Optional[float] = None,
        incE: float = 0.01,
        incT: typing.Optional[float] = None,
        incdomE: float = 0.01,
        incdomT: typing.Optional[float] = None,
        bit_cutoffs: typing.Optional[BIT_CUTOFFS] = None,
    ) -> TopHits: ...
    def iterate_seq(
        self,
        query: Sequence,
        db: int = 1,
        ranges: typing.Optional[typing.List[typing.Tuple[int, int]]] = None,
        builder: typing.Optional[Builder] = None,
        select_hits: typing.Optional[typing.Callable[[TopHits], None]] = None,
        *,
        bias_filter: bool = True,
        null2: bool = True,
        seed: int = 42,
        Z: typing.Optional[float] = None,
        domZ: typing.Optional[float] = None,
        F1: float = 0.02,
        F2: float = 1e-3,
        F3: float = 1e-5,
        E: float = 10.0,
        T: typing.Optional[float] = None,
        domE: float = 10.0,
        domT: typing.Optional[float] = None,
        incE: float = 0.001,
        incT: typing.Optional[float] = None,
        incdomE: float = 0.001,
        incdomT: typing.Optional[float] = None,
        bit_cutoffs: typing.Optional[BIT_CUTOFFS] = None,
    ) -> IterativeSearch: ...
    def iterate_hmm(
        self,
        query: HMM,
        db: int = 1,
        ranges: typing.Optional[typing.List[typing.Tuple[int, int]]] = None,
        builder: typing.Optional[Builder] = None,
        select_hits: typing.Optional[typing.Callable[[TopHits], None]] = None,
        *,
        bias_filter: bool = True,
        null2: bool = True,
        seed: int = 42,
        Z: typing.Optional[float] = None,
        domZ: typing.Optional[float] = None,
        F1: float = 0.02,
        F2: float = 1e-3,
        F3: float = 1e-5,
        E: float = 10.0,
        T: typing.Optional[float] = None,
        domE: float = 10.0,
        domT: typing.Optional[float] = None,
        incE: float = 0.001,
        incT: typing.Optional[float] = None,
        incdomE: float = 0.001,
        incdomT: typing.Optional[float] = None,
        bit_cutoffs: typing.Optional[BIT_CUTOFFS] = None,
    ) -> IterativeSearch: ...

class IterativeSearch(pyhmmer.plan7.IterativeSearch):
    client: Client
    db: int
    def __init__(
        self,
        client: Client,
        query: typing.Union[HMM, Sequence],
        db: int,
        builder: Builder,
        ranges: typing.Optional[typing.List[typing.Tuple[int, int]]] = None,
        select_hits: typing.Optional[typing.Callable[[TopHits], None]] = None,
        options: typing.Optional[typing.Dict[str, object]] = None,
    ) -> None: ...
