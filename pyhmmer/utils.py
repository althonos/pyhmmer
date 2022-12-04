import typing
import functools

_T = typing.TypeVar("_T")
_S = typing.TypeVar("_S")


class peekable(typing.Iterator[_T], typing.Generic[_T]):
    """Turn an iterable into a peekable iterable.

    Example:
        >>> it = peekable(range(3))
        >>> it.peek()
        0
        >>> next(it)
        0
        >>> next(it)
        1
        >>> it.peek()
        2

    """

    _sentinel = object()

    def __init__(self, iterable: typing.Iterable[_T]):
        self.it = iter(iterable)
        self.peeked = self._sentinel

    def __iter__(self) -> "peekable[_T]":
        return self

    def __next__(self) -> _T:
        if self.peeked is not self._sentinel:
            self.peeked, item = self._sentinel, self.peeked
        else:
            item = next(self.it)
        return item  # type: ignore

    def peek(self) -> _T:
        if self.peeked is self._sentinel:
            self.peeked = next(self.it)
        return self.peeked  # type: ignore


class singledispatchmethod(typing.Generic[_T]):
    """Single-dispatch generic method descriptor.

    Backported from Python 3.8 to support older Python versions. Additional
    type annotations were obtained from the `functools` stubs in the Python
    `typeshed <https://github.com/python/typeshed>`_ repository.

    """

    def __init__(self, func: typing.Callable[..., _T]) -> None:
        if not callable(func) and not hasattr(func, "__get__"):
            raise TypeError("{!r} is not callable or a descriptor".format(func))
        self.dispatcher = functools.singledispatch(func)
        self.func = func

    @property
    def __isabstractmethod__(self) -> bool:
        return getattr(self.func, "__isabstractmethod__", False)

    def __get__(
        self, obj: _S, cls: typing.Optional[typing.Type[_S]]
    ) -> typing.Callable[..., _T]:
        def _method(*args, **kwargs):  # type: ignore
            method = self.dispatcher.dispatch(args[0].__class__)
            return method.__get__(obj, cls)(*args, **kwargs)

        _method.__isabstractmethod__ = self.__isabstractmethod__  # type: ignore
        _method.register = self.register  # type: ignore
        functools.update_wrapper(_method, self.func)
        return _method

    @typing.overload
    def register(
        self, cls: typing.Type[typing.Any], method: None = ...
    ) -> typing.Callable[[typing.Callable[..., _T]], typing.Callable[..., _T]]:
        ...

    @typing.overload
    def register(
        self, cls: typing.Callable[..., _T], method: None = ...
    ) -> typing.Callable[..., _T]:
        ...

    @typing.overload
    def register(
        self, cls: typing.Type[typing.Any], method: typing.Callable[..., _T]
    ) -> typing.Callable[..., _T]:
        ...

    def register(
        self, cls: typing.Any, method: typing.Optional[typing.Callable[..., _T]] = None
    ) -> typing.Any:
        """Registers a new implementation for the given class."""
        return self.dispatcher.register(cls, func=method)


class SizedIterator(typing.Generic[_T], typing.Iterator[_T], typing.Sized):
    """An iterator with a known number of items.
    """

    def __init__(self, n: int, iterable: typing.Iterable[_T]) -> None:
        self._it = iter(iterable)
        self._n = n

    def __len__(self) -> int:
        return self._n

    def __iter__(self: _S) -> _S:
        return self

    def __next__(self) -> _T:
        if self._n > 0:
            self._n -= 1
        return next(self._it)