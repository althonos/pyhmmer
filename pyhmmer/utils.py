import typing

_Item = typing.TypeVar("_Item")

class peekable(typing.Iterator[_Item], typing.Generic[_Item]):
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

    def __init__(self, iterable: typing.Iterable[_Item]):
        self.it = iter(iterable)
        self.peeked = self._sentinel

    def __iter__(self) -> "peekable[_Item]":
        return self

    def __next__(self) -> _Item:
        if self.peeked is not self._sentinel:
            self.peeked, item = self._sentinel, self.peeked
        else:
            item = next(self.it)
        return item  # type: ignore

    def peek(self) -> _Item:
        if self.peeked is self._sentinel:
            self.peeked = next(self.it)
        return self.peeked  # type: ignore
