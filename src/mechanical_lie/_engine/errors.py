"""Explicit invariant failures, including in optimized Python runs."""

from functools import lru_cache, wraps


class ValidationError(ValueError):
    """A mathematical or serialized-data invariant did not hold."""


def cached_mapping(function):
    """Cache immutable item tuples; return a fresh dict to every caller."""

    @lru_cache(maxsize=4096)
    def cached(*args):
        return tuple(function(*args).items())

    @wraps(function)
    def wrapped(*args):
        return dict(cached(*args))

    wrapped.cache_clear = cached.cache_clear
    return wrapped
