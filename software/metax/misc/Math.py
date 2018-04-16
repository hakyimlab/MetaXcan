import numpy
from numpy.core import (product, asarray, dot, transpose, multiply, newaxis, maximum)

def _rc(s, tolerance):
    cutoff = tolerance * maximum.reduce(s)
    return cutoff

def _ac(s, tolerance):
    return tolerance

def crpinv(a, rcond=1e-15, epsilon=None):
    """Pseudo inverse, relative cutoff """
    return _inv(a, _rc, rcond, epsilon)

def capinv(a, rcond=1e-15, epsilon=None):
    """Pseudo inverse, absolute cutoff"""
    return _inv(a, _ac, rcond, epsilon)

def _inv(a, cf, rcond, epsilon):
    """
    modified pseudo inverse
    """

    def _assertNoEmpty2d(*arrays):
        for a in arrays:
            if a.size == 0 and product(a.shape[-2:]) == 0:
                raise RuntimeError("Arrays cannot be empty")

    def _makearray(a):
        new = asarray(a)
        wrap = getattr(a, "__array_prepare__", new.__array_wrap__)
        return new, wrap

    a, wrap = _makearray(a)
    _assertNoEmpty2d(a)

    if epsilon is not None:
        epsilon = numpy.repeat(epsilon, a.shape[0])
        epsilon = numpy.diag(epsilon)
        a = a + epsilon
    a = a.conjugate()
    #WARNING! the "s" eigenvalues might not equal the eigenvalues of eigh
    u, s, vt = numpy.linalg.svd(a, 0)
    m = u.shape[0]
    n = vt.shape[1]
    eigen = numpy.copy(s)

    # cutoff = rcond*maximum.reduce(s)
    cutoff = cf(s, rcond)
    for i in range(min(n, m)):
        # The first Singular Value will always be selected because we want at least one, and the first is the highest
        if s[i] >= cutoff or i==0:
            s[i] = 1. / s[i]
        else:
            s[i] = 0.

    n_indep = numpy.count_nonzero(s)
    res = dot(transpose(vt), multiply(s[:, newaxis], transpose(u)))
    return wrap(res), n_indep, eigen

def standardize(x):
    mean = numpy.mean(x)
    #follow R's convention, ddof=1
    scale = numpy.std(x, ddof=1)
    if scale == 0:
        return None
    x = x - mean
    x = x / scale
    return x