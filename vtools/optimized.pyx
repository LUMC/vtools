"""
vtools.optimized
~~~~~~~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""
cimport numpy as np

cpdef int amount_atleast(np.int64_t[::1] values, int atleast):
    """
    Return amount of values at least `atleast`
    :param values: Iterable of int
    :param atleast: int
    :return: int
    """
    cdef int passed = 0
    cdef int val
    cdef size_t i
    for i in range(values.shape[0]):
        val = values[i]
        if val >= atleast:
            passed += 1
    return passed
