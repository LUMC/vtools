"""
vtools.optimized
~~~~~~~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""


cpdef int amount_atleast(values, int atleast):
    """
    Return amount of values at least `atleast`
    :param values: Iterable of int
    :param atleast: int
    :return: int
    """
    cdef int passed = 0
    cdef int val
    for val in values:
        if val >= atleast:
            passed += 1
    return passed
