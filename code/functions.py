# -*- coding: utf-8 -*-

from numba import jit

# Source: http://stackoverflow.com/a/29799815/1501575
# Pre-compiled function to find first element of array greater than
@jit(nopython=True)
def find_first_item_greater_than(vec, value):
    for i in xrange(len(vec)):
        if vec[i] > value:
            return i
    return -1
