"""
Contains L0 to L6q from Rice-McMurdie
"""

from numpy import dot, zeros, power
from params import *


def L0(data):
    result = zeros(buffer_size - 255)
    for m in range(buffer_size - 255):
        result[m] = abs(dot(data[m : m + 256], preamble_template))
    return result

# A version of L0 that can correlate smaller sections of data
def L0_temp(data):
    result = zeros(data.size)
    for m in range(data.size - 255):
        result[m] = abs(dot(data[m : m + 256], preamble_template))
    return result

# A version of L0 that only returns one correlated sample at a time
def L0_single(data):
    result = abs(dot(data[0 : 256], preamble_template))
    return result

def L6_single(data):
    result = power(abs(dot(data[0 : 128], preamble_template[0 : 128])),2) + power(abs(dot(data[128 : 256], preamble_template[128 : 256])),2)
    return result
