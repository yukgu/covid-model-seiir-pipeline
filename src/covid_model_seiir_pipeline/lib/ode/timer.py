from numba import njit, objmode
from numba.core import types
from numba.typed import Dict
import numpy as np
import ctypes


CLOCK_MONOTONIC = 0x1
clock_gettime_proto = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_int,
                                       ctypes.POINTER(ctypes.c_long))
pybind = ctypes.CDLL(None)
clock_gettime_addr = pybind.clock_gettime
clock_gettime_fn_ptr = clock_gettime_proto(clock_gettime_addr)


@njit
def timenow():
    timespec = np.zeros(2, dtype=np.int64)
    clock_gettime_fn_ptr(CLOCK_MONOTONIC, timespec.ctypes)
    ts = timespec[0]
    tns = timespec[1]
    return np.float64(ts) + 1e-9 * np.float64(tns)


START = Dict.empty(
    key_type=types.unicode_type,
    value_type=types.float64,
)

RESULTS = Dict.empty(
    key_type=types.unicode_type,
    value_type=types.float64,
)


@njit
def log_start(name):
    START[name] = timenow()


@njit
def log_end(name):
    assert name in START
    if name in RESULTS:
        RESULTS[name] += timenow() - START[name]
    else:
        RESULTS[name] = timenow() - START[name]

@njit
def clear_results():
    for name in RESULTS:
        del RESULTS[name]
