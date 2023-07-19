# -*- python -*-

cdef CMR *cmr = NULL


cdef CMR_CALL(CMR_ERROR _cmr_error):
    if _cmr_error == CMR_OKAY:
        return
    if _cmr_error == CMR_ERROR_INPUT:
        raise RuntimeError("User input error")
    if _cmr_error == CMR_ERROR_MEMORY:
        raise RuntimeError("Memory (re)allocation failed")
    if _cmr_error == CMR_ERROR_INVALID:
        raise RuntimeError("Invalid input")
    if _cmr_error == CMR_ERROR_TIMEOUT:
        raise RuntimeError("Time limit exceeded")
    raise RuntimeError("Unknown error")
