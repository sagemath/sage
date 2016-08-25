"_separate_channels"

def _separate_channels(_data, _width, _nchannels):
    """
    Separates the channels. This is an internal helper method for
    precomputing some data while initializing the wav object.
    """
    cdef int nchannels = _nchannels
    cdef int width = _width
    cdef char* data = _data
    cdef int n
    cdef short int x

    cdef int l = len(_data) / (width)

    channel_data = [[] for i in xrange(nchannels)]
    if width == 1:
        # handle the one byte case

        for n from 0 <= n < l:
            channel_data[n % nchannels].append(ord(data[n])-127)

    elif width == 2:
        a = 32768
        for n from 0 <= n < l:
            # compute the value as an integer
            x = <int> (data[2*n]) + 256 * <int>(data[2*n + 1])
            #x -= 65536*(x > a)
            channel_data[n % nchannels].append(x)
    else:
        raise NotImplementedError("greater than 16-bit wavs not supported")

    return channel_data
