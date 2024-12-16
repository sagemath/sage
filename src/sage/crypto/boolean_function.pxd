cdef inline unsigned int hamming_weight(unsigned int x) noexcept:
    # valid for 32bits
    x -= (x>>1) & 0x55555555UL                        # 0-2 in 2 bits
    x = ((x>>2) & 0x33333333UL) + (x & 0x33333333UL)  # 0-4 in 4 bits
    x = ((x>>4) + x) & 0x0f0f0f0fUL                   # 0-8 in 8 bits
    x *= 0x01010101UL
    return x>>24

cdef walsh_hadamard(long *f, int ldn)
