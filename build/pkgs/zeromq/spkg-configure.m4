SAGE_SPKG_CONFIGURE([zeromq], [
    dnl Issue #31624: Avoid C++ ABI issues
    SAGE_SPKG_DEPCHECK([gcc], [
        AX_ZMQ([4.2.5], [], [sage_spkg_install_zeromq=yes])
    ])
])
