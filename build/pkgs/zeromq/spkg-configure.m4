SAGE_SPKG_CONFIGURE([zeromq], [
    AX_ZMQ([4.2.5], [], [sage_spkg_install_zeromq=yes])
])
