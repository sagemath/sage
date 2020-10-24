*   the nofltk is because giac is compiled without fltk, so some outputs differ
    in the make check

*   In the upstream source, the ifactor function is disabling the pari factorization under macos.
    The reason was interruptions problems when linking guis with the macos binaries provided by upstream 
    (xcas, qcas). http://xcas.e.ujf-grenoble.fr/XCAS/viewtopic.php?f=4&t=1555.
    The macos-ifactor patch enables pari in ifactor under osx because the problem can't be 
    reproduced with the spkg library.
