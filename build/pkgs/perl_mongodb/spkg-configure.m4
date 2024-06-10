SAGE_SPKG_CONFIGURE(
    [perl_mongodb], [
    m4_pushdef([MODULES], m4_include(build/pkgs/perl_mongodb/distros/cpan.txt))
    AX_PROG_PERL_MODULES(MODULES,
      [],
      [sage_spkg_install_perl_mongodb=yes
      ]
    )
])
