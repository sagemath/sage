SAGE_SPKG_CONFIGURE(
    [perl_cpan_polymake_prereq], [
    m4_pushdef([MODULES], m4_include(build/pkgs/perl_cpan_polymake_prereq/distros/cpan.txt))
    AX_PROG_PERL_MODULES(MODULES,
      [],
      [sage_spkg_install_perl_cpan_polymake_prereq=yes
      ]
    )
])
