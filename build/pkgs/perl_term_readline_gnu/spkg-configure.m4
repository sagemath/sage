SAGE_SPKG_CONFIGURE(
    [perl_term_readline_gnu], [dnl direct testing for Term::ReadLine::Gnu does not work
    AX_PROG_PERL_MODULES( Term::ReadLine, [dnl check that it's a GNU one
      AC_MSG_CHECKING( Term::ReadLine module...)
      # #29563 TERM needs to be set to a value, or Term::ReadLine::Gnu may refuse to load
      TERM=vt100 $PERL -e "use Term::ReadLine; Term::ReadLine->get_all_function_names" > /dev/null 2>&1
      if test $? -ne 0; then
        AC_MSG_RESULT(non-GNU)
        sage_spkg_install_perl_term_readline_gnu=yes
      else
        AC_MSG_RESULT(GNU)
      fi
     ],[sage_spkg_install_perl_term_readline_gnu=yes])
])
