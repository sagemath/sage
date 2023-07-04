# We don't use the "converseg" program, but we need to ensure that we
# only detect nauty >= 2.6 because we use the digraph6 format from
# that version -- and converseg was added in nauty-2.6.
AC_DEFUN([SAGE_TEST_NAUTY_PROGS], [
    m4_foreach([nautyprog], [directg, gentourng, geng, genbg, gentreeg, converseg], [
      AC_PATH_PROG([$2]nautyprog, [[$1]nautyprog])
      AS_IF([test x$[$2]nautyprog = x], [sage_spkg_install_nauty=yes])
    ])
    AC_SUBST(SAGE_NAUTY_BINS_PREFIX, ['$1'])
])
SAGE_SPKG_CONFIGURE([nauty], [
  AC_PATH_PROG([GENGCHECK],[geng])
  AS_IF([test x$GENGCHECK = x], [
    AC_PATH_PROG([GENGnautyCHECK],[nauty-geng])
     AS_IF([test x$GENGnautyCHECK = x], [sage_spkg_install_nauty=yes],
       [SAGE_TEST_NAUTY_PROGS(nauty-,nau)])
    ], [SAGE_TEST_NAUTY_PROGS(,foo)])
 ], [], [], [
 AS_IF([test x$sage_spkg_install_nauty = xyes], [
   AC_SUBST(SAGE_NAUTY_BINS_PREFIX, [''])
 ])
])
