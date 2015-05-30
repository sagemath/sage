dnl provided on autoconf mailing list by Russ Albery

AC_DEFUN([AX_PROG_PERL_VERSION],
[AC_CACHE_CHECK([for Perl version $1 or later], [ax_cv_prog_perl_version],
    [AS_IF(["$PERL" -e 'require $1;' >/dev/null 2>&1],
        [ax_cv_prog_perl_version=yes],
        [ax_cv_prog_perl_version=no])])
 AS_IF([test x"$ax_cv_prog_perl_version" = xyes], [$2], [$3])])
