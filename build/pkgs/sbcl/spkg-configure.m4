SAGE_SPKG_CONFIGURE([sbcl], [dnl
  m4_pushdef([SAGE_SBCL_MINVER], ["1.4.16"])
   AC_CACHE_CHECK([for sbcl >= SAGE_SBCL_MINVER], [ac_cv_path_SBCL], [
       AC_PATH_PROGS_FEATURE_CHECK([SBCL], [sbcl], [
            sbcl_version=`$ac_path_SBCL --version 2>&1 \
                | $SED -n -e 's/SBCL *\([[0-9]]*\.[[0-9]]*\.[[0-9]]*\)/\1/p'`
            AS_IF([test -n "$sbcl_version"], [
                AX_COMPARE_VERSION([$sbcl_version], [ge], [SAGE_SBCL_MINVER], [
                    ac_cv_path_SBCL="$ac_path_SBCL"
                    ac_path_SBCL_found=:
                ])
            ])
        ])
    ])
  AS_IF([test -z "$ac_cv_path_SBCL" ], [
        sage_spkg_install_sbcl=yes
        AC_SUBST(SAGE_FRICAS_LISP, ecl)], [
        AC_SUBST(SAGE_FRICAS_LISP, "sbcl --dynamic-space-size 4Gb")])
  m4_popdef([SAGE_SBCL_MINVER])
])
