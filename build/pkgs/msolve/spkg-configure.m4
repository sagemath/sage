SAGE_SPKG_CONFIGURE([msolve], [
    PKG_CHECK_MODULES([msolve], [msolve >= 0.6.5], [], [
       AC_CACHE_CHECK([for msolve], [ac_cv_path_MSOLVE],
         [AC_PATH_PROGS_FEATURE_CHECK([MSOLVE], [msolve],
           [msolvin=$(mktemp)
            msolvout=$(mktemp)
            msolvchk=$(mktemp)
            echo -e 'x,y\n0\nx-y,\nx*y-1' >$msolvin
	    changequote(<<, >>)dnl
            echo -e '[0, [1,\n[[[-1, -1], [-1, -1]], [[1, 1], [1, 1]]]\n]]:' >$msolvchk
	    changequote([, ])dnl
            $ac_path_MSOLVE -f $msolvin -o $msolvout
            AS_IF([test x$(diff $msolvout $msolvchk) = x],
                  [ac_cv_path_MSOLVE=$ac_path_MSOLVE
                   ac_path_MSOLVE_found=:],
                  [sage_spkg_install_msolve=yes
                   AC_MSG_RESULT([could not find working msolve])
            ])
            ], [sage_spkg_install_msolve=yes
                AC_MSG_RESULT([could not find working msolve])
	    ])])
        ])
])
