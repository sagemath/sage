SAGE_SPKG_CONFIGURE([sympow], [
   AC_PATH_PROG([SYMPOW], [sympow])
   AS_IF([test -z "$ac_cv_path_SYMPOW"], [sage_spkg_install_sympow=yes
      ], [
      AC_MSG_CHECKING([whether sympow works well (cf. :trac:30147 and :trac:30633)])
      sympow_rank_test1=`echo "@<:@1,-1,0,-79,289@:>@" | sympow -analrank | grep ^"Analytic Rank is 4"`
      sympow_rank_test2=`sympow -curve @<:@1,-1,0,-79,289@:>@ -analrank | grep ^"Analytic Rank is 4"`
      sympow_rank_test3=`sympow -curve @<:@0,1,1,-2,0@:>@ -analrank | grep ^"Analytic Rank is 2"`
      AS_IF([test x"$sympow_rank_test1" = x -o x"$sympow_rank_test2" = x -o x"$sympow_rank_test3" = x], [
          AC_MSG_RESULT([no; cannot use system sympow])
          sage_spkg_install_sympow=yes
          ], [
          AC_MSG_RESULT([yes; use system sympow])
      ])
   ])
])
