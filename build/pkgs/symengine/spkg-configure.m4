SAGE_SPKG_CONFIGURE([symengine], [
    m4_pushdef(SAGE_SYMENGINE_VERSION_MAJOR, [0])
    m4_pushdef(SAGE_SYMENGINE_VERSION_MINOR, [11])
    SAGE_SPKG_DEPCHECK([gmp ecm flint mpc mpfr], [
        AC_CHECK_HEADER([symengine/symengine_config.h], [], [sage_spkg_install_symengine=yes])
        AC_MSG_CHECKING([whether we can link a program using symengine])
        SYMENGINE_SAVED_LIBS=$LIBS
        LIBS="$LIBS -lsymengine"
        AC_LINK_IFELSE([
            AC_LANG_PROGRAM([[#include <symengine/expression.h>]],
            [[using SymEngine::Expression;
              Expression x("x");
              auto ex = pow(x+sqrt(Expression(2)), 6);]]
            )], [AC_MSG_RESULT([yes])], [
            AC_MSG_RESULT([no]); sage_spkg_install_symengine=yes
            LIBS=$SYMENGINE_SAVED_LIBS
        ])
        AC_MSG_CHECKING([symengine version >= ]SAGE_SYMENGINE_VERSION_MAJOR[.]SAGE_SYMENGINE_VERSION_MINOR)
        AC_RUN_IFELSE([
            AC_LANG_PROGRAM(
            [[#include <symengine/symengine_config.h>
              #include <stdio.h>
            ]], [[
              fprintf(stderr, "%s\n", SYMENGINE_VERSION);
              if (SYMENGINE_MAJOR_VERSION >]] SAGE_SYMENGINE_VERSION_MAJOR[[) return 0;
              else if (SYMENGINE_MAJOR_VERSION ==]] SAGE_SYMENGINE_VERSION_MAJOR[[ &&
                       SYMENGINE_MINOR_VERSION >=]] SAGE_SYMENGINE_VERSION_MINOR[[) return 0;
              else return 1;
            ]])], [
                AC_MSG_RESULT([yes])
            ], [
                AC_MSG_RESULT([no])
                sage_spkg_install_symengine=yes
            ], [
                dnl assume that the person running cross-compiling
                dnl knows what they are doing
                AC_MSG_RESULT([yes])
            ])
    ])

    m4_popdef([SAGE_SYMENGINE_VERSION_MAJOR])
    m4_popdef([SAGE_SYMENGINE_VERSION_MINOR])
], [], [], []
)

