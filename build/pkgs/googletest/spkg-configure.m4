SAGE_SPKG_CONFIGURE([googletest], [dnl
    AC_MSG_CHECKING([whether googletest is available])
    rm -rf conftest_srcdir
    mkdir conftest_srcdir
    cat > conftest_srcdir/CMakeLists.txt <<EOF
cmake_minimum_required (VERSION 3.11.0)
project(dummy)
find_package(GTest REQUIRED)
EOF
    AS_IF([cmake -S conftest_srcdir -B conftest_srcdir/build >& ]AS_MESSAGE_LOG_FD[ 2>&1], [dnl
        AC_MSG_RESULT([yes])
    ], [dnl
        AC_MSG_RESULT([no])
        AS_VAR_SET([sage_spkg_install_googletest], [yes])
    ])
])
