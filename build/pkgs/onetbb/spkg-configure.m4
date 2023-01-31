SAGE_SPKG_CONFIGURE([onetbb], [dnl
    AC_MSG_CHECKING([whether oneTBB >= 2018 is available])
    rm -rf conftest_srcdir
    mkdir conftest_srcdir
    cat > conftest_srcdir/CMakeLists.txt <<EOF
cmake_minimum_required (VERSION 3.11.0)
project(dummy)
# from https://github.com/scipopt/papilo/blob/master/CMakeLists.txt
find_package(TBB 2018 COMPONENTS tbb tbbmalloc REQUIRED)
EOF
    AS_IF([cmake -S conftest_srcdir -B conftest_srcdir/build >& ]AS_MESSAGE_LOG_FD[ 2>&1], [dnl
        AC_MSG_RESULT([yes])
    ], [dnl
        AC_MSG_RESULT([no])
        AS_VAR_SET([sage_spkg_install_onetbb], [yes])
    ])
])
