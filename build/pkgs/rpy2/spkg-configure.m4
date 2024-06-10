SAGE_SPKG_CONFIGURE([rpy2], [
    SAGE_PYTHON_PACKAGE_CHECK([rpy2])
], [dnl REQUIRED-CHECK
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_R])
    dnl rpy2 is only needed when there is a usable system R
    AS_VAR_IF([sage_spkg_install_r], [yes], [dnl
        AS_VAR_IF([sage_use_system_r], [installed], [dnl
            dnl Legacy SPKG installation of r
            AS_VAR_SET([SPKG_REQUIRE], [yes])
        ], [dnl No system package, no legacy SPKG installation
            AS_VAR_SET([SPKG_REQUIRE], [no])
        ])
    ])
])
