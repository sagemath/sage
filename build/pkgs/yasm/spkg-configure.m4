SAGE_SPKG_CONFIGURE(
    [yasm],
    # Yasm is only needed on x86(_64) systems; check also for system yasm which
    # must support "adox" (new Skylake instruction)
    [AC_MSG_CHECKING([for yasm supporting the adox instruction])
     AC_PATH_PROGS_FEATURE_CHECK([YASM], [yasm],
        [[{ echo "BITS 64"; echo "adox rax, rax"; } | ${ac_path_YASM} - -o /dev/null >/dev/null 2>/dev/null && ac_cv_path_YASM=${ac_path_YASM}]],
        [sage_spkg_install_yasm=yes; ac_cv_path_YASM=no])
     AC_MSG_RESULT($ac_cv_path_YASM)],
    [dnl REQUIRED-CHECK
     AS_CASE("$host_cpu", [i@<:@0-9@:>@86|x86_64|amd64], [], [sage_require_yasm=no])
     AC_REQUIRE([SAGE_SPKG_CONFIGURE_MPIR])
     AS_IF([test x$sage_spkg_install_mpir = xno], [
        sage_require_yasm=no
     ])
    ]
)
