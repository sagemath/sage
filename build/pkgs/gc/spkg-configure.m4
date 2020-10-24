SAGE_SPKG_CONFIGURE([gc], [
    SAGE_SPKG_DEPCHECK([libatomic_ops], [
        dnl #30629: Reject system libgc that breaks ECL build on WSL
        AC_MSG_CHECKING([whether we run on WSL])
        AS_IF([grep -q -E -i "wsl|microsoft" /proc/sys/kernel/osrelease 2>/dev/null], [
          AC_MSG_RESULT([yes; disabling using of system gc])
          sage_spkg_install_gc=yes
        ], [
          AC_MSG_RESULT([no])
          dnl  checking with pkg-config
          PKG_CHECK_MODULES([GC], [bdw-gc-threaded >= 7.6.4], [], [
            PKG_CHECK_MODULES([GC], [bdw-gc >= 7.6.4], [], [
              sage_spkg_install_gc=yes])])
        ])
    ])
])
