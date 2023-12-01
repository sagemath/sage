SAGE_SPKG_CONFIGURE([gap], [
  # Default to installing the SPKG, if the check is run at all.
  sage_spkg_install_gap=yes

  m4_pushdef([GAP_MINVER],["4.12.2"])

  SAGE_SPKG_DEPCHECK([ncurses readline zlib], [
    AC_PATH_PROG(GAP, gap)
    AS_IF([test -n "${GAP}"], [
      AC_MSG_CHECKING([for gap version GAP_MINVER or newer])

      # GAP will later add the "user" path to the list of root paths
      # so long as we don't initialize GAP with -r in Sage. But we
      # don't want to include it in the hard-coded list.
      GAPC="${GAP} -r -q --bare --nointeract -c"
      _cmd='Display(GAPInfo.KernelInfo.KERNEL_VERSION);'
      GAP_VERSION=$(${GAPC} "${_cmd}")
      AX_COMPARE_VERSION(["${GAP_VERSION}"], [ge], [GAP_MINVER], [
        AC_MSG_RESULT([yes])
        AC_MSG_CHECKING([the gap root paths])
        _cmd='Display(JoinStringsWithSeparator(GAPInfo.RootPaths,";"));'
        SYS_GAP_ROOT_PATHS=$(${GAPC} "${_cmd}")
        AC_MSG_RESULT([$SYS_GAP_ROOT_PATHS])
        AS_IF([test -n "${SYS_GAP_ROOT_PATHS}"], [
          AC_MSG_CHECKING([for the PrimGrp, SmallGrp, and TransGrp packages])
          # The crazy thing below is a "quadrigraph" for a square bracket
          _cmd="Display(@<:@"
          _cmd="${_cmd} TestPackageAvailability(\"PrimGrp\"),"
          _cmd="${_cmd} TestPackageAvailability(\"SmallGrp\"),"
          _cmd="${_cmd} TestPackageAvailability(\"TransGrp\")"
          _cmd="${_cmd} @:>@);"
          _output=$( ${GAPC} "${_cmd}" )
          AS_IF([test $? -eq 0], [
            AS_CASE([$_output],
              [*fail*],[AC_MSG_RESULT([no (at least one package missing)])],[
                # default case, i.e. no "fail"
                AC_MSG_RESULT([yes])
                sage_spkg_install_gap=no
            ])
          ], [
            # The gap command itself failed
            AC_MSG_RESULT([no (package check command failed)])
          ])
        ])
      ],[
        # Version too old
        AC_MSG_RESULT([no])
      ])
    ])
  ])

  m4_popdef([GAP_MINVER])
],[],[],[
  # This is the post-check phase, where we make sage-conf
  # substitutions, in this case of GAP_ROOT_PATHS. We begin with the
  # two root paths used by the sage distribution. The '${prefix}' is
  # a magic string that sage-conf will replace.
  GAP_ROOT_PATHS='${prefix}/lib/gap;${prefix}/share/gap';

  AS_IF([test "${sage_spkg_install_gap}" = "no"],[
    # If we're using the system GAP, append the system root
    # paths to the existing two sage paths.
    GAP_ROOT_PATHS="${GAP_ROOT_PATHS};${SYS_GAP_ROOT_PATHS}"
  ])

  AC_SUBST(GAP_ROOT_PATHS, "${GAP_ROOT_PATHS}")
])
