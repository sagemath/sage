SAGE_SPKG_CONFIGURE([gap], [
  # Default to installing the SPKG, if the check is run at all.
  sage_spkg_install_gap=yes

  m4_pushdef([GAP_MINVER],["4.12.2"])
  m4_pushdef([GAP_LTVER],["4.13.0"])

  SAGE_SPKG_DEPCHECK([ncurses readline zlib], [
    AC_PATH_PROG(GAP, gap)
    AS_IF([test -n "${GAP}"], [
      AC_MSG_CHECKING([for gap version >= GAP_MINVER, < GAP_LTVER])
      dnl GAP will later add the "user" path to the list of root paths
      dnl so long as we don't initialize GAP with -r in Sage. But we
      dnl don't want to include it in the hard-coded list.
      GAPRUN="${GAP} -r -q --bare --nointeract -c"
      _cmd='Display(GAPInfo.KernelInfo.KERNEL_VERSION);'
      GAP_VERSION=$( ${GAPRUN} "${_cmd}" 2>/dev/null )
      AX_COMPARE_VERSION(["${GAP_VERSION}"], [ge], [GAP_MINVER], [dnl
       AC_MSG_RESULT([yes])
       AX_COMPARE_VERSION(["${GAP_VERSION}"], [lt], [GAP_LTVER], [dnl
        AC_MSG_RESULT([yes])
        AC_MSG_CHECKING([for gap root paths])
        _cmd='Display(JoinStringsWithSeparator(GAPInfo.RootPaths,";"));'
        SYS_GAP_ROOT_PATHS=$( ${GAPRUN} "${_cmd}" 2>/dev/null )
        AC_MSG_RESULT([$SYS_GAP_ROOT_PATHS])
        AS_IF([test -n "${SYS_GAP_ROOT_PATHS}"], [
          AC_MSG_CHECKING([for the PrimGrp, SmallGrp, and TransGrp packages])
          # Check for a very minimal set of packages without which the
          # sage test suite will fail. The crazy thing below is a
          # "quadrigraph" for a square bracket.
          _cmd="Display(@<:@"
          _cmd="${_cmd} TestPackageAvailability(\"PrimGrp\"),"
          _cmd="${_cmd} TestPackageAvailability(\"SmallGrp\"),"
          _cmd="${_cmd} TestPackageAvailability(\"TransGrp\")"
          _cmd="${_cmd} @:>@);"
          _output=$( ${GAPRUN} "${_cmd}" 2>/dev/null )
          AS_IF([test $? -eq 0], [
            AS_CASE([$_output],
              [*fail*],[AC_MSG_RESULT([no (at least one package missing)])],[
                # default case, i.e. no "fail"
                AC_MSG_RESULT([yes])

                AC_MSG_CHECKING([if we can link against libgap])
                # That was all for the CLI. Now we check for libgap,
                # too. There's a long list of headers we need in
                # src/sage/libs/gap/gap_includes.pxd, but libgap-api.h
                # combined with the version test above should be
                # sufficient even on systems where the headers are
                # packaged separately.
                _old_libs=$LIBS
                LIBS="${LIBS} -lgap"
                AC_LANG_PUSH([C])
                AC_LINK_IFELSE([
                  AC_LANG_PROGRAM(
                    [[#include <gap/libgap-api.h>]],
                    [[
                      int main(int argc, char** argv) {
                        GAP_Initialize(0, 0, 0, 0, 0);
                        return 0;
                      }
                    ]])
                ],[
                  AC_MSG_RESULT([yes])
                  sage_spkg_install_gap=no
                ],[
                  AC_MSG_RESULT([no])
                ])
                AC_LANG_POP
                LIBS="${_old_libs}"
            ])
          ], [dnl The gap command itself failed
            AC_MSG_RESULT([no (package check command failed)])
          ])
        ])
      ], [dnl Version too new
        AC_MSG_RESULT([no])
      ])
     ], [dnl Version too old
       AC_MSG_RESULT([no])
     ])
    ])
  ])

  m4_popdef([GAP_LTVER])
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
