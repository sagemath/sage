SAGE_SPKG_CONFIGURE([gp2c], [
  # Default to installing the SPKG, if the check is run at all.
  sage_spkg_install_gp2c=yes

  AS_IF([test "x$USING_SYSTEM_PARI" = "xyes"], [
    # We're using the system pari, so we can use the system gp2c
    # if we find a suitable one.
    AC_PATH_PROG([GP2C], [gp2c])
    AS_IF([test -n "$GP2C"], [
      # We found gp2c on the system; use it.
      sage_spkg_install_gp2c=no
    ])
  ])
],[],[
  # Pre-check phase. We need to determine the location of pari.cfg
  # regardless of whether or not we're using the system copy of pari.
  #
  # Store the depcheck result to be reused in this macro's check and
  # post-check phases. We do this in pre-check because the "check"
  # phase itself may be skipped via --with-system-gp2c=no.
  USING_SYSTEM_PARI=no
  SAGE_SPKG_DEPCHECK([pari], [USING_SYSTEM_PARI=yes])
],[
  # Post-check phase. Here we may need to locate pari.cfg if we're using
  # the system's pari (which can put pari.cfg wherever it wants) but sage's
  # gp2c (which needs to know where pari.cfg lives).
  #
  # Can we avoid this if the user hasn't passed --enable-gp2c to ./configure?
  #
  AS_IF([test "x$sage_spkg_install_gp2c" = "xyes"], [
    AS_IF([test "x$USING_SYSTEM_PARI" = "xyes"], [
      # Installing the gp2c package but don't know where pari.cfg is. There's
      # no good way to do this, except to try every known location. If you
      # have two copies of pari.cfg in two locations, this loop will overwrite
      # the location of the earlier with the latter... but how on earth that
      # could happen is beyond my imagination. Since we don't take a
      # --with-pari flag, we guess at the prefix where pari/gp was installed
      # by assuming the executable is in e.g. /usr/bin or /usr/local/bin and
      # then stripping off the "bin" component of that.
      gp_prefix=$(dirname -- "$(dirname -- "$GP")")

      # During the great systemd /bin -> /usr/bin merge, there are
      # systems where /bin is a symlink to /usr/bin. The real "gp"
      # executable lives in /usr/bin, but (in at least one case, trac
      # 31051) we find the copy in /bin first. Why Fedora 32 has the
      # symlink /bin in front of /usr/bin in their PATH is a question
      # I cannot answer; nevertheless, it's pretty easy to add back
      # the "usr" here if we determine that the gp prefix is "/",
      # which it ain't. This is viable only so long as nobody really
      # puts "gp" in /bin, but... they should not.
      AS_IF([test "x$gp_prefix" = "x/"], [gp_prefix=/usr])

      AC_MSG_NOTICE([gp prefix is $gp_prefix])

      # Gentoo:     $gp_prefix/share/pari/pari.cfg
      # Arch/Conda: $gp_prefix/lib/pari/pari.cfg
      # Fedora:     $gp_prefix/share/doc/pari/pari.cfg
      # Fedora from 2.13.2-2: $gp_prefix/lib64/pari/pari.cfg
      m4_foreach([pari_cfg_path], [share/pari,lib64/pari,lib/pari,share/doc/pari], [
        AS_IF([test -f "${gp_prefix}/pari_cfg_path/pari.cfg"], [
          libpari_pari_cfg="${gp_prefix}/pari_cfg_path/pari.cfg"
          AC_MSG_NOTICE([found a pari.cfg at $libpari_pari_cfg])
        ])
      ])

      # Debian:     $gp_prefix/lib/<arch-tuple>/pari/pari.cfg
      #
      # See https://wiki.debian.org/Multiarch/Tuples for a list of valid
      # Debian arch tuples. We rely on "dpkg-architecture" to output the
      # right one. If it doesn't, the "-f" test below prevents anything
      # too bad from happening.
      debian_arch=$(dpkg-architecture -qDEB_BUILD_MULTIARCH 2>/dev/null)
      AS_IF([test -f "${gp_prefix}/lib/${debian_arch}/pari/pari.cfg"], [
        libpari_pari_cfg="${gp_prefix}/lib/${debian_arch}/pari/pari.cfg"
        AC_MSG_NOTICE([found a pari.cfg at $libpari_pari_cfg])
      ])

      # If we can't find pari.cfg, gp2c isn't going to work.
      AS_IF([test -z "$libpari_pari_cfg"], [
        AC_MSG_WARN([using system pari and unable to locate pari.cfg; building the optional pacakge gp2c will not work])
      ])
    ], [
      # Not using the system pari
      libpari_pari_cfg='$SAGE_LOCAL/lib/pari/pari.cfg'
    ])

    AC_MSG_NOTICE([pari.cfg is $libpari_pari_cfg])
  ])
  AC_SUBST(SAGE_PARI_CFG, [$libpari_pari_cfg])
])
