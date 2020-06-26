SAGE_SPKG_CONFIGURE([gp2c], [
  # Default to installing the SPKG
  sage_spkg_install_gp2c=yes

  # And to using sage's pari.cfg (from the pari SPKG).
  libpari_pari_cfg='$SAGE_LOCAL/lib/pari/pari.cfg'

  SAGE_SPKG_DEPCHECK([pari], [
    # We're using the system pari, so we can use the system gp2c
    # if we find a suitable one.
    AC_PATH_PROG([GP2C], [gp2c])
    if test -n "$GP2C"; then
      # We found gp2c on the system; use it.
      sage_spkg_install_gp2c=no
    fi

    # Either way, if we're using the system pari, we no longer want
    # to use sage's pari.cfg.
    libpari_pari_cfg=''
  ])
],[],[],[
  # Post-check phase. Here we may need to locate pari.cfg if we're using
  # the system's pari (which can put pari.cfg wherever it wants) but sage's
  # gp2c (which needs to know where pari.cfg lives).
  #
  # Can we avoid this if the user hasn't passed --enable-gp2c to ./configure?
  #
  if test "x$sage_spkg_install_gp2c" = "xyes"; then
    if test -z "$libpari_pari_cfg"; then
      # Installing the gp2c package but don't know where pari.cfg is. There's
      # no good way to do this, except to try every known location. If you
      # have two copies of pari.cfg in two locations, this loop will overwrite
      # the location of the earlier with the latter... but how on earth that
      # could happen is beyond my imagination. Since we don't take a
      # --with-pari flag, we guess at the prefix where pari/gp was installed
      # by assuming the executable is in e.g. /usr/bin or /usr/local/bin and
      # then stripping off the "bin" component of that.
      gp_prefix=$(dirname -- "$(dirname -- "$GP")")
      AC_MSG_NOTICE([gp prefix is $gp_prefix])

      # Gentoo:     $gp_prefix/share/pari/pari.cfg
      # Arch/Conda: $gp_prefix/lib/pari/pari.cfg
      # Fedora:     $gp_prefix/share/doc/pari/pari.cfg
      m4_foreach([pari_cfg_path], [share/pari,lib/pari,share/doc/pari], [
        if test -f "${gp_prefix}/pari_cfg_path/pari.cfg"; then
          libpari_pari_cfg="${gp_prefix}/pari_cfg_path/pari.cfg"
          AC_MSG_NOTICE([found a pari.cfg at $libpari_pari_cfg])
        fi
      ])

      # Debian:     $gp_prefix/lib/<arch-tuple>/pari/pari.cfg
      #
      # See https://wiki.debian.org/Multiarch/Tuples for a list of valid
      # Debian arch tuples. We rely on "dpkg-architecture" to output the
      # right one. If it doesn't, the "-f" test below prevents anything
      # too bad from happening.
      debian_arch=$(dpkg-architecture -qDEB_BUILD_MULTIARCH 2>/dev/null)
      if test -f "${gp_prefix}/lib/${debian_arch}/pari/pari.cfg"; then
        libpari_pari_cfg="${gp_prefix}/lib/${debian_arch}/pari/pari.cfg"
        AC_MSG_NOTICE([found a pari.cfg at $libpari_pari_cfg])
      fi

      # If we can't find pari.cfg, gp2c isn't going to work.
      if test -z "$libpari_pari_cfg"; then
        AC_MSG_ERROR([unable to locate pari.cfg])
      fi
    fi

    # Only print out the location of pari.cfg if we're not using the
    # system gp2c, because if we're using the system gp2c, who cares.
    AC_MSG_NOTICE([pari.cfg is $libpari_pari_cfg])
  fi
  AC_SUBST(SAGE_PARI_CFG, [$libpari_pari_cfg])
])
