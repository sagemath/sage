======================================================================
Current patches to PARI in Sage:
======================================================================

Patches to configuration files:
* get_config_options.patch (Leif Leonhardy): Catch invalid arguments to
  "--graphic" (and treat such as an error) since otherwise strange
  compilation errors might occur (cf. #9722, too).
* get_dlcflags.patch (Jeroen Demeyer): Add -fno-common to DLCFLAGS on
  Darwin. Submitted upstream, but upstream only applied it for PowerPC.
  Since this doesn't break anything and only improves performance, add
  the flag unconditionally.

C files:
* stackwarn.patch (Jeroen Demeyer, #19883): do not display warnings
  regarding the stack size (unless DEBUGMEM is set).

Building documentation
* perl_inc.patch (Francois Bissey, #21622): import upstream patch
to deal with perl install where "." is not part of @INC anymore
in perl.
