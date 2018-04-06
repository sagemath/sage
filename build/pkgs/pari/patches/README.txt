======================================================================
Current patches to PARI in Sage:
======================================================================

Patches to configuration files:
* get_dlcflags.patch (Jeroen Demeyer): Add -fno-common to DLCFLAGS on
  Darwin. Submitted upstream, but upstream only applied it for PowerPC.
  Since this doesn't break anything and only improves performance, add
  the flag unconditionally.

C files:
* stackwarn.patch (Jeroen Demeyer, #19883): do not display warnings
  regarding the stack size (unless DEBUGMEM is set).
