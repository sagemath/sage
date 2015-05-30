======================================================================
Current patches to PARI in Sage:
======================================================================

Patches to configuration files:
* get_ld.patch (Carl Witty): Disable -rpath.
* get_config_options.patch (Leif Leonhardy): Catch invalid arguments to
  "--graphic" (and treat such as an error) since otherwise strange
  compilation errors might occur (cf. #9722, too).
* get_dlcflags.patch (Jeroen Demeyer): Add -fno-common to DLCFLAGS on
  Darwin. Submitted upstream, but upstream only applied it for PowerPC.
  Since this doesn't break anything and only improves performance, add
  the flag unconditionally.
* KERNELCFLAGS.patch (Jeroen Demeyer): when SAGE_DEBUG=yes, compile
  kernel files with -O1 instead of -funroll-loops; -O0 gives a
  segmentation fault on some OS X systems when doing
  factor(10356613*10694706299664611221)
  See #13921, also reported upstream:
  - http://pari.math.u-bordeaux.fr/archives/pari-dev-1301/msg00000.html

C files:
* det_garbage.patch (Jeroen Demeyer, #15654): When computing a
  determinant(), only collect garbage once per outer loop iteration.
  Better increase PARI stack size instead of collecting garbage too
  often.
* public_memory_functions.patch (Jeroen Demeyer, #16997): Make some of
  PARI's private memory functions public to improve interface with Sage.
