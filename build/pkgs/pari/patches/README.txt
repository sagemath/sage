======================================================================
Current patches to PARI in Sage:
======================================================================

Patches to configuration files:
* get_ld.patch (Carl Witty): Disable -rpath.
* get_tests.patch (John Cremona): Disable testing of ellglobalred in
  "make test-all" in spkg-check, since it requires the elldata database
  which we do not include since it is about 18MB.
* get_config_options.patch (Leif Leonhardy): Catch invalid arguments to
  "--graphic" (and treat such as an error) since otherwise strange
  compilation errors might occur (cf. #9722, too).
* get_fltk.patch (Leif Leonhardy): Explicitly link against libstdc++
  when using FLTK (for plotting) to support Fedora 13, and do an extra
  check for the FLTK include dir (cf. #9722).
* get_X11.patch (Leif Leonhardy): Also search */lib64/* for X11
  libraries (since on some systems this is really a separate directory,
  i.e. neither a symbolic link to */lib/* nor the target of a symbolic
  link */lib/*; cf. #9722, too).
* get_dlcflags.patch (Jeroen Demeyer): Add -fno-common to DLCFLAGS on
  Darwin. Submitted upstream, but upstream only applied it for PowerPC.
  Since this doesn't break anything and only improves performance, add
  the flag unconditionally.
* install_doc_no_make.patch (Jeroen Demeyer): Do not *build* the
  documentation when doing install-doc or install-docpdf.  We must not
  build the documentation because that requires tex.  On the other hand,
  to have ? and ?? work within gp, we must install the .tex files (but
  not .dvi files).  So simply not doing install-doc doesn't work.
* perl_path.patch (Jeroen Demeyer): change first line of all perl
  scripts to "#!/usr/bin/env perl" (#10559).
* cygwin_dll_a.patch (Jean-Pierre Flori): copy libpari.dll.a on Cygwin
  (see #13333).
* KERNELCFLAGS.patch (Jeroen Demeyer): when SAGE_DEBUG=yes, compile
  kernel files with -O1 instead of -funroll-loops; -O0 gives a
  segmentation fault on some OS X systems when doing
  factor(10356613*10694706299664611221)
  See #13921, also reported upstream:
  - http://pari.math.u-bordeaux.fr/archives/pari-dev-1301/msg00000.html

C files:
* src/kernel/gmp/mp.c (Leif Leonhardy):
  Do not override GMP's memory functions.
  In addition, let PARI use "GMP internals" (access members of GMP
  structures directly) *conditionally*. (We *don't* disable that by
  default, since at least currently this is compatible with both GMP
  *and* MPIR. To disable the use, add "-DPARI_DONT_USE_GMP_INTERNALS"
  to CFLAGS. This is a trivial patch to mp.c only; see also the comment
  there.):
* GCC_PR49330.patch (Jeroen Demeyer): in pari_init_functions(), reorder
  code to work around a GCC bug concerning pointer arithmetic:
  http://gcc.gnu.org/bugzilla/show_bug.cgi?id=49330
  This bug manifests itself as a Bus Error on OS X 10.4 PPC with
  gcc-4.6.3. Discussed with upstream:
  - http://pari.math.u-bordeaux.fr/archives/pari-dev-1203/msg00001.html
  - http://pari.math.u-bordeaux.fr/archives/pari-dev-1202/msg00006.html
* trac_13902_determinant.patch (Jeroen Demeyer): patch backported from
  upstream git repository based on commits
  - 28ea998bc661f5bbde18b6d6b0f50111a10ae16c
  - 249432f7088bfa114ed5cd3a5d64ef51ee968e35
* polred.patch (Jeroen Demeyer): Fix polred(), add polredbest()
  backported from upstream commits
  - 2d00a24fbb1ffe8eba35b9a04763e36eef8a5f7b
  - a3596c56f9439144a0dbed4c47bd6ff9485e3fc8
  - 1a00ca416de4daebccaab2be1a4b8a061a9f2fde
  - ad550d9bbfee8113087407c3262bffc27a020c98
* q_eint1.patch (Jeroen Demeyer): Fix precision of eint1(), backported
  from upstream commits
  - b894e0171112493bb4945a6a0e4af039b54a1d1f
  - 54eb37055167dbd8d4856b1bd698cf7018db561b
  and completely backport new algorithm for mpveceint1(). The prefix
  q_ in the patch filename is to ensure it comes after polred.patch.
