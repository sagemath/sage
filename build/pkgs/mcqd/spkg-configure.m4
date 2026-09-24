SAGE_SPKG_CONFIGURE([mcqd], [
 SAGE_SPKG_DEPCHECK([cmake], [
   AS_TMPDIR([foobar42mcqd])
   cmake -Sbuild/pkgs/mcqd -B$tmp
   if test $? -ne 0; then
       sage_spkg_install_mcqd=yes
   fi
 ])
])
