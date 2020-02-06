SAGE_SPKG_CONFIGURE([boost], [
 SAGE_SPKG_DEPCHECK([boost_cropped], [
  AC_RUN_IFELSE([dnl an extra sanity check
  AC_LANG_PROGRAM(
    [[#include <boost/program_options/errors.hpp>
    ]], [[
      boost::program_options::error err("Error message");
      return 0;
    ]])], [], [sage_spkg_install_boost=yes])
 ])
])
