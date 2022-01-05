SAGE_SPKG_CONFIGURE([libbraiding], [
  # Since libbraiding is a C++ library with no pkg-config file,
  # the best we can do here is compile and run a test program
  # linked against it.
  AC_LANG_PUSH(C++)
  SAVED_LIBS=$LIBS
  LIBS="$LIBS -lbraiding"
  AC_MSG_CHECKING([if we can link against libbraiding])
  AC_RUN_IFELSE([
    AC_LANG_PROGRAM([
      #include <braiding.h>
      #include <list>
      using namespace Braiding;
    ],[
      // Mimic BraidGroup(2)([1,1]).thurston_type() in SageMath.
      // thurstontype == 1 corresponds to "periodic"
      if (thurstontype(2, {1,1}) == 1) { return 0; } else { return 1; }
    ])
  ],
  [
    AC_MSG_RESULT([yes])
    sage_spkg_install_libbraiding=no
  ],
  [
    AC_MSG_RESULT([no])
    sage_spkg_install_libbraiding=yes
  ],[
    AC_LINK_IFELSE([
      AC_LANG_PROGRAM([
        #include <braiding.h>
        #include <list>
        using namespace Braiding;
      ],[
        // Mimic BraidGroup(2)([1,1]).thurston_type() in SageMath.
        // thurstontype == 1 corresponds to "periodic"
        if (thurstontype(2, {1,1}) == 1) { return 0; } else { return 1; }
      ])
    ],
    [
      AC_MSG_RESULT([yes])
      sage_spkg_install_libbraiding=no
    ],
    [
      AC_MSG_RESULT([no])
      sage_spkg_install_libbraiding=yes
    ])
  ])
  LIBS=$SAVED_LIBS
  AC_LANG_POP
])
