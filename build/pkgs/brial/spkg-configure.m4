SAGE_SPKG_CONFIGURE([brial], [
  dnl Issue #31624: Avoid C++ ABI issues
  SAGE_SPKG_DEPCHECK([gcc boost_cropped m4ri], [
    # If we're using the system m4ri and boost, ensure that we can
    # compile and run an executable linked against both libbrial and
    # libbrial_groebner (both are used by SageMath).
    AC_LANG_PUSH(C++)
    SAVED_LIBS=$LIBS
    LIBS="$LIBS -lbrial -lbrial_groebner"
    AC_MSG_CHECKING([if we can link against brial libraries and run])
    AC_RUN_IFELSE([
      AC_LANG_PROGRAM([
        #include <polybori.h>
        #include <polybori/groebner/groebner_alg.h>
        USING_NAMESPACE_PBORI
        USING_NAMESPACE_PBORIGB

        class MyConstant : public BooleConstant{
          public: void negate() { this->m_value = !this->m_value; }
        };
      ],[
        BoolePolyRing r = BoolePolyRing(2, COrderEnums::dlex);
        ReductionStrategy rs = ReductionStrategy(r);
        rs.llReduceAll(); // uses groebner lib
        if (2 != r.nVariables()) { return 1; }
        if (r.constant(true) == r.constant(false)) { return 2; }
        MyConstant f = MyConstant();
        f.negate(); // ensures v1.1.0+ if m_value isn't const
        if (!f.isOne()) { return 3; }
        return 0;
      ])
    ],
    [
      dnl check we're not on Fedora 30 - more precisely, we reject version 1.2.5
      dnl for which the version is verifiable by the following code.
      AC_MSG_CHECKING([version not equal to 1.2.5])
      AC_RUN_IFELSE([
        AC_LANG_PROGRAM(
          [[#include <polybori.h>
            #include <polybori/config.h>
          ]], [[
             if (VERSION=="1.2.5") return 0;
             else return 1;
        ]])
      ], [
        AC_MSG_RESULT([found a possibly buggy 1.2.5. Rejecting])
        sage_spkg_install_brial=yes
      ], [
        AC_MSG_RESULT([yes])
        sage_spkg_install_brial=no
      ], [
        AC_MSG_RESULT([cross compiling. Assume yes])
        sage_spkg_install_brial=no
      ])
    ],
    [
      AC_MSG_RESULT([no])
      sage_spkg_install_brial=yes
    ],
    [
      AC_MSG_CHECKING([if we can link against brial libraries])
      AC_LINK_IFELSE([
        AC_LANG_PROGRAM([
          #include <polybori.h>
          #include <polybori/groebner/groebner_alg.h>
          USING_NAMESPACE_PBORI
          USING_NAMESPACE_PBORIGB

          class MyConstant : public BooleConstant{
            public: void negate() { this->m_value = !this->m_value; }
          };
        ],[
          BoolePolyRing r = BoolePolyRing(2, COrderEnums::dlex);
          ReductionStrategy rs = ReductionStrategy(r);
          rs.llReduceAll(); // uses groebner lib
          if (2 != r.nVariables()) { return 1; }
          if (r.constant(true) == r.constant(false)) { return 2; }
          MyConstant f = MyConstant();
          f.negate(); // ensures v1.1.0+ if m_value isn't const
          if (!f.isOne()) { return 3; }
          return 0;
        ])
      ],[
        AC_MSG_RESULT([yes])
        sage_spkg_install_brial=yes
      ],[
        AC_MSG_RESULT([no])
        sage_spkg_install_brial=no
      ])
    ])
    LIBS=$SAVED_LIBS
    AC_LANG_POP
  ],
  [ # If we're installing sage's boost or m4ri, then we have to
    # install its BRiAl, too.
    sage_spkg_install_brial=yes
  ])
])
