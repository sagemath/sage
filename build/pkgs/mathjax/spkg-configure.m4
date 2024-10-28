SAGE_SPKG_CONFIGURE([mathjax], [
  # Arch:   /usr/share/mathjax
  # Gentoo: /usr/share/mathjax
  # Void:   /usr/share/mathjax
  AC_MSG_CHECKING([for MathJax-3.x])
  m4_foreach([mathjax_dir], [/usr/share/mathjax], [
    # tex-chtml.hs is used in src/sage_docbuild/conf.py
    # and was not present in MathJax-2.x
    AS_IF([test -f "mathjax_dir/tex-chtml.js"], [
      SAGE_MATHJAX_DIR="mathjax_dir"
      AC_MSG_RESULT([mathjax_dir])
    ])
  ])
  AS_IF([test -z "${SAGE_MATHJAX_DIR}"], [
    AC_MSG_RESULT([no])
    sage_spkg_install_mathjax=yes
  ])
],[],[],[
  # post-check
  AS_IF([test x$sage_spkg_install_mathjax = xyes], [
    # Our spkg-src script adds an extra "mathjax"
    SAGE_MATHJAX_DIR='${prefix}'/share/mathjax/mathjax
  ])

  AC_SUBST(SAGE_MATHJAX_DIR, "${SAGE_MATHJAX_DIR}")
])
