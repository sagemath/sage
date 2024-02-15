headers = files(
    'L.h',
    'Lcomplex.h',
    'Ldokchitser.h',
    'Lexplicit_formula.h',
    'Lgamma.h',
    'Lgram.h',
    'Lnumbertheory.h',
    'Lprint.h',
    'Lriemannsiegel_blfi.h',
    'Lcommon.h',
    'Ldirichlet_series.h',
    'Lelliptic.h',
    'Lfind_zeros.h',
    'Lglobals.h',
    'Lmisc.h',
    'Lnumberzeros.h',
    'Lriemannsiegel.h',
    'Lvalue.h'
) + config_h
install_headers(headers, subdir : 'lcalc')
subdir('lcalc')

# TODO: if PRECISION_MULTIPLE
# pkginclude_HEADERS += mpreal.h Lmpreal.h
# endif

sources = [
    'Ldokchitser.cc',
    'Lelliptic.cc',
    'Lgamma.cc',
    'Lglobals.cc',
    'Lmisc.cc',
    'Lnumbertheory.cc',
    'Lriemannsiegel.cc'
]

flags = [
    '-version-info 1:1:0',
    # The -no-undefined flag is required to create a Windows DLL, because
    # that platform does not support shared libraries with undefined
    # symbols in them. libLfunction, however, does not rely on undefined
    # symbols in the first place. So, we append the flag unconditionally.
    '-no-undefined'
]

lfunction = shared_library('Lfunction', sources,
    #link_args : flags, # TODO: Implement? Remove?
    include_directories : include_directories('.'),
    dependencies : [config_h_dep],
    install : true
)
lcalc_dep = declare_dependency(
    link_with : lfunction,
    include_directories: include_directories('.'),
)
