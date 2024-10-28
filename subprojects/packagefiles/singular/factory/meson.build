# Dependencies
gmp = dependency('gmp', required: true)
flint = dependency('flint', required: true)
ntl = dependency('ntl', required: false)

# Configuration data
fs = import('fs')
gen_headers = []
conf = configuration_data()
# Use no streamio (C++)
conf.set('NOSTREAMIO', not cpp.has_header('iostream.h'))#, 'strstream.h', 'fstream.h', 'iostream', 'string', 'fstream', 'ctype.h'))
conf.set('FACTORYVERSION', meson.project_version())
factory_configuration = '@0@ in @1@'.format(meson.project_version(), meson.source_root())
conf.set_quoted('FACTORYCONFIGURATION', factory_configuration)

gen_headers += fs.copyfile('factory.template', 'factory.h')#, install: true, install_dir: 'factory')
gen_headers += fs.copyfile('factoryconf.template', 'factoryconf.h')#, install: true, install_dir: 'factory')

# https://mesonbuild.com/Wrap-best-practices-and-tips.html#do-not-put-configh-in-external-search-path
subdir('internal')
internal_inc = include_directories('internal')

# Include directories
inc_dirs = [include_directories('include', 'include/factory', '.'), src, internal_inc]

# Compiler flags
cpp_args = [
  '-DFACTORY_BUILDING_DLL',
]

# Source files
sources = files(
  'canonicalform.cc',
  'cf_algorithm.cc',
  'cf_char.cc',
  'cfCharSets.cc',
  'cfCharSetsUtil.cc',
  'cf_chinese.cc',
  'cf_cyclo.cc',
  'cf_eval.cc',
  'cfEzgcd.cc',
  'cf_factor.cc',
  'cf_factory.cc',
  'cf_gcd.cc',
  'cfGcdAlgExt.cc',
  'cfGcdUtil.cc',
  'cf_generator.cc',
  'cf_globals.cc',
  'cf_hnf.cc',
  'cf_inline.cc',
  'cf_irred.cc',
  'cf_iter.cc',
  'cf_iter_inline.cc',
  'cf_linsys.cc',
  'cf_map.cc',
  'cf_map_ext.cc',
  'cfModGcd.cc',
  'cfNewtonPolygon.cc',
  'cfNTLzzpEXGCD.cc',
  'cfModResultant.cc',
  'cf_ops.cc',
  'cf_primes.cc',
  'cf_random.cc',
  'cf_resultant.cc',
  'cf_reval.cc',
  'cfSubResGcd.cc',
  'cf_switches.cc',
  'cf_util.cc',
  'cf_binom.cc',
  'cfUnivarGcd.cc',
  'cf_roots.cc',
  'debug.cc',
  'DegreePattern.cc',
  'ExtensionInfo.cc',
  'facAbsBiFact.cc',
  'facAbsFact.cc',
  'facAlgExt.cc',
  'facAlgFunc.cc',
  'facAlgFuncUtil.cc',
  'facBivar.cc',
  'facFactorize.cc',
  'fac_sqrfree.cc',
  'fac_util.cc',
  'fac_berlekamp.cc',
  'fac_cantzass.cc',
  'fac_univar.cc',
  'fac_multivar.cc',
  'fac_multihensel.cc',
  'fac_iterfor.cc',
  'fac_distrib.cc',
  'facFqBivar.cc',
  'facFqBivarUtil.cc',
  'facFqFactorize.cc',
  'facFqFactorizeUtil.cc',
  'facFqSquarefree.cc',
  'facHensel.cc',
  'facIrredTest.cc',
  'facMul.cc',
  'facSparseHensel.cc',
  'ffops.cc',
  'FLINTconvert.cc',
  'gf_tabutil.cc',
  'gfops.cc',
  'imm.cc',
  'int_cf.cc',
  'int_int.cc',
  'int_intdiv.cc',
  'int_poly.cc',
  'int_rat.cc',
  'int_pp.cc',
  'variable.cc',
  'NTLconvert.cc',
  'singext.cc',
  'parseutil.cc',
  'ftmpl_inst.cc',
  'threadsupport.cc',
)

# Library
libfactory = library(
  'factory',
  sources, gen_headers,
  include_directories: inc_dirs,
  dependencies: [gmp, flint, ntl],
  cpp_args: cpp_args,
  install: true,
)

# Executable for generating GF(q)-tables
# gengftables = executable(
#   'gengftables',
#   'gengftables-conway.cc',
#   include_directories: inc_dirs,
#   dependencies: [gmp, flint, ntl],
#   install: true,
# )

# # Test program
# test_prog = executable(
#   'test_prog',
#   'test.cc',
#   include_directories: inc_dirs,
#   dependencies: [gmp, flint, ntl],
#   install: false,
# )

# test('factory test', test_prog)

# Install headers
install_headers(
  [
    'cf_assert.h',
    'canonicalform.h',
    'cf_algorithm.h',
    'cfCharSets.h',
    'cfCharSetsUtil.h',
    'cf_cyclo.h',
    'cf_defs.h',
    'cf_eval.h',
    'cfEzgcd.h',
    'cf_factory.h',
    'cf_generator.h',
    'cf_globals.h',
    'cfGcdAlgExt.h',
    'cfGcdUtil.h',
    'cf_hnf.h',
    'cf_irred.h',
    'cf_iter.h',
    'cf_map.h',
    'cf_map_ext.h',
    'cfModGcd.h',
    'cfNewtonPolygon.h',
    'cfNTLzzpEXGCD.h',
    'cfModResultant.h',
    'cf_primes.h',
    'cf_primetab.h',
    'cf_random.h',
    'cf_reval.h',
    'cfSubResGcd.h',
    'cf_switches.h',
    'cf_util.h',
    'cf_binom.h',
    'cfUnivarGcd.h',
    'cf_roots.h',
    'debug.h',
    'DegreePattern.h',
    'ExtensionInfo.h',
    'facAbsBiFact.h',
    'facAbsFact.h',
    'facAlgExt.h',
    'facAlgFunc.h',
    'facAlgFuncUtil.h',
    'facBivar.h',
    'facFactorize.h',
    'fac_sqrfree.h',
    'fac_util.h',
    'fac_berlekamp.h',
    'fac_cantzass.h',
    'fac_univar.h',
    'fac_multivar.h',
    'fac_multihensel.h',
    'fac_iterfor.h',
    'fac_distrib.h',
    'facFqBivar.h',
    'facFqBivarUtil.h',
    'facFqFactorize.h',
    'facFqFactorizeUtil.h',
    'facFqSquarefree.h',
    'facHensel.h',
    'facIrredTest.h',
    'facMul.h',
    'facSparseHensel.h',
    'ffops.h',
    'FLINTconvert.h',
    'gf_tabutil.h',
    'gfops.h',
    'gmpext.h',
    'imm.h',
    'int_cf.h',
    'int_int.h',
    'int_poly.h',
    'int_rat.h',
    'int_pp.h',
    'timing.h',
    'variable.h',
    'NTLconvert.h',
    'singext.h',
    'parseutil.h',
  ],
  subdir: 'factory',
)

# Install GF(q)-tables
gftables = [
  'gftables/10201', 'gftables/1024', 'gftables/10609', 'gftables/11449',
  'gftables/11881', 'gftables/121', 'gftables/12167', 'gftables/125',
  'gftables/12769', 'gftables/128', 'gftables/1331', 'gftables/1369',
  'gftables/14641', 'gftables/15625', 'gftables/16', 'gftables/16129',
  'gftables/16384', 'gftables/16807', 'gftables/1681', 'gftables/169',
  'gftables/17161', 'gftables/1849', 'gftables/18769', 'gftables/19321',
  'gftables/19683', 'gftables/2048', 'gftables/2187', 'gftables/2197',
  'gftables/2209', 'gftables/22201', 'gftables/22801', 'gftables/2401',
  'gftables/243', 'gftables/24389', 'gftables/24649', 'gftables/25',
  'gftables/256', 'gftables/26569', 'gftables/27', 'gftables/27889',
  'gftables/2809', 'gftables/28561', 'gftables/289', 'gftables/29791',
  'gftables/29929', 'gftables/3125', 'gftables/32', 'gftables/32041',
  'gftables/32761', 'gftables/32768', 'gftables/343', 'gftables/3481',
  'gftables/361', 'gftables/36481', 'gftables/3721', 'gftables/37249',
  'gftables/38809', 'gftables/39601', 'gftables/4', 'gftables/4096',
  'gftables/44521', 'gftables/4489', 'gftables/49', 'gftables/4913',
  'gftables/49729', 'gftables/5041', 'gftables/50653', 'gftables/512',
  'gftables/51529', 'gftables/52441', 'gftables/529', 'gftables/5329',
  'gftables/54289', 'gftables/57121', 'gftables/58081', 'gftables/59049',
  'gftables/6241', 'gftables/625', 'gftables/63001', 'gftables/64',
  'gftables/6561', 'gftables/6859', 'gftables/6889', 'gftables/729',
  'gftables/7921', 'gftables/8', 'gftables/81', 'gftables/8192',
  'gftables/841', 'gftables/9', 'gftables/9409', 'gftables/961',
]

install_data(gftables, install_dir: join_paths(get_option('datadir'), 'factory/gftables'))

# Pkg-config file
# pkgconfig = import('pkgconfig')
# pkgconfig.generate(
#   name: 'factory',
#   description: 'Factory library',
#   version: '0.1.0',
#   libraries: libfactory,
#   subdirs: 'factory',
#   install_dir: join_paths(get_option('libdir'), 'pkgconfig'),
# )

factory_dep = declare_dependency(
  include_directories: [include_directories('.', 'include'), src],
  link_with: libfactory,
  dependencies: [gmp, flint, ntl],
  sources: gen_headers,
)
