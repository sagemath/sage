project(
    'mpfi',
    'c',
    version: '1.5.4'
)

c = meson.get_compiler('c')

# Dependencies
gmp = dependency('gmp', version: '>= 4.1.0')
mpfr = dependency('mpfr', version: '>= 4.0.1')

# Configuration data
conf = configuration_data()
conf.set('PACKAGE_VERSION', '"' + meson.project_version() + '"')
# Check for functions
conf.set('HAVE_DUP2', c.has_function('dup2') ? 1 : 0)
conf.set('HAVE_GETTIMEOFDAY', c.has_function('gettimeofday') ? 1 : 0)
conf.set('HAVE_MPFR_Q_SUB', c.has_function('mpfr_q_sub') ? 1 : 0)
conf.set('HAVE_MPFR_Z_SUB', c.has_function('mpfr_z_sub') ? 1 : 0)
conf.set('HAVE_MPFR_Z_DIV', c.has_function('mpfr_z_div') ? 1 : 0)

configure_file(
  output: 'mpfi_config.h',
  configuration: conf
)

# Source files
src = files(
    'src/abs.c', 'src/acos.c', 'src/acosh.c', 'src/add.c',
    'src/add_d.c', 'src/add_fr.c', 'src/add_q.c', 'src/add_si.c',
    'src/add_ui.c', 'src/add_z.c', 'src/alea.c', 'src/asin.c',
    'src/asinh.c', 'src/atan2.c', 'src/atan.c', 'src/atanh.c',
    'src/bisect.c', 'src/blow.c', 'src/cbrt.c', 'src/clear.c',
    'src/clears.c', 'src/cmp.c', 'src/cmp_sym_pi.c', 'src/constants.c',
    'src/cos.c', 'src/cosh.c', 'src/csc.c', 'src/csch.c',
    'src/cot.c', 'src/coth.c', 'src/d_div.c', 'src/d_sub.c',
    'src/diam.c', 'src/div.c', 'src/div_2exp.c', 'src/div_2si.c',
    'src/div_2ui.c', 'src/div_d.c', 'src/div_ext.c', 'src/div_fr.c',
    'src/div_q.c', 'src/div_si.c', 'src/div_ui.c', 'src/div_z.c',
    'src/erandom.c', 'src/error.c', 'src/exp.c', 'src/exp2.c',
    'src/exp10.c', 'src/expm1.c', 'src/exp2m1.c', 'src/exp10m1.c',
    'src/fr_div.c', 'src/fr_sub.c', 'src/get_endpoints.c', 'src/get_fr.c',
    'src/get_d.c', 'src/get_prec.c', 'src/get_version.c', 'src/has_zero.c',
    'src/hypot.c', 'src/increase.c', 'src/init.c', 'src/init2.c',
    'src/inits.c', 'src/inits2.c', 'src/intersect.c', 'src/interv_d.c',
    'src/interv_fr.c', 'src/interv_q.c', 'src/interv_si.c', 'src/interv_ui.c',
    'src/interv_z.c', 'src/inp_str.c', 'src/inv.c', 'src/is_empty.c',
    'src/is_inside.c', 'src/log.c', 'src/log10.c', 'src/log1p.c',
    'src/log2.c', 'src/log2p1.c', 'src/log10p1.c', 'src/mag.c',
    'src/mid.c', 'src/mig.c', 'src/mul.c', 'src/mul_2exp.c',
    'src/mul_2si.c', 'src/mul_2ui.c', 'src/mul_d.c', 'src/mul_fr.c',
    'src/mul_q.c', 'src/mul_si.c', 'src/mul_ui.c', 'src/mul_z.c',
    'src/neg.c', 'src/nrandom.c', 'src/out_str.c', 'src/predicates.c',
    'src/print_binary.c', 'src/put.c', 'src/put_d.c', 'src/put_fr.c',
    'src/put_q.c', 'src/put_si.c', 'src/put_ui.c', 'src/put_z.c',
    'src/q_div.c', 'src/q_sub.c', 'src/quadrant.c', 'src/rec_sqrt.c',
    'src/revert_if_needed.c', 'src/round_prec.c', 'src/sec.c', 'src/sech.c',
    'src/set.c', 'src/set_d.c', 'src/set_flt.c', 'src/set_fr.c',
    'src/set_ld.c', 'src/set_prec.c', 'src/set_q.c', 'src/set_si.c',
    'src/set_str.c', 'src/set_ui.c', 'src/set_z.c', 'src/si_div.c',
    'src/si_sub.c', 'src/sign.c', 'src/sin.c', 'src/sinh.c',
    'src/sqr.c', 'src/sqrt.c', 'src/sub.c', 'src/sub_d.c',
    'src/sub_fr.c', 'src/sub_q.c', 'src/sub_si.c', 'src/sub_ui.c',
    'src/sub_z.c', 'src/swap.c', 'src/tan.c', 'src/tanh.c',
    'src/ui_div.c', 'src/ui_sub.c', 'src/union.c', 'src/urandom.c',
    'src/z_div.c', 'src/z_sub.c'
)

# Library
libmpfi = static_library(
    'mpfi',
    src,
    include_directories: include_directories('.', 'src'),
    dependencies: [gmp, mpfr],
    install: true,
)

# Install headers
install_headers('src/mpfi.h', 'src/mpfi_io.h', subdir: 'mpfi')

# Pkg-config file
# pkgconfig = import('pkgconfig')
# pkgconfig.generate(
#     name: 'mpfi',
#     description: 'MPFI library',
#     version: meson.project_version(),
#     libraries: libmpfi,
#     subdirs: 'mpfi',
#     install_dir: join_paths(get_option('libdir'), 'pkgconfig'),
# )

mpfi_dep = declare_dependency(
    include_directories: include_directories('src'),
    link_with: libmpfi,
    dependencies: [gmp, mpfr],
)

# Extra distribution files
# dist_files = [
#     'mpfi_config.h.in',
#     'AUTHORS',
#     'COPYING',
#     'COPYING.LESSER',
#     'NEWS',
#     'TODO',
# ]
# foreach file : dist_files
#     meson.add_dist_script('cp', file, meson.current_build_dir())
# endforeach

# Subdirectories
#subdir('doc')
#subdir('src')
#subdir('tests')
