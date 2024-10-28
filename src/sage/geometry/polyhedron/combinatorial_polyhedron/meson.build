py.install_sources(
  'all.py',
  'base.pxd',
  'combinatorial_face.pxd',
  'conversions.pxd',
  'face_data_structure.pxd',
  'face_iterator.pxd',
  'face_list_data_structure.pxd',
  'list_of_faces.pxd',
  'polyhedron_face_lattice.pxd',
  subdir: 'sage/geometry/polyhedron/combinatorial_polyhedron',
)

extension_data = {
  'base' : files('base.pyx'),
  'combinatorial_face' : files('combinatorial_face.pyx'),
  'conversions' : files('conversions.pyx'),
  'face_iterator' : files('face_iterator.pyx'),
  'face_list_data_structure' : files('face_list_data_structure.pyx'),
  'list_of_faces' : files('list_of_faces.pyx'),
  'polyhedron_face_lattice' : files('polyhedron_face_lattice.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/geometry/polyhedron/combinatorial_polyhedron',
    install: true,
    include_directories: [inc_cpython, inc_data_structures, inc_rings],
    dependencies: [py_dep, cysignals, gmp],
  )
endforeach

