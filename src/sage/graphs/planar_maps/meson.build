py.install_sources(
  '__init__.py',
  'all.py',
  'Banner.py',
  'CustomSwap.py',
  'CycleUtilsProvider.py',
  'CyclicChainedList.py',
  'DynamicPlanarMapShow.py',
  'example.py',
  'LabelledMap.py',
  'MapDecorator.py',
  'MapError.py',
  'MapGenerator.py',
  'MapPermutation.py',
  'MutableLabelledMap.py',
  'MutableMapExample.py',
  'MutableTopologicalDemiEdge.py',
  'PermutationUtilsAbstractor.py',
  'PrimitiveMutableTopologicalDemiEdge.py',
  'PrimitiveRotatingPermutation.py',
  'PrimitiveRotatingPermutationUtilsAbstractor.py',
  'RootedMap.py',
  'RotatingPermutation.py',
  'RotatingPermutationUtilsAbstractor.py',
  'SplayTree.py',
  'TopologicalDemiEdge.py',
  subdir: 'sage/graphs/planar_maps',
)

extension_data = {}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/graphs/planar_maps',
    install: true,
    include_directories: [],
    dependencies: [py_dep, cysignals, gmp],
  )
endforeach

