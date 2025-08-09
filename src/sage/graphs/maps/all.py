from sage.misc.lazy_import import lazy_import

lazy_import("sage.graphs.maps.labelled_map", "LabelledMap")
lazy_import("sage.graphs.maps.mutable_labelled_map", "MutableLabelledMap")
lazy_import("sage.graphs.maps.primitive_mutable_labelled_map", "PrimitiveMutableLabelledMap")
lazy_import("sage.graphs.maps.rooted_map", "RootedMap")
lazy_import("sage.graphs.maps.dynamic_planar_map_show", "DynamicPlanarMapShow")
lazy_import("sage.graphs.maps.example", "MapExample")
lazy_import("sage.graphs.maps.map_generator", "MapGenerator")

del lazy_import
