from sage.misc.lazy_import import lazy_import


lazy_import("sage.graphs.planar_maps.Banner", ["mapBanner", "bannerExampleStart", "bannerExampleEnd"])
lazy_import("sage.graphs.planar_maps.LabelledMap", "LabelledMap")
lazy_import("sage.graphs.planar_maps.MutableLabelledMap", "MutableLabelledMap")
lazy_import("sage.graphs.planar_maps.PrimitiveMutableLabelledMap", "PrimitiveMutableLabelledMap")
lazy_import("sage.graphs.planar_maps.RootedMap", "RootedMap")
lazy_import("sage.graphs.planar_maps.DynamicPlanarMapShow", "DynamicPlanarMapShow")
lazy_import("sage.graphs.planar_maps.example", "MapExample")
lazy_import("sage.graphs.planar_maps.MapGenerator", "MapGenerator")

del lazy_import
