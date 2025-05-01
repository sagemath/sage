from re import M
from numpy import kaiser
from sage.all import Permutation
from sage.graphs.planar_maps.MutableLabelledMap import MutableLabelledMap
from sage.graphs.planar_maps.MapGenerator import MapGenerator
from sage.graphs.planar_maps.Banner import bannerExampleStart, bannerExampleEnd, mapBanner
from sage.graphs.planar_maps.DynamicPlanarMapShow import DynamicPlanarMapShow


class MapExample:
    """
    Class to represents example of use of the library
    """

    def __init__(self) -> None:
        pass

    def run(self):

        print(bannerExampleStart)
        print("Starting to show some examples in action")
        self.showExample("Random rooted tree of with 4 edge",
                         self.exampleRandomRootedTree(4))

        self.showExample("Random rooted map with 4 edges",
                         self.exampleRandomRootedMap(4))

        self.showExample("Triangle ", self.exampleSimpleGone(3))

        self.showExample("A X", self.exampleX(2))

        self.showExample("3 triangle linked together",
                         self.exampleRepeatingPolygon(3, 3))
        print(mapBanner)
        print("Image : https://igor-kortchemski.perso.math.cnrs.fr/hdr.pdf")
        print(bannerExampleEnd)

    def showExample(self, name, myMap):
        print("="*100)
        print("Example : ", name)

        print("Number of node:", myMap.n)
        print("Number of edge:", myMap.m)
        print("Number of face", myMap.f)
        print("Genus", myMap.g)

        print("Permutations:")
        myMap.pretty_print()

        print("Showing the map...")
        print("With the basic show method...")
        myMap.show(show_halfedges=False)
        print("Using DynamicPlanarMapShow...")
        ds = DynamicPlanarMapShow(myMap)
        ds.start(show_halfedges=False)
        print("Done")

        print("Showing the dual of the map")
        print("With the basic show method...")
        myMapDual = myMap.dual()
        myMapDual.show(show_halfedges=False)
        print("Using DynamicPlanarMapShow...")
        ds = DynamicPlanarMapShow(myMapDual)
        ds.start(show_halfedges=False)
        print("Done")

        print("="*100)

    def exampleRandomRootedTree(self, m):
        """

        Random rooted tree of size m

        """
        mapGenerator = MapGenerator()

        return mapGenerator.getRandomTree(m)

    def exampleRandomRootedMap(self, m):
        """
        Random rooted map
        """
        mapGenerator = MapGenerator()

        return mapGenerator.getRandomPlanarMap(m)

    def exampleSimpleGone(self, n):
        """
        n-gone
        """
        sigma = Permutation([(1,), (2,)])
        alpha = Permutation([(1, 2)])
        myMap = MutableLabelledMap(sigma=sigma, alpha=alpha)

        A = myMap.X(1)
        B = A.c
        C = A
        for _ in range(n - 1):
            C = C.addEdgeAfter()
        U, _ = C.link(B)
        U.contract()
        return myMap

    def exampleX(self, n):
        """
        A x with each segment of the x  containing n edges
        """
        sigma = Permutation([(1,), (2,)])
        alpha = Permutation([(1, 2)])

        myMap = MutableLabelledMap(sigma=sigma, alpha=alpha)

        A = myMap.X(1)
        for _ in range(3):
            T = A
            for _ in range(n):
                T = T.addEdgeAfter()
            A = A.n

        return myMap

    def exampleRepeatingPolygon(self, n, p):
        """
        Repeat a n-gone p times
        """

        # First draw a n-gone
        polygon = self.exampleSimpleGone(n)

        myMap = polygon.copy()

        A = myMap.X(1)
        A = A.addEdgeBefore()
        B = polygon.X(1)
        for _ in range(p - 1):
            A = A.copyOn(B)[1][0]
            A = A.addEdgeAfter()
        A.delete()

        return myMap


if __name__ == "__main__":
    mapExample = MapExample()
    mapExample.run()
