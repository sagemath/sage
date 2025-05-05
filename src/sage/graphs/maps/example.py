from sage.all import Permutation
from sage.graphs.maps.mutable_labelled_map import MutableLabelledMap
from sage.graphs.maps.map_generator import MapGenerator
from sage.graphs.maps.banner import bannerExampleStart, bannerExampleEnd, mapBanner
from sage.graphs.maps.dynamic_planar_map_show import DynamicPlanarMapShow


class MapExample:
    """
    Class to represents example of use of the library
    """

    def __init__(self) -> None:
        """

        Init the MapExample object

        EXAMPLES::

            sage: MapExample()
            MapExample

        """
        pass

    def __repr__(self) -> str:
        """

        Return the string representation of self

        EXAMPLES::

            sage: MapExample()
            MapExample

        """

        return "MapExample"

    def run(self, false_run=False):
        """
        Run the example.

        INPUT:

            false_run indicate if it is a false run 
            default False
        EXAMPLES::

            sage: MapExample().run(false_run=True)

        """
        if false_run:
            return
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

    def showExample(self, name, myMap, false_run=False):
        """
        show the example.

        INPUT:
            -name, the name of the example
            -myMap, the map used as example
            -false_run, indicate if it is a false run 
            default False

        EXAMPLES::

            sage: m.showExample("A X",m.exampleX(2),false_run=True)

        """

        if false_run:
            return

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
        print("Using dynamicShow...")
        ds = DynamicPlanarMapShow(myMap)
        ds.start(show_halfedges=False)
        print("Done")

        print("Showing the dual of the map")
        print("With the basic show method...")
        myMapDual = myMap.dual()
        myMapDual.show(show_halfedges=False)
        print("Using dynamicShow...")
        ds = DynamicPlanarMapShow(myMapDual)
        ds.start(show_halfedges=False)
        print("Done")

        print("="*100)

    def exampleRandomRootedTree(self, m):
        """

        INPUT:  
        m the size of the tree

        OUTPUT:
        Random rooted tree of size m

        EXAMPLES::

            sage: MapExample().exampleRandomRootedTree(3)
            Rooted map | Sigma : [2, 1, 3, 5, 4, 6] Alpha : [3, 4, 1, 2, 6, 5]

        .. NOTE::
            O(m)
        """
        mapGenerator = MapGenerator()

        return mapGenerator.getRandomTree(m)

    def exampleRandomRootedMap(self, m):
        """

        INPUT:  
        m the size of the planar map

        OUTPUT:
        Random rooted planar map of size m

        EXAMPLES::

            sage: MapExample().exampleRandomRootedMap(3)
            Rooted map | Sigma : [2, 4, 1, 3, 5, 6] Alpha : [3, 5, 1, 6, 2, 4]

        .. NOTE::
            O(m)
        """
        mapGenerator = MapGenerator()

        return mapGenerator.getRandomPlanarMap(m)

    def exampleSimpleGone(self, n):
        """

        INPUT:  
        n the size of the n-gone

        OUTPUT:
        A n gone 

        EXAMPLES::

            sage: MapExample().exampleSimpleGone(3)
            Labelled map | Sigma : [3, 6, 1, 5, 4, 2], Alpha : [2, 1, 4, 3, 6, 5]

        .. NOTE::
            O(nlog(n))
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
        INPUT:
        The number of segment in each  segment of the x

        OUTPUT:
        A x with each segment of the x containing n edges

        EXAMPLES::

            sage: MapExample().exampleX(3)
            Labelled map | Sigma : [3, 2, 9, 5, 4, 7, 6, 8, 15, 11, 10, 13, 12, 14, 1, 17, 16, 19, 18, 20], Alpha : [2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17, 20, 19] 

        .. NOTE::
            O(nlog(n))
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
        INPUT:
            -n, the number of side of the n-gone
            -p, the number of time to repeat the n-gone 

        OUTPUT:

        A repeated  n-gone p times link by edge

        EXAMPLES::

            sage: MapExample().exampleX(3)
            Labelled map | Sigma : [3, 2, 9, 5, 4, 7, 6, 8, 15, 11, 10, 13, 12, 14, 1, 17, 16, 19, 18, 20], Alpha : [2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17, 20, 19] 

        .. NOTE::
            O(pn(log(n)+log(p)))

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
        A.contract()

        return myMap


if __name__ == "__main__":
    mapExample = MapExample()
    mapExample.run()
