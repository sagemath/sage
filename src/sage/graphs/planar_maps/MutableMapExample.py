
from MutableLabelledMap import *


def simpleGone(n):
    """
    n-gone
    """
    sigma = RotatingPermutation([(1,), (2,)])
    alpha = RotatingPermutation([(1, 2)])
    myMap = MutableLabelledMap(sigma=sigma, alpha=alpha)

    A = myMap.X(1)
    B = A.c
    C = A
    for i in range(n - 1):
        C = C.addEdgeAfter()
    U, _ = C.link(B)
    U.contract()
    myMap.pretty_print()
    return myMap


def X(n):
    """
    A x with each segment of the x  containing n edges
    """
    sigma = RotatingPermutation([(1,), (2,)])
    alpha = RotatingPermutation([(1, 2)])

    myMap = MutableLabelledMap(sigma=sigma, alpha=alpha)

    A = myMap.X(1)
    for _ in range(3):
        T = A
        for _ in range(n):
            T = T.addEdgeAfter()
        A = A.n
    myMap.show(use_sage_viewer=False,
               show_halfedges=False, show_vertices=False)


def repeatingPolygon(n, p):
    # First draw a n-gone
    triangle = simpleGone(n)

    myMap = triangle.copy()

    A = myMap.X(1)
    A = A.addEdgeBefore()
    B = triangle.X(1)
    for _ in range(p - 1):
        A = A.copyOn(B)[1][0]
        A = A.addEdgeAfter()
    A.delete()
    myMap.show(use_sage_viewer=False,
               show_halfedges=False, show_vertices=False)


if __name__ == "__main__":
    # Pentagone
    simpleGone(5).show()

    # X
#    X(2)

    # Repeating a triangle 4 times inside each other
#    repeatingPolygon(3, 4)
