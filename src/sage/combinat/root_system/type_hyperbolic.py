"""
Hyperbolic Coxeter matrices for hyperbolic Coxeter types.

These matrices are defined by their position in the Humphreys book "Reflection groups and Coxeter groups". 
The first number in the parenthesis is the page, the second number is the column and the third number is the row.
"""

from sage.combinat.root_system.coxeter_matrix import CoxeterMatrix
from sage.combinat.root_system.coxeter_type import CoxeterType

hyperbolic_coxeter_matrices = {
    (141, 1, 1): CoxeterMatrix([
                [1, 4, 2, 2],
                [4, 1, 3, 2],
                [2, 3, 1, 5],
                [2, 2, 5, 1]
            ]),
    (141, 1, 2): CoxeterMatrix([
                [1, 3, 2, 2],
                [3, 1, 5, 2],
                [2, 5, 1, 3],
                [2, 2, 3, 1]
            ]),
    (141, 1, 3): CoxeterMatrix([
                [1, 5, 2, 2],
                [5, 1, 3, 2],
                [2, 3, 1, 5],
                [2, 2, 5, 1]
            ]),
    (141, 1, 4): CoxeterMatrix([
                [1, 5 ,2 ,2],
                [5, 1, 3, 3],
                [2, 3, 1, 2],
                [2, 3, 2, 1]
            ]),
    (141, 2, 1): CoxeterMatrix([
                [1, 3, 3, 2],
                [3, 1, 2, 4],
                [3, 2, 1, 3],
                [2, 4, 3, 1]
            ]),
    (141, 2, 2): CoxeterMatrix([
                [1, 3, 4, 2],
                [3, 1, 2, 4],
                [4, 2, 1, 3],
                [2, 4, 3, 1]
            ]),
    (141, 2, 3): CoxeterMatrix([
                [1, 3, 5, 2],
                [3, 1, 2, 4],
                [5, 2, 1, 3],
                [2, 4, 3, 1]
            ]),
    (141, 2, 4): CoxeterMatrix([
                [1, 3, 3, 2],
                [3, 1, 2, 5],
                [3, 2, 1, 3],
                [2, 5, 3, 1]
            ]),
    (141, 2, 5): CoxeterMatrix([
                [1, 3, 5, 2],
                [3, 1, 2, 5],
                [5, 2, 1, 3],
                [2, 5, 3, 1]
            ]),
    (141, 3, 1): CoxeterMatrix([
                [1, 4, 2, 2, 2],
                [4, 1, 3, 2, 2],
                [2, 3, 1, 3, 2],
                [2, 2, 3, 1, 5],
                [2, 2, 2, 5, 1]
            ]),
    (141, 3, 2): CoxeterMatrix([
                [1, 3, 2, 2, 2],
                [3, 1, 3, 2, 2],
                [2, 3, 1, 3, 2],
                [2, 2, 3, 1, 5],
                [2, 2, 2, 5, 1]
            ]),
    (141, 3, 3): CoxeterMatrix([
                [1, 5, 2, 2, 2],
                [5, 1, 3, 2, 2],
                [2, 3, 1, 3, 2],
                [2, 2, 3, 1, 5],
                [2, 2, 2, 5, 1]
            ]),
    (141, 3, 4): CoxeterMatrix([
                [1, 5, 2, 2, 2],
                [5, 1, 3, 2, 2],
                [2, 3, 1, 3, 3],
                [2, 2, 3, 1, 2],
                [2, 2, 3, 2, 1]
            ]),
    (141, 3, 5): CoxeterMatrix([
                [1, 3, 3, 2, 2],
                [3, 1, 2, 4, 2],
                [3, 2, 1, 2, 3],
                [2, 4, 2, 1, 3],
                [2, 2, 3, 3, 1]
            ]),
    (142, 1, 1): CoxeterMatrix([
                [1, 4, 3, 2],
                [4, 1, 2, 4],
                [3, 2, 1, 3],
                [2, 4, 3 ,1]
            ]),
    (142, 1, 2): CoxeterMatrix([
                [1, 4, 3, 2],
                [4, 1, 2, 4],
                [3, 2, 1, 4],
                [2, 4, 4, 1]
            ]),
    (142, 1, 3): CoxeterMatrix([
                [1, 4, 4, 2],
                [4, 1, 2, 4],
                [4, 2, 1, 4],
                [2, 4, 4, 1]
            ]),
    (142, 1, 4): CoxeterMatrix([
                [1, 3, 3, 2],
                [3, 1, 2, 6],
                [3, 2, 1, 3],
                [2, 6, 3, 1]
            ]),
    (142, 1, 5): CoxeterMatrix([
                [1, 3, 6, 2],
                [3, 1, 2, 4],
                [6, 2, 1, 3],
                [2, 4, 3, 1]
            ]),
    (142, 1, 6): CoxeterMatrix([
                [1, 3, 6, 2],
                [3, 1, 2, 5],
                [6, 2, 1, 3],
                [2, 5, 3, 1]
            ]),
    (142, 1, 7): CoxeterMatrix([
                [1, 3, 6, 2],
                [3, 1, 2, 6],
                [6, 2, 1, 3],
                [2, 6, 3, 1]
            ]),
    (142, 1, 8): CoxeterMatrix([
                [1, 3, 3, 2],
                [3, 1, 3, 3],
                [3, 3, 1, 3],
                [2, 3, 3, 1]
            ]),
    (142, 1, 9): CoxeterMatrix([
                [1, 3, 2, 2],
                [3, 1, 3, 3],
                [2, 3, 1, 3],
                [2, 3, 3, 1]
            ]),
    (142, 1, 10): CoxeterMatrix([
                [1, 4, 2, 2],
                [4, 1, 3, 3],
                [2, 3, 1, 3],
                [2, 3, 3, 1]
            ]),
    (142, 1, 11): CoxeterMatrix([
                [1, 5, 2, 2],
                [5, 1, 3, 3],
                [2, 3, 1, 3],
                [2, 3, 3, 1]
            ]),
    (142, 1, 12): CoxeterMatrix([
                [1, 6, 2, 2],
                [6, 1, 3, 3],
                [2, 3, 1, 3],
                [2, 3, 3, 1]
            ]),
    (142, 2, 1): CoxeterMatrix([
                [1, 4, 2, 2],
                [4, 1, 4, 2],
                [2, 4, 1, 3],
                [2, 2, 3, 1]
            ]),
    (142, 2, 2): CoxeterMatrix([
                [1, 4, 2, 2],
                [4, 1, 4, 2],
                [2, 4, 1, 4],
                [2, 2, 4, 1]
            ]),
    (142, 2, 3): CoxeterMatrix([
                [1, 4, 2, 2],
                [4, 1, 3, 2],
                [2, 3, 1, 6],
                [2, 2, 6, 1]
            ]),
    (142, 2, 4): CoxeterMatrix([
                [1, 5, 2, 2],
                [5, 1, 3, 2],
                [2, 3, 1, 6],
                [2, 2, 6, 1]
            ]),
    (142, 2, 5): CoxeterMatrix([
                [1, 3, 2, 2],
                [3, 1, 3, 2],
                [2, 3, 1, 6],
                [2, 2, 6, 1]
            ]),
    (142, 2, 6): CoxeterMatrix([
                [1, 3, 2, 2],
                [3, 1, 6, 2],
                [2, 6, 1, 3],
                [2, 2, 3, 1]
            ]),
    (142, 2, 7): CoxeterMatrix([
                [1, 6, 2, 2],
                [6, 1, 3, 2],
                [2, 3, 1, 6],
                [2, 2, 6, 1]
            ]),
    (142, 2, 8): CoxeterMatrix([
                [1, 6, 2, 2],
                [6, 1, 3, 3],
                [2, 3, 1, 2],
                [2, 3, 2, 1]
            ]),
    (142, 2, 9): CoxeterMatrix([
                [1, 3, 2, 2],
                [3, 1, 4, 4],
                [2, 4, 1, 2],
                [2, 4, 2, 1]
            ]),
    (142, 2, 10): CoxeterMatrix([
                [1, 4, 2, 2],
                [4, 1, 4, 4],
                [2, 4, 1, 2],
                [2, 4, 2, 1]
            ]),
    (142, 2, 11): CoxeterMatrix([
                [1, 3, 3, 3],
                [3, 1, 3, 3],
                [3, 3, 1, 3],
                [3, 3, 3, 1]
            ]),
    (142, 3, 1): CoxeterMatrix([
                [1, 3, 2, 2, 2],
                [3, 1, 4, 2, 2],
                [2, 4, 1, 3, 2],
                [2, 2, 3, 1, 4],
                [2, 2, 2, 4, 1]
            ]),
    (142, 3, 2): CoxeterMatrix([
                [1, 3, 2, 2, 2],
                [3, 1, 3, 2, 2],
                [2, 3, 1, 4, 3],
                [2, 2, 4, 1, 2],
                [2, 2, 3, 2, 1]
            ]),
    (142, 3, 3): CoxeterMatrix([
                [1, 3, 2, 2, 2],
                [3, 1, 4, 2, 2],
                [2, 4, 1, 3, 3],
                [2, 2, 3, 1, 2],
                [2, 2, 3, 2, 1]
            ]),
    (142, 3, 4): CoxeterMatrix([
                [1, 4, 2, 2, 2],
                [4, 1, 3, 2, 2],
                [2, 3, 1, 4, 3],
                [2, 2, 4, 1, 2],
                [2, 2, 3, 2, 1]
            ]),
    (142, 3, 5): CoxeterMatrix([
                [1, 3, 2, 2, 2],
                [3, 1, 3, 3, 2],
                [2, 3, 1, 2, 3],
                [2, 3, 2, 1, 3],
                [2, 2, 3, 3, 1]
            ]),
    (142, 3, 6): CoxeterMatrix([
                [1, 4, 2, 2, 2],
                [4, 1, 3, 3, 2],
                [2, 3, 1, 2, 3],
                [2, 3, 2, 1, 3],
                [2, 2, 3, 3, 1]
            ]),
    (142, 3, 7): CoxeterMatrix([
                [1, 2, 4, 2, 2],
                [2, 1, 3, 2, 2],
                [4, 3, 1, 3, 3],
                [2, 2, 3, 1, 2],
                [2, 2, 3, 2, 1]
            ]),
    (142, 3, 8): CoxeterMatrix([
                [1, 3, 2, 3, 2],
                [3, 1, 3, 2, 3],
                [2, 3, 1, 3, 2],
                [3, 2, 3, 1, 3],
                [2, 3, 2, 3, 1]
            ]),
    (142, 3, 9): CoxeterMatrix([
                [1, 4, 3, 2, 2],
                [4, 1, 2, 3, 2],
                [3, 2, 1, 2, 3],
                [2, 3, 2, 1, 4],
                [2, 2, 3, 4, 1]
            ]),
    (143, 1, 1): CoxeterMatrix([
                [1, 3, 2, 2, 2, 2],
                [3, 1, 4, 2, 2, 2],
                [2, 4, 1, 3, 2, 2],
                [2, 2, 3, 1, 3, 2],
                [2, 2, 2, 3, 1, 3],
                [2, 2, 2, 2, 3, 1]
            ]),
    (143, 1, 2): CoxeterMatrix([
                [1, 3, 2, 2, 2, 2],
                [3, 1, 3, 2, 2, 2],
                [2, 3, 1, 4, 2, 2],
                [2, 2, 4, 1, 3, 2],
                [2, 2, 2, 3, 1, 3],
                [2, 2, 2, 2, 3, 1]
            ]),
    (143, 1, 3): CoxeterMatrix([
                [1, 3, 2, 2, 2, 2],
                [3, 1, 4, 2, 2, 2],
                [2, 4, 1, 3, 2, 2],
                [2, 2, 3, 1, 3, 2],
                [2, 2, 2, 3, 1, 4],
                [2, 2, 2, 2, 4, 1]
            ]),
    (143, 1, 4): CoxeterMatrix([
                [1, 3, 2, 2, 2, 2],
                [3, 1, 4, 2, 2, 2],
                [2, 4, 1, 3, 2, 2],
                [2, 2, 3, 1, 3, 3],
                [2, 2, 2, 3, 1, 2],
                [2, 2, 2, 3, 2, 1]
            ]),
    (143, 1, 5): CoxeterMatrix([
                [1, 3, 2, 2, 2, 2],
                [3, 1, 2, 3, 2, 2],
                [2, 2, 1, 3, 2, 2],
                [2, 3, 3, 1, 3, 2],
                [2, 2, 2, 3, 1, 4],
                [2, 2, 2, 2, 4, 1]
            ]),
    (143, 1, 6): CoxeterMatrix([
                [1, 4, 2, 2, 2, 2],
                [4, 1, 2, 3, 2, 2],
                [2, 2, 1, 3, 2, 2],
                [2, 3, 3, 1, 3, 2],
                [2, 2, 2, 3, 1, 4],
                [2, 2, 2, 2, 4, 1]
            ]),
    (143, 1, 7): CoxeterMatrix([
                [1, 3, 2, 2, 2, 2],
                [3, 1, 2, 3, 2, 2],
                [2, 2, 1, 3, 2, 2],
                [2, 3, 3, 1, 3, 3],
                [2, 2, 2, 3, 1, 2],
                [2, 2, 2, 3, 2, 1]
            ]),
    (143, 1, 8): CoxeterMatrix([
                [1, 4, 2, 2, 2, 2],
                [4, 1, 2, 3, 2, 2],
                [2, 2, 1, 3, 2, 2],
                [2, 3 ,3, 1, 3, 3],
                [2, 2, 2, 3, 1, 2],
                [2, 2, 2, 3, 2, 1]
            ]),
    (143, 1, 9): CoxeterMatrix([
                [1, 2, 2, 3, 2, 2],
                [2, 1, 2, 3, 2, 2],
                [2, 2, 1, 3, 2, 2],
                [3, 3, 3, 1, 3, 3],
                [2, 2, 2, 3, 1, 2],
                [2, 2, 2, 3, 2, 1]
            ]),
    (143, 1, 10): CoxeterMatrix([
                [1, 3, 2, 2, 2, 2, 2],
                [3, 1, 2, 3, 2, 2, 2],
                [2, 2, 1, 3, 2, 2, 2],
                [2, 3, 3, 1, 3, 2, 2],
                [2, 2, 2, 3, 1, 3, 2],
                [2, 2, 2, 2, 3, 1, 4],
                [2, 2, 2, 2, 2, 4, 1]
            ]),
    (143, 1, 11): CoxeterMatrix( [
                [1, 2, 3, 2, 2, 2, 2],
                [2, 1, 3, 2, 2, 2, 2],
                [3, 3, 1, 2, 3, 2, 2],
                [2, 2, 2, 1, 3, 2, 2],
                [2, 2, 3, 3, 1, 3, 2],
                [2, 2, 2, 2, 3, 1, 3],
                [2, 2, 2, 2, 2, 3, 1]
            ]),
    (143, 1, 12): CoxeterMatrix([
                [1, 3, 2, 2, 2, 2, 2],
                [3, 1, 3, 3, 2, 2, 2],
                [2, 3, 1, 2, 3, 2, 2],
                [2, 3, 2, 1, 2, 3, 2],
                [2, 2, 3, 2, 1, 2, 3],
                [2, 2, 2, 3, 2, 1, 3],
                [2, 2, 2, 2, 3, 3, 1]
            ]),
    (143, 2, 1): CoxeterMatrix( [
                [1, 3, 3, 2, 2, 2],
                [3, 1, 2, 3, 2, 2],
                [3, 2, 1, 2, 4, 2],
                [2, 3, 2, 1, 2, 3],
                [2, 2, 4, 2, 1, 3],
                [2, 2, 2, 3, 3, 1]
            ]),
    (143, 2, 2): CoxeterMatrix([
                [1, 3, 3, 2, 2, 2],
                [3, 1, 2, 4, 2, 2],
                [3, 2, 1, 2, 4, 2],
                [2, 4, 2, 1, 2, 3],
                [2, 2, 4, 2, 1, 3],
                [2, 2, 2, 3, 3, 1]
            ]),
    (143, 2, 3): CoxeterMatrix([
                [1, 3, 2, 2, 2, 2],
                [3, 1, 3, 3, 2, 2],
                [2, 3, 1, 2, 3, 2],
                [2, 3, 2, 1, 2, 3],
                [2, 2, 3, 2, 1, 3],
                [2, 2, 2, 3, 3, 1]
            ]),
    (144, 1, 1): CoxeterMatrix([
                [1, 3, 2, 2, 2, 2, 2, 2],
                [3, 1, 2, 3, 2, 2, 2, 2],
                [2, 2, 1, 3, 2, 2, 2, 2],
                [2, 3, 3, 1, 3, 2, 2, 2],
                [2, 2, 2, 3, 1, 3, 2, 2],
                [2, 2, 2, 2, 3, 1, 3, 2],
                [2, 2, 2, 2, 2, 3, 1, 4],
                [2, 2, 2, 2, 2, 2, 4, 1]
            ]),
    (144, 1, 2): CoxeterMatrix([
                [1, 2, 3, 2, 2, 2, 2, 2],
                [2, 1, 3, 2, 2, 2, 2, 2],
                [3, 3, 1, 3, 2, 2, 2, 2],
                [2, 2, 3, 1 , 2, 3, 2, 2],
                [2, 2, 2, 2, 1, 3, 2, 2],
                [2, 2, 2, 3, 3, 1, 3, 2],
                [2, 2, 2, 2, 2, 3, 1, 3],
                [2, 2, 2, 2, 2, 2, 3, 1]
            ]),
    (144, 1, 3): CoxeterMatrix([
                [1, 3, 2, 2, 2, 2, 2, 2],
                [3, 1, 3, 2, 2, 2, 2, 2],
                [2, 3, 1, 3, 2, 2, 2, 2],
                [2, 2, 3, 1, 3, 3, 2, 2],
                [2, 2, 2, 3, 1, 2, 3, 2],
                [2, 2, 2, 3, 2, 1, 2, 3],
                [2, 2, 2, 2, 3, 2, 1, 2],
                [2, 2, 2, 2, 2, 3, 2, 1]
            ]),
    (144, 1, 4): CoxeterMatrix([
                [1, 3, 2, 2, 2, 2, 2, 2],
                [3, 1, 3, 3, 2, 2, 2, 2],
                [2, 3, 1, 2, 3, 2, 2, 2],
                [2, 3, 2, 1, 2, 3, 2, 2],
                [2, 2, 3, 2, 1, 2, 3, 2],
                [2, 2, 2, 3, 2, 1, 2, 3],
                [2, 2, 2, 2, 3, 2, 1, 3],
                [2, 2, 2, 2, 2, 3, 3, 1]
            ]),
    (144, 1, 5): CoxeterMatrix([
                [1, 3, 2, 2, 2, 2, 2, 2, 2],
                [3, 1, 3, 2, 2, 2, 2, 2, 2],
                [2, 3, 1, 2, 3, 2, 2, 2, 2],
                [2, 2, 2, 1, 3, 2, 2, 2, 2],
                [2, 2, 3, 3, 1, 3, 2, 2, 2],
                [2, 2, 2, 2, 3, 1, 3, 2, 2],
                [2, 2, 2, 2, 2, 3, 1, 3, 2],
                [2, 2, 2, 2, 2, 2, 3, 1, 3],
                [2, 2, 2, 2, 2, 2, 2, 3, 1]
            ]),
    (144, 1, 6): CoxeterMatrix([
                [1, 3, 2, 2, 2, 2, 2, 2, 2],
                [3, 1, 2, 3, 2, 2, 2, 2, 2],
                [2, 2, 1, 3, 2, 2, 2, 2, 2],
                [2, 3, 3, 1, 3, 2, 2, 2, 2],
                [2, 2, 2, 3, 1, 3, 2, 2, 2],
                [2, 2, 2, 2, 3, 1, 3 ,2 ,2],
                [2, 2, 2, 2, 2, 3, 1, 3, 2],
                [2, 2, 2, 2, 2, 2, 3, 1, 4],
                [2, 2, 2, 2, 2, 2, 2, 4, 1]
            ]),
    (144, 1, 7): CoxeterMatrix([
                [1, 2, 3, 2, 2, 2, 2, 2, 2],
                [2, 1, 3, 2, 2, 2, 2, 2, 2],
                [3, 3, 1, 3, 2, 2, 2, 2, 2],
                [2, 2, 3, 1, 3, 2, 2, 2, 2],
                [2, 2, 2, 3, 1, 2, 3, 2, 2],
                [2, 2, 2, 2, 2, 1, 3, 2, 2],
                [2, 2, 2, 2, 3, 3, 1, 3, 2],
                [2, 2, 2, 2, 2, 2, 3, 1, 3],
                [2, 2, 2, 2, 2, 2, 2, 3, 1]
            ]),
    (144, 1, 8): CoxeterMatrix([
                [1, 3, 2, 2, 2, 2, 2, 2, 2],
                [3, 1, 3, 3, 2, 2, 2, 2, 2],
                [2, 3, 1, 2, 3, 2, 2, 2, 2],
                [2, 3, 2, 1, 2, 3, 2, 2, 2],
                [2, 2, 3, 2, 1, 2, 3, 2, 2],
                [2, 2, 2, 3, 2, 1, 2, 3, 2],
                [2, 2, 2, 2, 3, 2, 1, 2, 3],
                [2, 2, 2, 2, 2, 3, 2, 1, 3],
                [2, 2, 2, 2, 2, 2, 3, 3, 1]
            ]),
    (144, 1, 9): CoxeterMatrix([
                [1, 3, 2, 2, 2, 2, 2, 2, 2, 2],
                [3, 1, 2, 3, 2, 2, 2, 2, 2, 2],
                [2, 2, 1, 3, 2, 2, 2, 2, 2, 2],
                [2, 3, 3, 1, 3, 2, 2, 2, 2, 2],
                [2, 2, 2, 3, 1, 3, 2, 2, 2, 2],
                [2, 2, 2, 2, 3, 1, 3, 2, 2, 2],
                [2, 2, 2, 2, 2, 3, 1, 3, 2, 2],
                [2, 2, 2, 2, 2, 2, 3, 1, 3, 2],
                [2, 2, 2, 2, 2, 2, 2, 3, 1, 3],
                [2, 2, 2, 2, 2, 2, 2, 2, 3, 1]
            ]),
    (144, 1, 10):CoxeterMatrix([
                [1, 3, 2, 2, 2, 2, 2, 2, 2, 2],
                [3, 1, 2, 3, 2, 2, 2, 2, 2, 2],
                [2, 2, 1, 3, 2, 2, 2, 2, 2, 2],
                [2, 3, 3, 1, 3, 2, 2, 2, 2, 2],
                [2, 2, 2, 3, 1, 3, 2, 2, 2, 2],
                [2, 2, 2, 2, 3, 1, 3, 2, 2, 2],
                [2, 2, 2, 2, 2, 3, 1, 3, 2, 2],
                [2, 2, 2, 2, 2, 2, 3, 1, 3, 2],
                [2, 2, 2, 2, 2, 2, 2, 3, 1, 4],
                [2, 2, 2, 2, 2, 2, 2, 2, 4, 1]
            ]),
    (144, 1, 11): CoxeterMatrix([
                [1, 2, 3, 2, 2, 2, 2, 2, 2],
                [2, 1, 3, 2, 2, 2, 2, 2, 2],
                [3, 3, 1, 3, 2, 2, 2, 2, 2],
                [2, 2, 3, 1, 3, 2, 2, 2, 2],
                [2, 2, 2, 3, 1, 2, 3, 2, 2],
                [2, 2, 2, 2, 2, 1, 3, 2, 2],
                [2, 2, 2, 2, 3, 3, 1, 3, 2],
                [2, 2, 2, 2, 2, 2, 3, 1, 3],
                [2, 2, 2, 2, 2, 2, 2, 3, 1]
            ]),
}


class CoxeterType_Hyperbolic(CoxeterType):
    r"""
    Coxeter type hyperbolic.
    """
    def __init__(self, accronym, position):
        """
        EXAMPLES::

            sage: C = CoxeterType(["Hyp", (142, 1, 3)])
            sage: C
            Coxeter type with Humphrey's datum (Page : 142, Column : 1, Row : 3)

        TESTS::

            sage: TestSuite(C).run()
        """
        self._acronym = accronym
        self._position = tuple(position)

        if position not in hyperbolic_coxeter_matrices:
            raise ValueError(f"position {position} is not a valid hyperbolic Coxeter type position")
        super().__init__()

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: CoxeterType(['Hyp', (142, 1, 3)])
            Coxeter type with Humphrey's datum (Page : 142, Column : 1, Row : 3)
        """
        a, b, c = self._position
        return f"Coxeter type with Humphrey's datum (Page : {a}, Column : {b}, Row : {c})"

    def rank(self):
        """
        Return the rank of ``self``.

        This is the number of nodes of the associated Coxeter graph.

        EXAMPLES::

            sage: CoxeterType(['Hyp', (144,1,10)]).rank()
            10
            sage: CoxeterType(['Hyp', (142,2,2)]).rank()
            4
        """
        return hyperbolic_coxeter_matrices[self._position].rank()

    def coxeter_matrix(self):
        """
        Return the Coxeter matrix of ``self``.

        EXAMPLES::

            sage: ct = CartanType(["Hyp", (142, 2, 1)])
            sage: ct.coxeter_matrix()
            [1 3 3 2]
            [3 1 2 4]
            [3 2 1 3]
            [2 4 3 1]
        """
        return hyperbolic_coxeter_matrices[self._position]

    def humphreys_reference(self):
        """
        Return a string with the reference to Humphreys' Reflection groups and Coxeter groups.

        The reference is given by the page, column and row of the table in the book.

        EXAMPLES::

            sage: C = CoxeterType(["Hyp", (142, 1, 3)])
            sage: C.humphreys_reference()
            'Page : 142, Column : 1, Row : 3'
        """
        return self._position

    def coxeter_graph(self):
        """
        Return the Coxeter graph associated to ``self``.

        EXAMPLES::

            sage: C = CoxeterType(["Hyp",(141, 2, 1)])
            sage: C.coxeter_graph()
            Graph on 4 vertices

            sage: C = CoxeterType(["Hyp",(144, 1, 1)])
            sage: C.coxeter_graph()
            Graph on 8 vertices
        """
        return self.coxeter_matrix().coxeter_graph()

    def is_hyperbolic(self):
        """
        Return ``True`` since ``self`` is a hyperbolic Coxeter type.

        EXAMPLES::
            sage: C = CoxeterType(["Hyp", (142, 1, 3)])
            sage: C.is_hyperbolic()
            True

            sage: C = CoxeterType(["A", 2])
            sage: C.is_hyperbolic()
            False
        """
        return True

    def index_set(self):
        """
        Return the index set for ``self``.

        This is the list of the nodes of the associated Coxeter graph.

        EXAMPLES::

            sage: C = CoxeterType(["Hyp",(144, 1, 1)])
            sage: C.index_set()
            (1, 2, 3, 4, 5, 6, 7, 8)

            sage: C = CoxeterType(["Hyp",(143, 1, 6)])
            sage: C.index_set()
            (1, 2, 3, 4, 5, 6)
        """
        return self.coxeter_matrix().index_set()

    def is_affine(self):
        """
        Return ``False`` because hyperbolic Coxeter graphs are never affine.

        EXAMPLES::

            sage: C = CoxeterType(["Hyp",(142, 3, 4)])
            sage: C.is_affine()
            False
        """
        return False

    def is_finite(self):
        """
        Return ``False`` because hyperbolic Coxeter graphs are never finite.

        EXAMPLES::

            sage: C = CoxeterType(["Hyp",(142, 3, 4)])
            sage: C.is_finite()
            False
        """
        return False

    def is_crystallographic(self):
        """
        Return whether ``self`` is crystallographic.

        EXAMPLES::

            sage: C = CoxeterType(["Hyp",(142, 3, 4)])
            sage: C.is_crystallographic()
            True

            sage: C = CoxeterType(["Hyp",(141, 1, 1)])
            sage: C.is_crystallographic()
            False
        """
        return self.coxeter_matrix().is_crystallographic()

    def __eq__(self, other):
        """
        Return whether ``self`` is equal to ``other``.

        EXAMPLES::

            sage: C1 = CoxeterType(["Hyp", (142, 1, 3)])
            sage: C2 = CoxeterType(["Hyp", (142, 1, 3)])
            sage: C1 == C2
            True

            sage: C3 = CoxeterType(["Hyp", (141, 1, 1)])
            sage: C1 == C3
            False
        """
        return isinstance(other, CoxeterType_Hyperbolic) and self._position == other._position
