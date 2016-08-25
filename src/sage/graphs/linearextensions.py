"""
Linear Extensions of Directed Acyclic Graphs.

A linear extension of a directed acyclic graph is a total (linear) ordering on
the vertices that is compatible with the graph in the following sense:
if there is a path from x to y in the graph, the x appears before y in the
linear extension.

The algorithm implemented in this module is from "Generating Linear Extensions
Fast" by Preusse and Ruskey, which can be found at
http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.52.3057 .  The algorithm
generates the extensions in constant amortized time (CAT) -- a constant amount
of time per extension generated, or linear in the number of extensions
generated.

EXAMPLES:

Here we generate the 5 linear extensions of the following directed
acyclic graph::

    sage: from sage.graphs.linearextensions import LinearExtensions
    sage: D = DiGraph({ 0:[1,2], 1:[3], 2:[3,4] })
    sage: D.is_directed_acyclic()
    True
    sage: LinearExtensions(D).list()
    [[0, 1, 2, 3, 4],
     [0, 1, 2, 4, 3],
     [0, 2, 1, 3, 4],
     [0, 2, 1, 4, 3],
     [0, 2, 4, 1, 3]]

Notice how all of the total orders are compatible with the ordering
induced from the graph.

We can also get at the linear extensions directly from the graph.  From
the graph, the linear extensions are known as topological sorts ::

    sage: D.topological_sort_generator()
    [[0, 1, 2, 3, 4],
     [0, 1, 2, 4, 3],
     [0, 2, 1, 3, 4],
     [0, 2, 1, 4, 3],
     [0, 2, 4, 1, 3]]


"""
#*****************************************************************************
#      Copyright (C) 2008 Mike Hansen <mhansen@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

import sys
from copy import copy
from sage.combinat.combinat import CombinatorialClass


class LinearExtensions(CombinatorialClass):
    def __init__(self, dag):
        r"""
        Creates an object representing the class of all linear extensions
        of the directed acyclic graph \code{dag}.

        Note that upon construction of this object some pre-computation is
        done.  This is the "preprocessing routine" found in Figure 7 of
        "Generating Linear Extensions Fast" by Preusse and Ruskey.

        This is an in-place algorithm and the list self.le keeps track
        of the current linear extensions.  The boolean variable self.is_plus
        keeps track of the "sign".

        EXAMPLES::

            sage: from sage.graphs.linearextensions import LinearExtensions
            sage: D = DiGraph({ 0:[1,2], 1:[3], 2:[3,4] })
            sage: l = LinearExtensions(D)
            sage: l == loads(dumps(l))
            True

        """
        ################
        #Precomputation#
        ################
        dag_copy = copy(dag)
        le = []
        a  = []
        b  = []

        #The preprocessing routine found in Figure 7 of
        #"Generating Linear Extensions Fast" by
        #Pruesse and Ruskey
        while dag_copy.num_verts() != 0:
            #Find all the minimal elements of dag_copy
            minimial_elements = []
            for node in dag_copy.vertices():
                if len(dag_copy.incoming_edges(node)) == 0:
                    minimial_elements.append(node)
            if len(minimial_elements) == 1:
                le.append(minimial_elements[0])
                dag_copy.delete_vertex(minimial_elements[0])
            else:
                ap = minimial_elements[0]
                bp = minimial_elements[1]
                a.append(ap)
                b.append(bp)
                le.append(ap)
                le.append(bp)
                dag_copy.delete_vertex(ap)
                dag_copy.delete_vertex(bp)
        self.max_pair = len(a) - 1

        self.le = le
        self.a  = a
        self.b  = b
        self.dag = dag
        self.mrb = 0
        self.mra = 0
        self.is_plus = True
        self.linear_extensions = None
        self._name = "Linear extensions of %s"%dag


    def switch(self, i):
        """
        This implements the Switch procedure described on page 7
        of "Generating Linear Extensions Fast" by Pruesse and Ruskey.

        If i == -1, then the sign is changed.  If i > 0, then self.a[i]
        and self.b[i] are transposed.

        Note that this meant to be called by the generate_linear_extensions
        method and is not meant to be used directly.

        EXAMPLES::

            sage: from sage.graphs.linearextensions import LinearExtensions
            sage: D = DiGraph({ 0:[1,2], 1:[3], 2:[3,4] })
            sage: l = LinearExtensions(D)
            sage: _ = l.list()
            sage: l.le = [0, 1, 2, 3, 4]
            sage: l.is_plus
            True
            sage: l.switch(-1)
            sage: l.is_plus
            False
            sage: l.a
            [1, 4]
            sage: l.b
            [2, 3]
            sage: l.switch(0)
            sage: l.le
            [0, 2, 1, 3, 4]
            sage: l.a
            [2, 4]
            sage: l.b
            [1, 3]


        """
        if i == -1:
            self.is_plus = not self.is_plus
        if i >= 0:
            a_index = self.le.index(self.a[i])
            b_index = self.le.index(self.b[i])
            self.le[a_index] = self.b[i]
            self.le[b_index] = self.a[i]

            self.b[i], self.a[i] = self.a[i], self.b[i]

        if self.is_plus:
            self.linear_extensions.append(self.le[:])


    def move(self, element, direction):
        """
        This implements the Move procedure described on page 7
        of "Generating Linear Extensions Fast" by Pruesse and Ruskey.

        If direction is "left", then this transposes element with the
        element on its left.  If the direction is "right", then this
        transposes element with the element on its right.

        Note that this is meant to be called by the generate_linear_extensions
        method and is not meant to be used directly.

        EXAMPLES::

            sage: from sage.graphs.linearextensions import LinearExtensions
            sage: D = DiGraph({ 0:[1,2], 1:[3], 2:[3,4] })
            sage: l = LinearExtensions(D)
            sage: _ = l.list()
            sage: l.le = [0, 1, 2, 3, 4]
            sage: l.move(1, "left")
            sage: l.le
            [1, 0, 2, 3, 4]
            sage: l.move(1, "right")
            sage: l.le
            [0, 1, 2, 3, 4]

        """
        index = self.le.index(element)
        if direction == "right":
            self.le[index] = self.le[index+1]
            self.le[index+1] = element
        elif direction == "left":
            self.le[index] = self.le[index-1]
            self.le[index-1] = element
        else:
            print("Bad direction!")
            sys.exit()
        if self.is_plus:
            self.linear_extensions.append(self.le[:])


    def right(self, i, letter):
        """
        If letter =="b", then this returns True if and only if
        self.b[i] is incomparable with the elements to its right
        in self.le.  If letter == "a", then it returns True if
        and only if self.a[i] is incomparable with the element to its
        right in self.le and the element to the right is not
        self.b[i]

        This is the Right function described on page 8 of
        "Generating Linear Extensions Fast" by Pruesse and Ruskey.

        Note that this is meant to be called by the generate_linear_extensions
        method and is not meant to be used directly.

        EXAMPLES::

            sage: from sage.graphs.linearextensions import LinearExtensions
            sage: D = DiGraph({ 0:[1,2], 1:[3], 2:[3,4] })
            sage: l = LinearExtensions(D)
            sage: _ = l.list()
            sage: l.le
            [0, 1, 2, 4, 3]
            sage: l.a
            [1, 4]
            sage: l.b
            [2, 3]
            sage: l.right(0, "a")
            False
            sage: l.right(1, "a")
            False
            sage: l.right(0, "b")
            False
            sage: l.right(1, "b")
            False

        """
        if letter == "a":
            x = self.a[i]
            yindex = self.le.index(x) + 1
            if yindex >= len(self.le):
                return False
            y = self.le[ yindex ]
            if self.incomparable(x,y) and y != self.b[i]:
                return True
            return False
        elif letter == "b":
            x = self.b[i]
            yindex = self.le.index(x) + 1
            if yindex >= len(self.le):
                return False
            y = self.le[ yindex ]
            if self.incomparable(x,y):
                return True
            return False
        else:
            raise ValueError("Bad letter!")

    def generate_linear_extensions(self, i):
        """
        This a Python version of the GenLE routine found in Figure 8
        of "Generating Linear Extensions Fast" by Pruesse and Ruskey.

        Note that this is meant to be called by the list
        method and is not meant to be used directly.

        EXAMPLES::

            sage: from sage.graphs.linearextensions import LinearExtensions
            sage: D = DiGraph({ 0:[1,2], 1:[3], 2:[3,4] })
            sage: l = LinearExtensions(D)
            sage: l.linear_extensions = []
            sage: l.linear_extensions.append(l.le[:])
            sage: l.generate_linear_extensions(l.max_pair)
            sage: l.linear_extensions
            [[0, 1, 2, 3, 4], [0, 2, 1, 3, 4]]

        """
        if i >= 0:
            self.generate_linear_extensions(i-1)
            mrb = 0
            typical = False
            while self.right(i, "b"):
                mrb += 1
                self.move(self.b[i], "right")
                self.generate_linear_extensions(i-1)
                mra = 0
                if self.right(i, "a"):
                    typical = True
                    cont = True
                    while cont:
                        mra += 1
                        self.move(self.a[i], "right")
                        self.generate_linear_extensions(i-1)
                        cont = self.right(i, "a")
                if typical:
                    self.switch(i-1)
                    self.generate_linear_extensions(i-1)
                    if mrb % 2 == 1:
                        mla = mra -1
                    else:
                        mla = mra + 1
                    for x in range(mla):
                        self.move(self.a[i], "left")
                        self.generate_linear_extensions(i-1)

            if typical and (mrb % 2 == 1):
                self.move(self.a[i], "left")
            else:
                self.switch(i-1)
            self.generate_linear_extensions(i-1)
            for x in range(mrb):
                self.move(self.b[i], "left")
                self.generate_linear_extensions(i-1)

    def list(self):
        """
        Returns a list of the linear extensions of the directed acyclic graph.

        Note that once they are computed, the linear extensions are
        cached in this object.

        EXAMPLES::

            sage: from sage.graphs.linearextensions import LinearExtensions
            sage: D = DiGraph({ 0:[1,2], 1:[3], 2:[3,4] })
            sage: LinearExtensions(D).list()
            [[0, 1, 2, 3, 4],
             [0, 1, 2, 4, 3],
             [0, 2, 1, 3, 4],
             [0, 2, 1, 4, 3],
             [0, 2, 4, 1, 3]]
        """
        if self.linear_extensions is not None:
            return self.linear_extensions[:]

        self.linear_extensions = []
        self.linear_extensions.append(self.le[:])
        self.generate_linear_extensions(self.max_pair)
        self.switch(self.max_pair)
        self.generate_linear_extensions(self.max_pair)
        self.linear_extensions.sort()
        return self.linear_extensions[:]


    def incomparable(self, x, y):
        """
        Returns True if vertices x and y are incomparable in the directed
        acyclic graph when thought of as a poset.

        EXAMPLES::

            sage: from sage.graphs.linearextensions import LinearExtensions
            sage: D = DiGraph({ 0:[1,2], 1:[3], 2:[3,4] })
            sage: l = LinearExtensions(D)
            sage: l.incomparable(0,1)
            False
            sage: l.incomparable(1,2)
            True
        """
        if (not self.dag.shortest_path(x, y)) and (not self.dag.shortest_path(y, x)):
            return True
        return False
