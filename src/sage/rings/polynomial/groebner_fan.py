r"""
Groebner Fans

Sage provides much of the functionality of ``gfan``, which is a
software package whose main function is to enumerate all reduced
Groebner bases of a polynomial ideal. The reduced Groebner bases
yield the maximal cones in the Groebner fan of the ideal. Several
subcomputations can be issued and additional tools are included.
Among these the highlights are:

-  Commands for computing tropical varieties.

-  Interactive walks in the Groebner fan of an ideal.

-  Commands for graphical renderings of Groebner fans and monomial
   ideals.


AUTHORS:

- Anders Nedergaard Jensen: Wrote the ``gfan`` C++ program, which
  implements algorithms many of which were invented by Jensen, Komei
  Fukuda, and Rekha Thomas. All the underlying hard work of the
  Groebner fans functionality of Sage depends on this C++ program.

- William Stein (2006-04-20): Wrote first version of the Sage code
  for working with Groebner fans.

- Tristram Bogart: the design of the Sage interface
  to ``gfan`` is joint work with Tristram Bogart, who also supplied
  numerous examples.

- Marshall Hampton (2008-03-25): Rewrote various functions to use
  ``gfan-0.3``. This is still a work in progress, comments are
  appreciated on sage-devel@googlegroups.com (or personally at
  hamptonio@gmail.com).

EXAMPLES::

    sage: x,y = QQ['x,y'].gens()
    sage: i = ideal(x^2 - y^2 + 1)
    sage: g = i.groebner_fan()
    sage: g.reduced_groebner_bases()
    [[x^2 - y^2 + 1], [-x^2 + y^2 - 1]]

TESTS::

    sage: x,y = QQ['x,y'].gens()
    sage: i = ideal(x^2 - y^2 + 1)
    sage: g = i.groebner_fan()
    sage: g == loads(dumps(g))
    True

REFERENCES:

- Anders N. Jensen; *Gfan, a software system for Groebner fans*;
  http://home.math.au.dk/jensen/software/gfan/gfan.html
"""

import string
import pexpect
from subprocess import PIPE, Popen

from sage.misc.sage_eval import sage_eval

from sage.structure.sage_object import SageObject
from sage.interfaces.gfan import gfan
from .multi_polynomial_ideal import is_MPolynomialIdeal
from .polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.modules.free_module_element import vector
from sage.plot.all import line, Graphics, polygon
from sage.plot.plot3d.shapes2 import line3d
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.geometry.fan import Fan
from sage.matrix.constructor import matrix
from sage.misc.misc_c import prod
from sage.cpython.string import str_to_bytes, bytes_to_str


def prefix_check(str_list):
    """
    Check if any strings in a list are prefixes of another string in
    the list.

    EXAMPLES::

        sage: from sage.rings.polynomial.groebner_fan import prefix_check
        sage: prefix_check(['z1','z1z1'])
        False
        sage: prefix_check(['z1','zz1'])
        True
    """
    for index1, string1 in enumerate(str_list):
        for index2, string2 in enumerate(str_list):
            if index1 != index2:
                if string1[:len(string2)] == string2:
                    return False
    return True


def max_degree(list_of_polys):
    """
    Compute the maximum degree of a list of polynomials

    EXAMPLES::

        sage: from sage.rings.polynomial.groebner_fan import max_degree
        sage: R.<x,y> = PolynomialRing(QQ,2)
        sage: p_list = [x^2-y,x*y^10-x]
        sage: max_degree(p_list)
        11.0
    """
    return float(max(qf.degree() for qf in list_of_polys))


def _cone_parse(fan_dict_cone):
    """
    Utility function that parses cone information into a dict indexed by
    dimension.

    INPUT:

    -  ``fan_dict_cone`` - the value of a fan_dict with
       key 'CONES'

    EXAMPLES::

        sage: R.<x,y,z,w> = PolynomialRing(QQ,4)
        sage: I = R.ideal([x^2+y^2+z^2-1,x^4-y^2-z-1,x+y+z,w+x+y])
        sage: GF = I.groebner_fan()
        sage: TI = GF.tropical_intersection()
        sage: from sage.rings.polynomial.groebner_fan import _cone_parse
        sage: _cone_parse(TI.fan_dict['CONES'])
        {1: [[0], [1], [2], [3], [4]], 2: [[2, 3]]}
    """
    cone_dict = {}
    cur_dim = 0
    for item in fan_dict_cone:
        temp = item.split('{')[-1]
        if temp.split('}') == '':
            temp = ['']
        else:
            temp = temp.split('}')[0]
            temp = temp.split(' ')
        if item.find('Dimension') != -1:
            cur_dim = Integer(item.split(' ')[-1])
            if cur_dim > 0:
                cone_dict[cur_dim] = [[Integer(q) for q in temp if q != '']]
        else:
            if cur_dim > 0:
                cone_dict[cur_dim] += [[Integer(q) for q in temp if q != '']]
    return cone_dict


class PolyhedralCone(SageObject):

    def __init__(self, gfan_polyhedral_cone, ring=QQ):
        """
        Convert polymake/gfan data on a polyhedral cone into a sage class.

        Currently (18-03-2008) needs a lot of work.

        EXAMPLES::

            sage: R3.<x,y,z> = PolynomialRing(QQ,3)
            sage: gf = R3.ideal([x^8-y^4,y^4-z^2,z^2-2]).groebner_fan()
            sage: a = gf[0].groebner_cone()
            sage: a.facets()
            [[0, 0, 1], [0, 1, 0], [1, 0, 0]]
        """
        cone_keys = ['AMBIENT_DIM', 'DIM', 'IMPLIED_EQUATIONS',
                     'LINEALITY_DIM', 'LINEALITY_SPACE', 'FACETS',
                     'RELATIVE_INTERIOR_POINT']
        poly_lines = gfan_polyhedral_cone.split('\n')
        self.cone_dict = {}
        cur_key = None
        for ting in poly_lines:
            if cone_keys.count(ting):
                cur_key = ting
                self.cone_dict[cur_key] = []
            elif cur_key and ting:
                self.cone_dict[cur_key].append(ting)
        self._facets = []
        for facet in self.cone_dict['FACETS']:
            temp_facet = facet.split('\t')[0]
            temp_facet = temp_facet.split(' ')
            temp_facet = [int(x) for x in temp_facet]
            self._facets.append(temp_facet)
        self._ambient_dim = int(self.cone_dict['AMBIENT_DIM'][0])
        self._dim = int(self.cone_dict['DIM'][0])
        self._lineality_dim = int(self.cone_dict['LINEALITY_DIM'][0])
        rel_int_pt_str = self.cone_dict['RELATIVE_INTERIOR_POINT'][0]
        self._relative_interior_point = [int(q) for q in rel_int_pt_str.split(' ')]

    def _repr_(self):
        """
        Return a basic description of the polyhedral cone.

        EXAMPLES::

            sage: R3.<x,y,z> = PolynomialRing(QQ,3)
            sage: gf = R3.ideal([x^8-y^4,y^4-z^2,z^2-2]).groebner_fan()
            sage: a = gf[0].groebner_cone()
            sage: a # indirect doctests
            Polyhedral cone in 3 dimensions of dimension 3
        """
        return "Polyhedral cone in {} dimensions of dimension {}".format(self.ambient_dim(), self.dim())

    def facets(self):
        """
        Return the inward facet normals of the Groebner cone.

        EXAMPLES::

            sage: R3.<x,y,z> = PolynomialRing(QQ,3)
            sage: gf = R3.ideal([x^8-y^4,y^4-z^2,z^2-2]).groebner_fan()
            sage: a = gf[0].groebner_cone()
            sage: a.facets()
            [[0, 0, 1], [0, 1, 0], [1, 0, 0]]
        """
        return self._facets

    def ambient_dim(self):
        """
        Return the ambient dimension of the Groebner cone.

        EXAMPLES::

            sage: R3.<x,y,z> = PolynomialRing(QQ,3)
            sage: gf = R3.ideal([x^8-y^4,y^4-z^2,z^2-2]).groebner_fan()
            sage: a = gf[0].groebner_cone()
            sage: a.ambient_dim()
            3
        """
        return self._ambient_dim

    def dim(self):
        """
        Return the dimension of the Groebner cone.

        EXAMPLES::

            sage: R3.<x,y,z> = PolynomialRing(QQ,3)
            sage: gf = R3.ideal([x^8-y^4,y^4-z^2,z^2-2]).groebner_fan()
            sage: a = gf[0].groebner_cone()
            sage: a.dim()
            3
        """
        return self._dim

    def lineality_dim(self):
        """
        Return the lineality dimension of the Groebner cone. This is
        just the difference between the ambient dimension and the dimension
        of the cone.

        EXAMPLES::

            sage: R3.<x,y,z> = PolynomialRing(QQ,3)
            sage: gf = R3.ideal([x^8-y^4,y^4-z^2,z^2-2]).groebner_fan()
            sage: a = gf[0].groebner_cone()
            sage: a.lineality_dim()
            0
        """
        return self._lineality_dim

    def relative_interior_point(self):
        """
        Return a point in the relative interior of the Groebner cone.

        EXAMPLES::

            sage: R3.<x,y,z> = PolynomialRing(QQ,3)
            sage: gf = R3.ideal([x^8-y^4,y^4-z^2,z^2-2]).groebner_fan()
            sage: a = gf[0].groebner_cone()
            sage: a.relative_interior_point()
            [1, 1, 1]
        """
        return self._relative_interior_point


class PolyhedralFan(SageObject):
    def __init__(self, gfan_polyhedral_fan, parameter_indices=None):
        """
        Convert polymake/gfan data on a polyhedral fan into a sage class.

        INPUT:

        -  ``gfan_polyhedral_fan`` - output from gfan of a
           polyhedral fan.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: i2 = ideal(x*z + 6*y*z - z^2, x*y + 6*x*z + y*z - z^2, y^2 + x*z + y*z)
            sage: gf2 = i2.groebner_fan(verbose=False)
            sage: pf = gf2.polyhedralfan()
            sage: pf.rays()
            [[-1, 0, 1], [-1, 1, 0], [1, -2, 1], [1, 1, -2], [2, -1, -1]]
        """
        if parameter_indices is None:
            parameter_indices = []
        fan_keys = ['AMBIENT_DIM', 'DIM', 'LINEALITY_DIM', 'RAYS', 'N_RAYS',
                    'LINEALITY_SPACE', 'ORTH_LINEALITY_SPACE', 'F_VECTOR',
                    'CONES', 'MAXIMAL_CONES', 'PURE', 'SIMPLICIAL',
                    'MULTIPLICITIES']
        poly_lines = gfan_polyhedral_fan.split('\n')
        self.fan_dict = {}
        cur_key = None
        for ting in poly_lines:
            if fan_keys.count(ting):
                cur_key = ting
                self.fan_dict[cur_key] = []
            elif cur_key and ting != '':
                self.fan_dict[cur_key].append(ting)
        self._ambient_dim = int(self.fan_dict['AMBIENT_DIM'][0])
        self._dim = int(self.fan_dict['DIM'][0])
        self._lineality_dim = int(self.fan_dict['LINEALITY_DIM'][0])
        self._rays = []
        for ray in self.fan_dict['RAYS']:
            temp_ray = ray.split('\t')[0]
            temp_ray = temp_ray.split(' ')
            temp_ray = [int(x) for x in temp_ray]
            if parameter_indices != []:
                for q in parameter_indices:
                    temp_ray = temp_ray[0:q] + [0] + temp_ray[q:]
            self._rays.append(temp_ray)
        self._cone_dict = _cone_parse(self.fan_dict['CONES'])
        self._maximal_cone_dict = _cone_parse(self.fan_dict['MAXIMAL_CONES'])
        self._str = gfan_polyhedral_fan

    def _repr_(self):
        """
        Return a basic description of the polyhedral fan.

        EXAMPLES::

            sage: R3.<x,y,z> = PolynomialRing(QQ,3)
            sage: gf = R3.ideal([x^8-y^4,y^4-z^2,z^2-2]).groebner_fan()
            sage: pf = gf.polyhedralfan()
            sage: pf # indirect doctest
            Polyhedral fan in 3 dimensions of dimension 3
        """
        return "Polyhedral fan in {} dimensions of dimension {}".format(self.ambient_dim(), self.dim())

    def _str_(self):
        r"""
        Return the raw output of gfan as a string.

        This should only be needed internally as all relevant output
        is converted to sage objects.

        EXAMPLES::

            sage: R3.<x,y,z> = PolynomialRing(QQ,3)
            sage: gf = R3.ideal([x^8-y^4,y^4-z^2,z^2-2]).groebner_fan()
            sage: pf = gf.polyhedralfan()
            sage: pf._str_()
            '_application fan\n_version 2.2\n_type SymmetricFan\n\nAMBIENT_DIM\n3\n\nDIM\n3\n\nLINEALITY_DIM\n0\n\nRAYS\n0 0 1\t# 0\n0 1 0\t# 1\n1 0 0\t# 2\n\nN_RAYS\n3\n\nLINEALITY_SPACE\n\nORTH_LINEALITY_SPACE\n1 0 0\n0 1 0\n0 0 1\n\nF_VECTOR\n1 3 3 1\n\nSIMPLICIAL\n1\n\nPURE\n1\n\nCONES\n{}\t# Dimension 0\n{0}\t# Dimension 1\n{1}\n{2}\n{0 1}\t# Dimension 2\n{0 2}\n{1 2}\n{0 1 2}\t# Dimension 3\n\nMAXIMAL_CONES\n{0 1 2}\t# Dimension 3\n'
        """
        return self._str

    def ambient_dim(self):
        """
        Return the ambient dimension of the Groebner fan.

        EXAMPLES::

            sage: R3.<x,y,z> = PolynomialRing(QQ,3)
            sage: gf = R3.ideal([x^8-y^4,y^4-z^2,z^2-2]).groebner_fan()
            sage: a = gf.polyhedralfan()
            sage: a.ambient_dim()
            3
        """
        return self._ambient_dim

    def dim(self):
        """
        Return the dimension of the Groebner fan.

        EXAMPLES::

            sage: R3.<x,y,z> = PolynomialRing(QQ,3)
            sage: gf = R3.ideal([x^8-y^4,y^4-z^2,z^2-2]).groebner_fan()
            sage: a = gf.polyhedralfan()
            sage: a.dim()
            3
        """
        return self._dim

    def lineality_dim(self):
        """
        Return the lineality dimension of the fan. This is the
        dimension of the largest subspace contained in the fan.

        EXAMPLES::

            sage: R3.<x,y,z> = PolynomialRing(QQ,3)
            sage: gf = R3.ideal([x^8-y^4,y^4-z^2,z^2-2]).groebner_fan()
            sage: a = gf.polyhedralfan()
            sage: a.lineality_dim()
            0
        """
        return self._lineality_dim

    def rays(self):
        """
        A list of rays of the polyhedral fan.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: i2 = ideal(x*z + 6*y*z - z^2, x*y + 6*x*z + y*z - z^2, y^2 + x*z + y*z)
            sage: gf2 = i2.groebner_fan(verbose=False)
            sage: pf = gf2.polyhedralfan()
            sage: pf.rays()
            [[-1, 0, 1], [-1, 1, 0], [1, -2, 1], [1, 1, -2], [2, -1, -1]]
        """
        return sorted(self._rays)

    def cones(self):
        """
        A dictionary of cones in which the keys are the cone dimensions.
        For each dimension, the value is a list of the cones,
        where each element consists of a list of ray indices.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: f = 1+x+y+x*y
            sage: I = R.ideal([f+z*f, 2*f+z*f, 3*f+z^2*f])
            sage: GF = I.groebner_fan()
            sage: PF = GF.tropical_intersection()
            sage: PF.cones()
            {1: [[0], [1], [2], [3], [4], [5]], 2: [[0, 1], [0, 2], [0, 3], [0, 4], [1, 2], [1, 3], [2, 4], [3, 4], [1, 5], [2, 5], [3, 5], [4, 5]]}
        """
        return self._cone_dict

    def maximal_cones(self):
        """
        A dictionary of the maximal cones in which the keys are the
        cone dimensions.  For each dimension, the value is a list of
        the maximal cones, where each element consists of a list of ray indices.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: f = 1+x+y+x*y
            sage: I = R.ideal([f+z*f, 2*f+z*f, 3*f+z^2*f])
            sage: GF = I.groebner_fan()
            sage: PF = GF.tropical_intersection()
            sage: PF.maximal_cones()
            {2: [[0, 1], [0, 2], [0, 3], [0, 4], [1, 2], [1, 3], [2, 4], [3, 4], [1, 5], [2, 5], [3, 5], [4, 5]]}
        """
        return self._maximal_cone_dict

    def f_vector(self):
        """
        The f-vector of the fan.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: f = 1+x+y+x*y
            sage: I = R.ideal([f+z*f, 2*f+z*f, 3*f+z^2*f])
            sage: GF = I.groebner_fan()
            sage: PF = GF.tropical_intersection()
            sage: PF.f_vector()
            [1, 6, 12]
        """
        str_data = self.fan_dict['F_VECTOR'][0]
        fv = [Integer(x) for x in str_data.split(' ')]
        return fv

    def is_simplicial(self):
        """
        Whether the fan is simplicial or not.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: f = 1+x+y+x*y
            sage: I = R.ideal([f+z*f, 2*f+z*f, 3*f+z^2*f])
            sage: GF = I.groebner_fan()
            sage: PF = GF.tropical_intersection()
            sage: PF.is_simplicial()
            True
        """
        return bool(int(self.fan_dict['SIMPLICIAL'][0]))

    def to_RationalPolyhedralFan(self):
        """
        Converts to the RationalPolyhedralFan class, which is more actively
        maintained.  While the information in each class is essentially the
        same, the methods and implementation are different.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: f = 1+x+y+x*y
            sage: I = R.ideal([f+z*f, 2*f+z*f, 3*f+z^2*f])
            sage: GF = I.groebner_fan()
            sage: PF = GF.tropical_intersection()
            sage: fan = PF.to_RationalPolyhedralFan()
            sage: [tuple(q.facet_normals()) for q in fan]
            [(M(0, -1, 0), M(-1, 0, 0)), (M(0, 0, -1), M(-1, 0, 0)), (M(0, 0, 1), M(-1, 0, 0)), (M(0, 1, 0), M(-1, 0, 0)), (M(0, 0, -1), M(0, -1, 0)), (M(0, 0, 1), M(0, -1, 0)), (M(0, 1, 0), M(0, 0, -1)), (M(0, 1, 0), M(0, 0, 1)), (M(1, 0, 0), M(0, -1, 0)), (M(1, 0, 0), M(0, 0, -1)), (M(1, 0, 0), M(0, 0, 1)), (M(1, 0, 0), M(0, 1, 0))]

        Here we use the RationalPolyhedralFan's Gale_transform method on a tropical
        prevariety.

        .. link

        ::

            sage: fan.Gale_transform()
            [ 1  0  0  0  0  1 -2]
            [ 0  1  0  0  1  0 -2]
            [ 0  0  1  1  0  0 -2]

        """
        try:
            return self._fan
        except AttributeError:
            cdnt = []
            cones = self.cones()
            for x in cones:
                if x > 1:
                    cdnt += cones[x]
            fan = Fan(cones=cdnt, rays=self.rays(), discard_faces=True)
            self._fan = fan
            return self._fan


class InitialForm(SageObject):
    def __init__(self, cone, rays, initial_forms):
        """
        A system of initial forms from a polynomial system.

        To each form is associated a cone and a list of
        polynomials (the initial form system itself).

        This class is intended for internal use inside of the
        :class:`TropicalPrevariety` class.

        EXAMPLES::

            sage: from sage.rings.polynomial.groebner_fan import InitialForm
            sage: R.<x,y> = QQ[]
            sage: inform = InitialForm([0], [[-1, 0]], [y^2 - 1, y^2 - 2, y^2 - 3])
            sage: inform._cone
            [0]
        """
        self._cone = cone
        self._rays = rays
        self._initial_forms = initial_forms

    def cone(self):
        """
        The cone associated with the initial form system.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: I = R.ideal([(x+y)^2-1,(x+y)^2-2,(x+y)^2-3])
            sage: GF = I.groebner_fan()
            sage: PF = GF.tropical_intersection()
            sage: pfi0 = PF.initial_form_systems()[0]
            sage: pfi0.cone()
            [0]
        """
        return self._cone

    def rays(self):
        """
        The rays of the cone associated with the initial form system.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: I = R.ideal([(x+y)^2-1,(x+y)^2-2,(x+y)^2-3])
            sage: GF = I.groebner_fan()
            sage: PF = GF.tropical_intersection()
            sage: pfi0 = PF.initial_form_systems()[0]
            sage: pfi0.rays()
            [[-1, 0]]
        """
        return self._rays

    def internal_ray(self):
        """
        A ray internal to the cone associated with the initial form
        system.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: I = R.ideal([(x+y)^2-1,(x+y)^2-2,(x+y)^2-3])
            sage: GF = I.groebner_fan()
            sage: PF = GF.tropical_intersection()
            sage: pfi0 = PF.initial_form_systems()[0]
            sage: pfi0.internal_ray()
            (-1, 0)
        """
        return sum([vector(q) for q in self.rays()])

    def initial_forms(self):
        """
        The initial forms (polynomials).

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: I = R.ideal([(x+y)^2-1,(x+y)^2-2,(x+y)^2-3])
            sage: GF = I.groebner_fan()
            sage: PF = GF.tropical_intersection()
            sage: pfi0 = PF.initial_form_systems()[0]
            sage: pfi0.initial_forms()
            [y^2 - 1, y^2 - 2, y^2 - 3]
        """
        return self._initial_forms


def verts_for_normal(normal, poly):
    """
    Return the exponents of the vertices of a Newton polytope
    that make up the supporting hyperplane for the given outward
    normal.

    EXAMPLES::

        sage: from sage.rings.polynomial.groebner_fan import verts_for_normal
        sage: R.<x,y,z> = PolynomialRing(QQ,3)
        sage: f1 = x*y*z - 1
        sage: f2 = f1*(x^2 + y^2 + 1)
        sage: verts_for_normal([1,1,1],f2)
        [(3, 1, 1), (1, 3, 1)]
    """
    exps = [tuple(x) for x in poly.exponents()]
    expmat = matrix(exps)
    vals = expmat * vector(QQ, normal)
    maxval = max(vals)
    outverts = []
    for i in range(len(exps)):
        if vals[i] == maxval:
            outverts.append(exps[i])
    return outverts


class TropicalPrevariety(PolyhedralFan):

    def __init__(self, gfan_polyhedral_fan, polynomial_system, poly_ring,
                 parameters=None):
        """
        This class is a subclass of the PolyhedralFan class,
        with some additional methods for tropical prevarieties.

        INPUT:

        - ``gfan_polyhedral_fan`` -- output from ``gfan`` of a
          polyhedral fan.
        - ``polynomial_system`` -- a list of polynomials
        - ``poly_ring`` -- the polynomial ring of the list of polynomials
        - ``parameters`` (optional) -- a list of variables to be considered
           as parameters

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: I = R.ideal([(x+y+z)^2-1,(x+y+z)-x,(x+y+z)-3])
            sage: GF = I.groebner_fan()
            sage: TI = GF.tropical_intersection()
            sage: TI._polynomial_system
            [x^2 + 2*x*y + y^2 + 2*x*z + 2*y*z + z^2 - 1, y + z, x + y + z - 3]
        """
        parameter_indices = []
        if parameters is not None:
            allvars = poly_ring.gens()
            parameter_indices = [allvars.index(q) for q in parameters]
        PolyhedralFan.__init__(self, gfan_polyhedral_fan,
                               parameter_indices=parameter_indices)
        self._polynomial_system = polynomial_system
        self._parameters = parameters

    def initial_form_systems(self):
        """
        Return a list of systems of initial forms for each cone
        in the tropical prevariety.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: I = R.ideal([(x+y)^2-1,(x+y)^2-2,(x+y)^2-3])
            sage: GF = I.groebner_fan()
            sage: PF = GF.tropical_intersection()
            sage: pfi = PF.initial_form_systems()
            sage: for q in pfi:
            ....:     print(q.initial_forms())
            [y^2 - 1, y^2 - 2, y^2 - 3]
            [x^2 - 1, x^2 - 2, x^2 - 3]
            [x^2 + 2*x*y + y^2, x^2 + 2*x*y + y^2, x^2 + 2*x*y + y^2]
        """
        try:
            return self._initial_form_systems
        except AttributeError:
            initial_form_systems = []
            pvars = self._polynomial_system[0].parent().gens()
            nvars = len(pvars)
            for dcone in self.cones():
                for acone in self.cones()[dcone]:
                    rays = [self.rays()[i] for i in acone]
                    repray = sum([vector(q) for q in rays])
                    iforms = []
                    for poly in self._polynomial_system:
                        verts = verts_for_normal(repray, poly)
                        nform = 0
                        for x in verts:
                            factorlist = [pvars[i]**x[i] for i in range(nvars)]
                            temp_monomial = prod(factorlist)
                            nform += poly.monomial_coefficient(temp_monomial) * temp_monomial
                        iforms.append(nform)
                    initial_form_systems.append(InitialForm(acone, rays, iforms))
            self._initial_form_systems = initial_form_systems
            return self._initial_form_systems


def ring_to_gfan_format(input_ring):
    """
    Converts a ring to gfan's format.

    EXAMPLES::

        sage: R.<w,x,y,z> = QQ[]
        sage: from sage.rings.polynomial.groebner_fan import ring_to_gfan_format
        sage: ring_to_gfan_format(R)
        'Q[w, x, y, z]'
        sage: R2.<x,y> = GF(2)[]
        sage: ring_to_gfan_format(R2)
        'Z/2Z[x, y]'
    """
    gens = str(input_ring.gens()).replace('(', '[').replace(')', ']')
    if input_ring.base_ring() is QQ:
        return 'Q' + gens
    elif input_ring.base_ring() is ZZ:
        return 'Z' + gens
    else:
        return 'Z/{}Z'.format(input_ring.characteristic()) + gens


def ideal_to_gfan_format(input_ring, polys):
    """
    Return the ideal in gfan's notation.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: polys = [x^2*y - z, y^2*z - x, z^2*x - y]
            sage: from sage.rings.polynomial.groebner_fan import ideal_to_gfan_format
            sage: ideal_to_gfan_format(R, polys)
            'Q[x, y, z]{x^2*y-z,y^2*z-x,x*z^2-y}'

        TESTS:

        Test that :trac:`20146` is fixed::

            sage: P = PolynomialRing(QQ,"x11,x12,x13,x14,x15,x21,x22,x23,x24,x25,x31,x32,x33,x34,x35"); x = P.gens(); M = Matrix(3,x)
            sage: I = P.ideal(M.minors(2))
            sage: ideal_to_gfan_format(P,I.gens())
            'Q[x11, x12, x13, x14, x15, x21, x22, x23, x24, x25, x31, x32, x33, x34, x35]{-x12*x21+x11*x22,-x13*x21+x11*x23,-x14*x21+x11*x24,-x15*x21+x11*x25,-x13*x22+x12*x23,-x14*x22+x12*x24,-x15*x22+x12*x25,-x14*x23+x13*x24,-x15*x23+x13*x25,-x15*x24+x14*x25,-x12*x31+x11*x32,-x13*x31+x11*x33,-x14*x31+x11*x34,-x15*x31+x11*x35,-x13*x32+x12*x33,-x14*x32+x12*x34,-x15*x32+x12*x35,-x14*x33+x13*x34,-x15*x33+x13*x35,-x15*x34+x14*x35,-x22*x31+x21*x32,-x23*x31+x21*x33,-x24*x31+x21*x34,-x25*x31+x21*x35,-x23*x32+x22*x33,-x24*x32+x22*x34,-x25*x32+x22*x35,-x24*x33+x23*x34,-x25*x33+x23*x35,-x25*x34+x24*x35}'
    """
    ideal_gen_str = "{" + ",".join(str(poly).replace(" ", "").replace("'", "")
                                   for poly in polys) + "}"
    ring_str = ring_to_gfan_format(input_ring)
    output = ring_str + ideal_gen_str
    return output


class GroebnerFan(SageObject):

    def __init__(self, I, is_groebner_basis=False, symmetry=None, verbose=False):
        """
        This class is used to access capabilities of the program ``Gfan``.

        In addition to computing Groebner fans, ``Gfan`` can compute
        other things in tropical geometry such as tropical prevarieties.

        INPUT:

        -  ``I`` -- ideal in a multivariate polynomial ring

        -  ``is_groebner_basis`` -- bool (default ``False``). if
           ``True``, then I.gens() must be a Groebner basis with respect to the
           standard degree lexicographic term order.

        -  ``symmetry`` -- default: ``None``; if not ``None``, describes
           symmetries of the ideal

        -  ``verbose`` -- default: ``False``; if ``True``, printout
           useful info during computations

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: I = R.ideal([x^2*y - z, y^2*z - x, z^2*x - y])
            sage: G = I.groebner_fan(); G
            Groebner fan of the ideal:
            Ideal (x^2*y - z, y^2*z - x, x*z^2 - y) of Multivariate Polynomial Ring in x, y, z over Rational Field

        Here is an example of the use of the tropical_intersection command, and then using the RationalPolyhedralFan
        class to compute the Stanley-Reisner ideal of the tropical prevariety::

            sage: R.<x,y,z> = QQ[]
            sage: I = R.ideal([(x+y+z)^3-1,(x+y+z)^3-x,(x+y+z)-3])
            sage: GF = I.groebner_fan()
            sage: PF = GF.tropical_intersection()
            sage: PF.rays()
            [[-1, 0, 0], [0, -1, 0], [0, 0, -1], [1, 1, 1]]
            sage: RPF = PF.to_RationalPolyhedralFan()
            sage: RPF.Stanley_Reisner_ideal(PolynomialRing(QQ,4,'A, B, C, D'))
            Ideal (A*B, A*C, B*C*D) of Multivariate Polynomial Ring in A, B, C, D over Rational Field

        """
        self.__is_groebner_basis = is_groebner_basis
        self.__symmetry = symmetry
        if symmetry:
            print("WARNING! Symmetry option not yet implemented!!")
        self.__verbose = verbose
        if not is_MPolynomialIdeal(I):
            raise TypeError("I must be a multivariate polynomial ideal")
        if not prefix_check([str(R_gen) for R_gen in I.ring().gens()]):
            raise RuntimeError("Ring variables cannot contain each other as prefixes")
        S = I.ring()
        R = S.base_ring()
        # todo: add support for ZZ, which only works for bases computation, not tropical intersections
        if not R.is_field():
            raise NotImplementedError("Groebner fan computation only implemented over fields")
        if not (R is QQ or (R.is_finite() and R.is_prime_field() and R.order() <= 32749)):
            # 32749 is previous_prime(2^15)
            raise NotImplementedError("Groebner fan computation only implemented over Q or GF(p) for p <= 32749.")
        if S.ngens() > 52:
            raise NotImplementedError("Groebner fan computation only implemented for rings in at most 52 variables.")

        self.__ideal = I
        self.__ring = S

    def _repr_(self):
        """
        Describes the Groebner fan and its corresponding ideal.

        EXAMPLES::

            sage: R.<q,u> = PolynomialRing(QQ,2)
            sage: gf = R.ideal([q-u,u^2-1]).groebner_fan()
            sage: gf # indirect doctest
            Groebner fan of the ideal:
            Ideal (q - u, u^2 - 1) of Multivariate Polynomial Ring in q, u over Rational Field
        """
        return "Groebner fan of the ideal:\n{}".format(self.__ideal)

    def __eq__(self, right):
        """
        Tests equality of Groebner fan objects.

        EXAMPLES::

            sage: R.<q,u> = PolynomialRing(QQ,2)
            sage: gf = R.ideal([q^2-u,u^2-q]).groebner_fan()
            sage: gf2 = R.ideal([u^2-q,q^2-u]).groebner_fan()
            sage: gf.__eq__(gf2)
            True
        """
        return type(self) is type(right) and self.ideal() == right.ideal()

    def ideal(self):
        """
        Return the ideal the was used to define this Groebner fan.

        EXAMPLES::

            sage: R.<x1,x2> = PolynomialRing(QQ,2)
            sage: gf = R.ideal([x1^3-x2,x2^3-2*x1-2]).groebner_fan()
            sage: gf.ideal()
            Ideal (x1^3 - x2, x2^3 - 2*x1 - 2) of Multivariate Polynomial Ring in x1, x2 over Rational Field
        """
        return self.__ideal

    def _gfan_maps(self):
        """
        INPUT: none

        OUTPUT:

        - map from Sage ring to ``gfan`` ring

        - map from ``gfan`` ring to Sage ring

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: G = R.ideal([x^2*y - z, y^2*z - x, z^2*x - y]).groebner_fan()
            sage: G._gfan_maps()
            (Ring morphism:
              From: Multivariate Polynomial Ring in x, y, z over Rational Field
              To:   Multivariate Polynomial Ring in a, b, c over Rational Field
              Defn: x |--> a
                    y |--> b
                    z |--> c,
             Ring morphism:
              From: Multivariate Polynomial Ring in a, b, c over Rational Field
              To:   Multivariate Polynomial Ring in x, y, z over Rational Field
              Defn: a |--> x
                    b |--> y
                    c |--> z)
        """
        try:
            return self.__gfan_maps
        except AttributeError:
            S = self.__ring
            n = S.ngens()

            # Define a polynomial ring in n variables
            # that are named a,b,c,d, ..., z, A, B, C, ...
            R = S.base_ring()
            T = PolynomialRing(R, n, string.ascii_letters[:n])

            # Define the homomorphism that sends the
            # generators of S to the generators of T.
            phi = S.hom(T.gens())

            # Define the homomorphism that sends the
            # generators of T to the generators of S.
            psi = T.hom(S.gens())
            self.__gfan_maps = (phi, psi)
            return self.__gfan_maps

    def _gfan_ring(self):
        """
        Return the ring in ``gfan`` notation.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: G = R.ideal([x^2*y - z, y^2*z - x, z^2*x - y]).groebner_fan()
            sage: G._gfan_ring()
            'Q[x, y, z]'
        """
        return ring_to_gfan_format(self.ring())

    def _gfan_ideal(self):
        """
        Return the ideal in ``gfan`` notation.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: G = R.ideal([x^2*y - z, y^2*z - x, z^2*x - y]).groebner_fan()
            sage: G._gfan_ideal()
            'Q[x, y, z]{x^2*y-z,y^2*z-x,x*z^2-y}'
        """
        try:
            return self.__gfan_ideal
        except AttributeError:
            self.__gfan_ideal = ideal_to_gfan_format(self.ring(),
                                                     self.__ideal.gens())
            return self.__gfan_ideal

    def weight_vectors(self):
        """
        Return the weight vectors corresponding to the reduced Groebner
        bases.

        EXAMPLES::

            sage: r3.<x,y,z> = PolynomialRing(QQ,3)
            sage: g = r3.ideal([x^3+y,y^3-z,z^2-x]).groebner_fan()
            sage: g.weight_vectors()
            [(3, 7, 1), (5, 1, 2), (7, 1, 4), (5, 1, 4), (1, 1, 1), (1, 4, 8), (1, 4, 10)]
            sage: r4.<x,y,z,w> = PolynomialRing(QQ,4)
            sage: g4 = r4.ideal([x^3+y,y^3-z,z^2-x,z^3 - w]).groebner_fan()
            sage: len(g4.weight_vectors())
            23
        """
        gfan_processes = Popen(['gfan', '_weightvector', '-m'],
                               stdin=PIPE, stdout=PIPE, stderr=PIPE)
        ans, err = gfan_processes.communicate(input=str_to_bytes(self.gfan()))
        ans = bytes_to_str(ans)
        ans = sage_eval(ans.replace('{', '').replace('}', '').replace('\n', ''))
        return [vector(QQ, x) for x in ans]

    def ring(self):
        """
        Return the multivariate polynomial ring.

        EXAMPLES::

            sage: R.<x1,x2> = PolynomialRing(QQ,2)
            sage: gf = R.ideal([x1^3-x2,x2^3-x1-2]).groebner_fan()
            sage: gf.ring()
            Multivariate Polynomial Ring in x1, x2 over Rational Field
        """
        return self.__ring

    def _gfan_reduced_groebner_bases(self):
        """
        A string of the reduced Groebner bases of the ideal as output by
        ``gfan``.

        EXAMPLES::

            sage: R.<a,b> = PolynomialRing(QQ,2)
            sage: gf = R.ideal([a^3-b^2,b^2-a-1]).groebner_fan()
            sage: gf._gfan_reduced_groebner_bases()
            'Q[a,b]{{b^6-1+2*b^2-3*b^4,a+1-b^2},{a^3-1-a,b^2-1-a}}'
        """
        try:
            return self.__gfan_reduced_groebner_bases
        except AttributeError:
            B = self.gfan(cmd='bases')
            B = B.replace('\n', '')
            self.__gfan_reduced_groebner_bases = B
            return B

    def characteristic(self):
        """
        Return the characteristic of the base ring.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: i1 = ideal(x*z + 6*y*z - z^2, x*y + 6*x*z + y*z - z^2, y^2 + x*z + y*z)
            sage: gf = i1.groebner_fan()
            sage: gf.characteristic()
            0
        """
        return self.__ring.characteristic()

    def reduced_groebner_bases(self):
        """
        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ, 3, order='lex')
            sage: G = R.ideal([x^2*y - z, y^2*z - x, z^2*x - y]).groebner_fan()
            sage: X = G.reduced_groebner_bases()
            sage: len(X)
            33
            sage: X[0]
            [z^15 - z, x - z^9, y - z^11]
            sage: X[0].ideal()
            Ideal (z^15 - z, x - z^9, y - z^11) of Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: X[:5]
            [[z^15 - z, x - z^9, y - z^11],
            [y^2 - z^8, x - z^9, y*z^4 - z, -y + z^11],
            [y^3 - z^5, x - y^2*z, y^2*z^3 - y, y*z^4 - z, -y^2 + z^8],
            [y^4 - z^2, x - y^2*z, y^2*z^3 - y, y*z^4 - z, -y^3 + z^5],
            [y^9 - z, y^6*z - y, x - y^2*z, -y^4 + z^2]]
            sage: R3.<x,y,z> = PolynomialRing(GF(2477),3)
            sage: gf = R3.ideal([300*x^3-y,y^2-z,z^2-12]).groebner_fan()
            sage: gf.reduced_groebner_bases()
            [[z^2 - 12, y^2 - z, x^3 + 933*y],
            [y^4 - 12, x^3 + 933*y, -y^2 + z],
            [x^6 - 1062*z, z^2 - 12, -300*x^3 + y],
            [x^12 + 200, -300*x^3 + y, -828*x^6 + z]]
        """
        try:
            return self.__reduced_groebner_bases
        except AttributeError:
            G = self._gfan_reduced_groebner_bases()
            if G.find(']') != -1:
                G = G.split(']')[1]
            G = G.replace('{{', '').replace('}}', '').split('},{')
            S = self.__ring
            X = [ReducedGroebnerBasis(self, [S(f) for f in G[i].split(',')],
                                      G[i]) for i in range(len(G))]
            self.__reduced_groebner_bases = X
            return X

    def _gfan_mod(self):
        """
        Return the extra options to the ``gfan`` command that are used
        by this object to account for working modulo a prime or in the
        presence of extra symmetries.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: gf = R.ideal([x^3-y,y^3-x-1]).groebner_fan()
            sage: gf._gfan_mod()
            ''
        """
        try:
            return self.__gfan_mod
        except AttributeError:
            mod = ''
            # p = self.characteristic()
            # if p:
            #     mod += ' --mod %s' % p
            # else:
            #     mod += ''

            if self.__is_groebner_basis:
                mod += ' -g'

            if self.__symmetry:
                mod += ' --symmetry'

            self.__gfan_mod = mod
            return self.__gfan_mod

    def gfan(self, cmd='bases', I=None, format=True):
        r"""
        Return the ``gfan`` output as a string given an input ``cmd``.

        The default is to produce the list of reduced Groebner bases
        in ``gfan`` format.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: gf = R.ideal([x^3-y,y^3-x-1]).groebner_fan()
            sage: gf.gfan()
            'Q[x,y]\n{{\ny^9-1-y+3*y^3-3*y^6,\nx+1-y^3}\n,\n{\nx^3-y,\ny^3-1-x}\n,\n{\nx^9-1-x,\ny-x^3}\n}\n'
        """
        if I is None:
            I = self._gfan_ideal()
        # todo -- put something in here (?) when self.__symmetry isn't None...
        cmd += self._gfan_mod()
        s = gfan(I, cmd, verbose=self.__verbose, format=format)
        if s.strip() == '{':
            raise RuntimeError("Error running gfan command %s on %s" % (cmd, self))
        return s

    def __iter__(self):
        """
        Return an iterator for the reduced Groebner bases.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: gf = R.ideal([x^3-y,y^3-x-1]).groebner_fan()
            sage: a = gf.__iter__()
            sage: next(a)
            [y^9 - 3*y^6 + 3*y^3 - y - 1, -y^3 + x + 1]
        """
        for x in self.reduced_groebner_bases():
            yield x

    def __getitem__(self, i):
        """
        Get a reduced groebner basis

        EXAMPLES::

            sage: R4.<w1,w2,w3,w4> = PolynomialRing(QQ,4)
            sage: gf = R4.ideal([w1^2-w2,w2^3-1,2*w3-w4^2,w4^2-w1]).groebner_fan()
            sage: gf[0]
            [w4^12 - 1, -w4^4 + w2, -w4^2 + w1, -1/2*w4^2 + w3]
        """
        return self.reduced_groebner_bases()[i]

    def buchberger(self):
        """
        Return a lexicographic reduced Groebner basis for the ideal.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: G = R.ideal([x - z^3, y^2 - x + x^2 - z^3*x]).groebner_fan()
            sage: G.buchberger()
            [-z^3 + y^2, -z^3 + x]
        """
        try:
            return self.__buchberger
        except AttributeError:
            B = self.gfan(cmd='buchberger')
            if B.find(']') != -1:
                B = B.split(']')[1]
            B = B.replace('}', '').replace('{', '')
            S = self.__ring
            B = [S(f) for f in B.split(',')]
            self.__buchberger = B
            return B

    def polyhedralfan(self):
        """
        Return a polyhedral fan object corresponding to the reduced
        Groebner bases.

        EXAMPLES::

            sage: R3.<x,y,z> = PolynomialRing(QQ,3)
            sage: gf = R3.ideal([x^8-y^4,y^4-z^2,z^2-1]).groebner_fan()
            sage: pf = gf.polyhedralfan()
            sage: pf.rays()
            [[0, 0, 1], [0, 1, 0], [1, 0, 0]]
        """
        try:
            return self.__polyhedralfan
        except AttributeError:
            f = self.gfan(cmd='topolyhedralfan', I=self._gfan_reduced_groebner_bases())
            self.__polyhedralfan = PolyhedralFan(f)
            return self.__polyhedralfan

    def homogeneity_space(self):
        """
        Return the homogeneity space of a the list of polynomials that
        define this Groebner fan.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: G = R.ideal([y^3 - x^2, y^2 - 13*x]).groebner_fan()
            sage: H = G.homogeneity_space()
        """
        try:
            return self.__homogeneity_space
        except AttributeError:
            h = self.gfan(cmd='homogeneityspace', format=False)
            self.__homogeneity_space = h
            return h

    def render(self, file=None, larger=False, shift=0, rgbcolor=(0, 0, 0),
               polyfill=max_degree, scale_colors=True):
        """
        Render a Groebner fan as sage graphics or save as an xfig file.

        More precisely, the output is a drawing of the Groebner fan
        intersected with a triangle. The corners of the triangle are
        (1,0,0) to the right, (0,1,0) to the left and (0,0,1) at the top.
        If there are more than three variables in the ring we extend these
        coordinates with zeros.

        INPUT:

        -  ``file`` - a filename if you prefer the output
           saved to a file. This will be in xfig format.

        -  ``shift`` - shift the positions of the variables in
           the drawing. For example, with shift=1, the corners will be b
           (right), c (left), and d (top). The shifting is done modulo the
           number of variables in the polynomial ring. The default is 0.

        -  ``larger`` - bool (default: ``False``); if ``True``, make
           the triangle larger so that the shape of the Groebner region
           appears. Affects the xfig file but probably not the sage graphics
           (?)

        -  ``rgbcolor`` - This will not affect the saved xfig
           file, only the sage graphics produced.

        -  ``polyfill`` - Whether or not to fill the cones with
           a color determined by the highest degree in each reduced Groebner
           basis for that cone.

        -  ``scale_colors`` - if True, this will normalize
           color values to try to maximize the range


        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: G = R.ideal([y^3 - x^2, y^2 - 13*x,z]).groebner_fan()
            sage: test_render = G.render()

        ::

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: G = R.ideal([x^2*y - z, y^2*z - x, z^2*x - y]).groebner_fan()
            sage: test_render = G.render(larger=True)

        TESTS:

        Testing the case where the number of generators is < 3. Currently,
        this should raise a ``NotImplementedError`` error.

        ::

            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: R.ideal([y^3 - x^2, y^2 - 13*x]).groebner_fan().render()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        S = self.__ring
        if S.ngens() < 3:
            print("For 2-D fan rendering the polynomial ring must have 3 variables (or more, which are ignored).")
            raise NotImplementedError
        cmd = 'render'
        if shift:
            cmd += ' --shiftVariables %s' % shift
        if larger:
            cmd += ' -L'
        s = self.gfan(cmd, I=self._gfan_reduced_groebner_bases().replace(' ', ','), format=False)
        if file is not None:
            with open(file, 'w') as f:
                f.write(s)
        sp = s.split('\n')
        sp2 = []
        for x in sp[9:]:
            xs = x.split(' ')
            y = []
            if x[0:3] != '2 3' and len(xs) > 1:
                for q in xs:
                    if q != '':
                        y.append(q)
                sp2.append(y)
        sp3 = []
        for j in range(len(sp2)):
            temp = []
            for i in range(0, len(sp2[j]) - 1, 2):
                temp.append([float(sp2[j][i]) / 1200.0,
                             float(sp2[j][i + 1]) / 1200.0])
            sp3.append(temp)
        r_lines = Graphics()
        for x in sp3:
            r_lines = r_lines + line(x, rgbcolor=rgbcolor)
        if polyfill:
            vals = [polyfill(q) for q in self.reduced_groebner_bases()]
            if isinstance(vals[0], list) or isinstance(vals[0], tuple):
                if scale_colors:
                    vmins = [min([q[i] for q in vals]) for i in (0, 1, 2)]
                    vmaxs = [max([q[i] for q in vals]) for i in (0, 1, 2)]
                    for i in (0, 1, 2):
                        if vmaxs[i] == vmins[i]:
                            vmaxs[i] = vmins[i] + .01
                    for index in range(len(sp3)):
                        col = [1 - (vals[index][i] - vmins[i]) / (vmaxs[i] - vmins[i]) for i in (0, 1, 2)]
                        r_lines += polygon(sp3[index], rgbcolor=col)
                else:
                    for index in range(len(sp3)):
                        r_lines += polygon(sp3[index], rgbcolor=vals[i])
            else:
                if scale_colors:
                    vmin = min(vals)
                    vmax = max(vals)
                    if vmin == vmax:
                        vmax = vmin + .01
                    for index in range(len(sp3)):
                        r_lines += polygon(sp3[index], hue=.1 + .6 * (vals[index] - vmin) / (vmax - vmin))
                else:
                    for index in range(len(sp3)):
                        r_lines += polygon(sp3[index], hue=vals[i])
        return r_lines

    def _cone_to_ieq(self, facet_list):
        """
        A simple utility function for converting a facet normal to an
        inequality form.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ,2) # dummy stuff to get a gfan object
            sage: gf = R.ideal([x^2+y,y^2]).groebner_fan()
            sage: gf._cone_to_ieq([[1,2,3,4]])
            [[0, 1, 2, 3, 4]]
        """
        ieq_list = []
        for q in facet_list:
            ieq_list.append([0] + q)
        return ieq_list

    def _embed_tetra(self, fpoint):
        """
        Takes a 4-d vector and projects it onto the plane perpendicular to
        (1,1,1,1). Stretches by a factor of 2 as well, since this is only
        for graphical display purposes.

        INPUT:

        -  ``fpoint`` - a list of four numbers

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ,1) # dummy stuff to get a gfan object
            sage: gf = R.ideal([x^2]).groebner_fan()
            sage: gf._embed_tetra([1/2,1/2,1/2,1/2])
            [0, 0, 0]
        """
        v1 = [1, 1, -1, -1]
        v2 = [1, -1, 1, -1]
        v3 = [-1, 1, 1, -1]
        x1 = sum([fpoint[ind] * v1[ind] for ind in range(4)])
        x2 = sum([fpoint[ind] * v2[ind] for ind in range(4)])
        x3 = sum([fpoint[ind] * v3[ind] for ind in range(4)])
        return [x1, x2, x3]

    def _4d_to_3d(self, polyhedral_data):
        """
        A utility function that takes a list of 4d polytopes, projects them
        to 3d, and returns a list of edges.

        INPUT:

        - ``polyhedral_data`` -- an object with 4d vertex and adjacency
          information

        OUTPUT:

        - ``edges`` -- a list of edges in 3d - each list item is a pair of
          points

        EXAMPLES::

            sage: R4.<w,x,y,z> = PolynomialRing(QQ,4)
            sage: gf = R4.ideal([w^2-x,x^2-y,y^2-z,z^2-1]).groebner_fan()
            sage: g_cone = gf[0].groebner_cone()
            sage: g_cone_facets = g_cone.facets()
            sage: g_cone_ieqs = gf._cone_to_ieq(g_cone_facets)
            sage: cone_data = Polyhedron(ieqs = g_cone_ieqs, eqns = [[1,-1,-1,-1,-1]])
            sage: cone_lines = gf._4d_to_3d(cone_data)
            sage: cone_lines
            [[[1, 1, -1], [1, -1/3, 1/3]],
            [[1, 1, -1], [-1/7, 3/7, 5/7]],
            [[1, 1, -1], [-3/5, -1/3, -1/5]],
            [[1, -1/3, 1/3], [-1/7, 3/7, 5/7]],
            [[1, -1/3, 1/3], [-3/5, -1/3, -1/5]],
            [[-1/7, 3/7, 5/7], [-3/5, -1/3, -1/5]]]
        """
        fpoints = polyhedral_data.vertices()
        tpoints = [self._embed_tetra(q) for q in fpoints]
        edges = []
        for vertex in polyhedral_data.vertices():
            i = vertex.index()
            for adjacent_vertex in vertex.adjacent():
                j = adjacent_vertex.index()
                if j > i:
                    try:
                        edges.append([tpoints[i], tpoints[j]])
                    except Exception:
                        print('tpoints: ' + str(tpoints))
                        print('fpoints: ' + str(fpoints))
                        print(adjacent_vertex)
                        print(polyhedral_data.ieqs())
                        raise RuntimeError
        return edges

    def render3d(self, verbose=False):
        """
        For a Groebner fan of an ideal in a ring with four variables, this
        function intersects the fan with the standard simplex perpendicular
        to (1,1,1,1), creating a 3d polytope, which is then projected into
        3 dimensions. The edges of this projected polytope are returned as
        lines.

        EXAMPLES::

            sage: R4.<w,x,y,z> = PolynomialRing(QQ,4)
            sage: gf = R4.ideal([w^2-x,x^2-y,y^2-z,z^2-x]).groebner_fan()
            sage: three_d = gf.render3d()

        TESTS:

        Now test the case where the number of generators is not 4. Currently,
        this should raise a ``NotImplementedError`` error.

        ::

            sage: P.<a,b,c> = PolynomialRing(QQ, 3, order="lex")
            sage: sage.rings.ideal.Katsura(P, 3).groebner_fan().render3d()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        S = self.__ring
        if S.ngens() != 4:
            print("For 3-D fan rendering the polynomial ring must have 4 variables")
            raise NotImplementedError
        g_cones = [q.groebner_cone() for q in self.reduced_groebner_bases()]
        g_cones_facets = [q.facets() for q in g_cones]
        g_cones_ieqs = [self._cone_to_ieq(q) for q in g_cones_facets]
        # Now the cones are intersected with a plane:
        cone_info = [Polyhedron(ieqs=q, eqns=[[1, -1, -1, -1, -1]])
                     for q in g_cones_ieqs]
        # This is really just for debugging
        if verbose:
            for x in cone_info:
                print(x.inequalities() + ([1, 1, 0, 0, 0], [1, 0, 1, 0, 0],
                                          [1, 0, 0, 1, 0], [1, 0, 0, 0, 1]))
                print(x.equations())
                print()
        cone_info = [Polyhedron(ieqs=x.inequalities() +
                                ([1, 1, 0, 0, 0], [1, 0, 1, 0, 0],
                                 [1, 0, 0, 1, 0], [1, 0, 0, 0, 1]),
                                eqns=x.equations()) for x in cone_info]
        all_lines = []
        for cone_data in cone_info:
            try:
                cone_lines = self._4d_to_3d(cone_data)
            except Exception:
                print(cone_data._rays)
                raise RuntimeError
            for a_line in cone_lines:
                all_lines.append(a_line)
        return sum([line3d(a_line) for a_line in all_lines])

    def _gfan_stats(self):
        """
        Return various statistics about this Groebner fan.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: G = R.ideal([y^3 - x^2, y^2 - 13*x]).groebner_fan()
            sage: G._gfan_stats()
            {'Dimension of homogeneity space': 0,
             'Maximal number of polynomials in Groebner basis': 3,
             'Maximal number of terms in Groebner basis': 6,
             'Maximal total degree of a Groebner basis': 4,
             'Minimal total degree of a Groebner basis': 2,
             'Number of reduced Groebner bases': 3,
             'Number of variables': 2}
        """
        try:
            return self.__stats
        except AttributeError:
            s = self.gfan(cmd='stats', I=self._gfan_reduced_groebner_bases().replace(' ', ','), format=False)
            d = {}
            for v in s.split('\n'):
                if v:
                    a, b = v.split(':')
                    d[a] = ZZ(b)
            self.__stats = d
            return d

    def dimension_of_homogeneity_space(self):
        """
        Return the dimension of the homogeneity space.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: G = R.ideal([y^3 - x^2, y^2 - 13*x]).groebner_fan()
            sage: G.dimension_of_homogeneity_space()
            0
        """
        return self._gfan_stats()['Dimension of homogeneity space']

    def maximal_total_degree_of_a_groebner_basis(self):
        """
        Return the maximal total degree of any Groebner basis.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: G = R.ideal([y^3 - x^2, y^2 - 13*x]).groebner_fan()
            sage: G.maximal_total_degree_of_a_groebner_basis()
            4
        """
        return self._gfan_stats()['Maximal total degree of a Groebner basis']

    def minimal_total_degree_of_a_groebner_basis(self):
        """
        Return the minimal total degree of any Groebner basis.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: G = R.ideal([y^3 - x^2, y^2 - 13*x]).groebner_fan()
            sage: G.minimal_total_degree_of_a_groebner_basis()
            2
        """
        return self._gfan_stats()['Minimal total degree of a Groebner basis']

    def number_of_reduced_groebner_bases(self):
        """
        Return the number of reduced Groebner bases.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: G = R.ideal([y^3 - x^2, y^2 - 13*x]).groebner_fan()
            sage: G.number_of_reduced_groebner_bases()
            3
        """
        return self._gfan_stats()['Number of reduced Groebner bases']

    def number_of_variables(self):
        """
        Return the number of variables.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: G = R.ideal([y^3 - x^2, y^2 - 13*x]).groebner_fan()
            sage: G.number_of_variables()
            2

        ::

            sage: R = PolynomialRing(QQ,'x',10)
            sage: R.inject_variables(globals())
            Defining x0, x1, x2, x3, x4, x5, x6, x7, x8, x9
            sage: G = ideal([x0 - x9, sum(R.gens())]).groebner_fan()
            sage: G.number_of_variables()
            10
        """
        return self.__ring.ngens()

    def tropical_basis(self, check=True, verbose=False):
        """
        Return a tropical basis for the tropical curve associated to this
        ideal.

        INPUT:

        -  ``check`` - bool (default: True); if True raises a
           ValueError exception if this ideal does not define a tropical curve
           (i.e., the condition that R/I has dimension equal to 1 + the
           dimension of the homogeneity space is not satisfied).

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ,3, order='lex')
            sage: G = R.ideal([y^3-3*x^2, z^3-x-y-2*y^3+2*x^2]).groebner_fan()
            sage: G
            Groebner fan of the ideal:
            Ideal (-3*x^2 + y^3, 2*x^2 - x - 2*y^3 - y + z^3) of Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: G.tropical_basis()
            [-3*x^2 + y^3, 2*x^2 - x - 2*y^3 - y + z^3, 3/4*x + y^3 + 3/4*y - 3/4*z^3]
        """
        try:
            return self.__tropical_basis
        except AttributeError:
            pass
        cmd = 'tropicalbasis'

        I = self.ideal()
        if not I.is_homogeneous():
            cmd += ' -h'
        if check:
            if I.dimension() != 1 + self.dimension_of_homogeneity_space():
                raise ValueError("The ideal does not define a tropical curve.")

        B = self.gfan(cmd)
        if B.find(']') != -1:
            B = B.split(']')[1]
        S = self.__ring
        B = B.replace('\n', '')
        B = B.replace('{', '').replace('}', '').split(',')
        if verbose:
            print(S, B)
        X = [S(f) for f in B]
        self.__tropical_basis = X
        return X

    def interactive(self, *args, **kwds):
        """
        See the documentation for self[0].interactive(). This does not work
        with the notebook.

        EXAMPLES::

            sage: print("This is not easily doc-testable; please write a good one!")
            This is not easily doc-testable; please write a good one!
        """
        self[0].interactive(*args, **kwds)

    def tropical_intersection(self, parameters=[], symmetry_generators=[], *args, **kwds):
        """
        Return information about the tropical intersection of the
        polynomials defining the ideal.

        This is the common refinement of the outward-pointing normal
        fans of the Newton polytopes of the generators of the
        ideal. Note that some people use the inward-pointing normal
        fans.

        INPUT:

        - ``parameters`` (optional) - a list of variables to be
          considered as parameters
        - ``symmetry_generators`` (optional) - generators of the symmetry group

        OUTPUT: a TropicalPrevariety object

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: I = R.ideal(x*z + 6*y*z - z^2, x*y + 6*x*z + y*z - z^2, y^2 + x*z + y*z)
            sage: gf = I.groebner_fan()
            sage: pf = gf.tropical_intersection()
            sage: pf.rays()
            [[-2, 1, 1]]

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: f1 = x*y*z - 1
            sage: f2 = f1*(x^2 + y^2 + z^2)
            sage: f3 = f2*(x + y + z - 1)
            sage: I = R.ideal([f1,f2,f3])
            sage: gf = I.groebner_fan()
            sage: pf = gf.tropical_intersection(symmetry_generators = '(1,2,0),(1,0,2)')
            sage: pf.rays()
            [[-2, 1, 1], [1, -2, 1], [1, 1, -2]]

            sage: R.<x,y,z> = QQ[]
            sage: I = R.ideal([(x+y+z)^2-1,(x+y+z)-x,(x+y+z)-3])
            sage: GF = I.groebner_fan()
            sage: TI = GF.tropical_intersection()
            sage: TI.rays()
            [[-1, 0, 0], [0, -1, -1], [1, 1, 1]]
            sage: GF = I.groebner_fan()
            sage: TI = GF.tropical_intersection(parameters=[y])
            sage: TI.rays()
            [[-1, 0, 0]]
        """
        try:
            return self.__tropical_intersection
        except AttributeError:
            cmd = 'tropicalintersection'
            id_str = self._gfan_ideal()
            if parameters != []:
                allvars = self.ring().gens()
                truevars = [q for q in allvars if q not in parameters]
                base_ring = self.ring().base_ring()
                new_ring = PolynomialRing(base_ring, len(truevars),
                                          ",".join(str(q) for q in truevars))
                old_polys = self.ideal().gens()
                new_polys = []
                sub = dict([[v, 1] for v in parameters])
                for apoly in old_polys:
                    mons = apoly.monomials()
                    mons = [m.subs(sub) for m in mons]
                    new_polys.append(sum(mons))
                id_str = ideal_to_gfan_format(new_ring, new_polys)
            if symmetry_generators:
                cmd = cmd + ' --symmetryExploit'
                id_str = id_str + '{' + symmetry_generators + '}'
            f = self.gfan(cmd=cmd, I=id_str)
            pf = TropicalPrevariety(f, self.ideal().gens(), self.ring(),
                                    parameters=parameters)
            pf._gfan_output = f
            self.__tropical_intersection = pf
            return pf

    def mixed_volume(self):
        """
        Return the mixed volume of the generators of this ideal.

        This is not really an ideal property, it can depend on the
        generators used.

        The generators must give a square system (as many polynomials
        as variables).

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: example_ideal = R.ideal([x^2-y-1,y^2-z-1,z^2-x-1])
            sage: gf = example_ideal.groebner_fan()
            sage: mv = gf.mixed_volume()
            sage: mv
            8

            sage: R2.<x,y> = QQ[]
            sage: g1 = 1 - x + x^7*y^3 + 2*x^8*y^4
            sage: g2 = 2 + y + 3*x^7*y^3 + x^8*y^4
            sage: example2 = R2.ideal([g1,g2])
            sage: example2.groebner_fan().mixed_volume()
            15
        """
        if len(self.ring().variable_names()) != self.ideal().ngens():
            raise ValueError('not a square system')
        return Integer(self.gfan(cmd='mixedvolume'))


class ReducedGroebnerBasis(SageObject, list):
    def __init__(self, groebner_fan, gens, gfan_gens):
        """
        A class for representing reduced Groebner bases as produced by
        ``gfan``.

        INPUT:

        -  ``groebner_fan`` - a GroebnerFan object from an
           ideal

        -  ``gens`` - the generators of the ideal

        -  ``gfan_gens`` - the generators as a gfan string

        EXAMPLES::

            sage: R.<a,b> = PolynomialRing(QQ,2)
            sage: gf = R.ideal([a^2-b^2,b-a-1]).groebner_fan()
            sage: from sage.rings.polynomial.groebner_fan import ReducedGroebnerBasis
            sage: ReducedGroebnerBasis(gf,gf[0],gf[0]._gfan_gens())
            [b - 1/2, a + 1/2]
        """
        self.__groebner_fan = groebner_fan
        list.__init__(self, gens)
        self.__gfan_gens = '{' + gfan_gens.replace(' ', ',') + '}'
        self.__ring = groebner_fan._gfan_ring()

    def _repr_(self):
        """
        Return the reduced Groebner basis as a string.

        EXAMPLES::

            sage: R.<z1,zz1> = PolynomialRing(QQ,2)
            sage: gf = R.ideal([z1^2*zz1-1,zz1-2]).groebner_fan()
            sage: rgb1 = gf.reduced_groebner_bases()[0]
            sage: rgb1 # indirect doctest
            [zz1 - 2, z1^2 - 1/2]
        """
        return list.__repr__(self)

    def _gfan_gens(self):
        """
        Return the reduced Groebner basis as a string in ``gfan`` format.

        EXAMPLES::

            sage: R.<z1,zz1> = PolynomialRing(QQ,2)
            sage: gf = R.ideal([z1^2*zz1-1,zz1-2]).groebner_fan()
            sage: rgb1 = gf.reduced_groebner_bases()[0]
            sage: rgb1._gfan_gens()
            '{zz1-2,z1^2-1/2}'
        """
        return self.__gfan_gens

    def _gfan(self):
        """
        Return a description of the Groebner fan this basis was derived
        from.

        EXAMPLES::

            sage: R.<z1,zz1> = PolynomialRing(QQ,2)
            sage: gf = R.ideal([z1^2*zz1-1,zz1-2]).groebner_fan()
            sage: rgb1 = gf.reduced_groebner_bases()[0]
            sage: rgb1._gfan()
            Groebner fan of the ideal:
            Ideal (z1^2*zz1 - 1, zz1 - 2) of Multivariate Polynomial Ring in z1, zz1 over Rational Field
        """
        return self.__groebner_fan

    def interactive(self, latex=False, flippable=False, wall=False,
                    inequalities=False, weight=False):
        """
        Do an interactive walk of the Groebner fan starting at this reduced
        Groebner basis.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: G = R.ideal([y^3 - x^2, y^2 - 13*x]).groebner_fan()
            sage: G[0].interactive()      # not tested
            Initializing gfan interactive mode
            *********************************************
            *     Press control-C to return to Sage     *
            *********************************************
            ....
        """
        cmd = 'gfan_interactive'
        if latex:
            cmd += ' -L'
        if flippable:
            cmd += ' -f'
        if wall:
            cmd += ' -w'
        if inequalities:
            cmd += ' -i'
        if weight:
            cmd += ' -W'
        cmd += self.__groebner_fan._gfan_mod()
        E = pexpect.spawn(cmd)
        print("Initializing gfan interactive mode")
        # E.sendline(self._gfan_ideal())
        E.sendline(self.__gfan_gens)
        print("*" * 45)
        print("*     Press control-C to return to Sage     *")
        print("*" * 45)
        try:
            E.interact()
        except OSError:
            print("Returning to Sage.")

    def groebner_cone(self, restrict=False):
        """
        Return defining inequalities for the full-dimensional Groebner cone
        associated to this marked minimal reduced Groebner basis.

        INPUT:

        -  ``restrict`` - bool (default: False); if True, add
           an inequality for each coordinate, so that the cone is restricted
           to the positive orthant.

        OUTPUT: tuple of integer vectors

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: G = R.ideal([y^3 - x^2, y^2 - 13*x]).groebner_fan()
            sage: poly_cone = G[1].groebner_cone()
            sage: poly_cone.facets()
            [[-1, 2], [1, -1]]
            sage: [g.groebner_cone().facets() for g in G]
            [[[0, 1], [1, -2]], [[-1, 2], [1, -1]], [[-1, 1], [1, 0]]]
            sage: G[1].groebner_cone(restrict=True).facets()
            [[-1, 2], [1, -1]]
        """
        try:
            return self.__groebner_cone[restrict]
        except AttributeError:
            self.__groebner_cone = {}
        except KeyError:
            pass
        cmd = 'groebnercone'
        if restrict:
            cmd += ' --restrict'
        gf = self.__groebner_fan
        c = gf.gfan(cmd=cmd, I=self.__ring + self.__gfan_gens)
        return PolyhedralCone(c)

    def ideal(self):
        """
        Return the ideal generated by this basis.

        EXAMPLES::

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: G = R.ideal([x - z^3, y^2 - 13*x]).groebner_fan()
            sage: G[0].ideal()
            Ideal (-13*z^3 + y^2, -z^3 + x) of Multivariate Polynomial Ring in x, y, z over Rational Field
        """
        return self.__groebner_fan.ring().ideal(self)
