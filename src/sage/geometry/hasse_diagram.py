"""
Construction of finite atomic and coatomic lattices from incidences

This module provides the function :func:`lattice_from_incidences` for
computing finite atomic and coatomic lattices in the sense of
partially ordered sets where any two elements have meet and joint. For
example, the face lattice of a polyhedron.
"""

#*****************************************************************************
#       Copyright (C) 2010 Andrey Novoseltsev <novoselt@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.graphs.digraph import DiGraph
from sage.combinat.posets.lattices import FiniteLatticePoset


def lattice_from_incidences(atom_to_coatoms, coatom_to_atoms,
                            face_constructor=None,
                            required_atoms=None,
                            key=None,
                            **kwds):
    r"""
    Compute an atomic and coatomic lattice from the incidence between
    atoms and coatoms.

    INPUT:

    - ``atom_to_coatoms`` -- list, ``atom_to_coatom[i]`` should list all
      coatoms over the ``i``-th atom;

    - ``coatom_to_atoms`` -- list, ``coatom_to_atom[i]`` should list all
      atoms under the ``i``-th coatom;

    - ``face_constructor`` -- function or class taking as the first two
      arguments sorted :class:`tuple` of integers and any keyword arguments.
      It will be called to construct a face over atoms passed as the first
      argument and under coatoms passed as the second argument. Default
      implementation will just return these two tuples as a tuple;

    - ``required_atoms`` -- list of atoms (default:None). Each
      non-empty "face" requires at least one of the specified atoms
      present. Used to ensure that each face has a vertex.

    - ``key`` -- any hashable value (default: None). It is passed down
      to :class:`~sage.combinat.posets.posets.FinitePoset`.

    - all other keyword arguments will be passed to ``face_constructor`` on
      each call.

    OUTPUT:

    - :class:`finite poset <sage.combinat.posets.posets.FinitePoset>` with
      elements constructed by ``face_constructor``.

    .. NOTE::

        In addition to the specified partial order, finite posets in Sage have
        internal total linear order of elements which extends the partial one.
        This function will try to make this internal order to start with the
        bottom and atoms in the order corresponding to ``atom_to_coatoms`` and
        to finish with coatoms in the order corresponding to
        ``coatom_to_atoms`` and the top. This may not be possible if atoms and
        coatoms are the same, in which case the preference is given to the
        first list.

    ALGORITHM:

    The detailed description of the used algorithm is given in [KP2002]_.

    The code of this function follows the pseudo-code description in the
    section 2.5 of the paper, although it is mostly based on frozen sets
    instead of sorted lists - this makes the implementation easier and should
    not cost a big performance penalty. (If one wants to make this function
    faster, it should be probably written in Cython.)

    While the title of the paper mentions only polytopes, the algorithm (and
    the implementation provided here) is applicable to any atomic and coatomic
    lattice if both incidences are given, see Section 3.4.

    In particular, this function can be used for strictly convex cones and
    complete fans.

    REFERENCES: [KP2002]_

    AUTHORS:

    - Andrey Novoseltsev (2010-05-13) with thanks to Marshall Hampton for the
      reference.

    EXAMPLES:

    Let us construct the lattice of subsets of {0, 1, 2}.
    Our atoms are {0}, {1}, and {2}, while our coatoms are {0,1}, {0,2}, and
    {1,2}. Then incidences are ::

        sage: atom_to_coatoms = [(0,1), (0,2), (1,2)]
        sage: coatom_to_atoms = [(0,1), (0,2), (1,2)]

    and we can compute the lattice as ::

        sage: from sage.geometry.cone import lattice_from_incidences
        sage: L = lattice_from_incidences(
        ....:                     atom_to_coatoms, coatom_to_atoms)
        sage: L
        Finite lattice containing 8 elements with distinguished linear extension
        sage: for level in L.level_sets(): print(level)
        [((), (0, 1, 2))]
        [((0,), (0, 1)), ((1,), (0, 2)), ((2,), (1, 2))]
        [((0, 1), (0,)), ((0, 2), (1,)), ((1, 2), (2,))]
        [((0, 1, 2), ())]

    For more involved examples see the *source code* of
    :meth:`sage.geometry.cone.ConvexRationalPolyhedralCone.face_lattice` and
    :meth:`sage.geometry.fan.RationalPolyhedralFan._compute_cone_lattice`.
    """

    def default_face_constructor(atoms, coatoms, **kwds):
        return (atoms, coatoms)
    if face_constructor is None:
        face_constructor = default_face_constructor
    atom_to_coatoms = [frozenset(atc) for atc in atom_to_coatoms]
    A = frozenset(range(len(atom_to_coatoms)))  # All atoms
    coatom_to_atoms = [frozenset(cta) for cta in coatom_to_atoms]
    C = frozenset(range(len(coatom_to_atoms)))  # All coatoms
    # Comments with numbers correspond to steps in Section 2.5 of the article
    L = DiGraph(1)       # 3: initialize L
    faces = {}
    atoms = frozenset()
    coatoms = C
    faces[atoms, coatoms] = 0
    next_index = 1
    Q = [(atoms, coatoms)]              # 4: initialize Q with the empty face
    while Q:                            # 5
        q_atoms, q_coatoms = Q.pop()    # 6: remove some q from Q
        q = faces[q_atoms, q_coatoms]
        # 7: compute H = {closure(q+atom) : atom not in atoms of q}
        H = {}
        candidates = set(A.difference(q_atoms))
        for atom in candidates:
            coatoms = q_coatoms.intersection(atom_to_coatoms[atom])
            atoms = A
            for coatom in coatoms:
                atoms = atoms.intersection(coatom_to_atoms[coatom])
            H[atom] = (atoms, coatoms)
        # 8: compute the set G of minimal sets in H
        minimals = set([])
        while candidates:
            candidate = candidates.pop()
            atoms = H[candidate][0]
            if atoms.isdisjoint(candidates) and atoms.isdisjoint(minimals):
                minimals.add(candidate)
        # Now G == {H[atom] : atom in minimals}
        for atom in minimals:   # 9: for g in G:
            g_atoms, g_coatoms = H[atom]
            if required_atoms is not None:
                if g_atoms.isdisjoint(required_atoms):
                    continue
            if (g_atoms, g_coatoms) in faces:
                g = faces[g_atoms, g_coatoms]
            else:               # 11: if g was newly created
                g = next_index
                faces[g_atoms, g_coatoms] = g
                next_index += 1
                Q.append((g_atoms, g_coatoms))  # 12
            L.add_edge(q, g)                    # 14

    # End of algorithm, now construct a FiniteLatticePoset.

    # In principle, it is recommended to use Poset or in this case perhaps
    # even LatticePoset, but it seems to take several times more time
    # than the above computation, makes unnecessary copies, and crashes.
    # So for now we will mimic the relevant code from Poset.

    # Enumeration of graph vertices must be a linear extension of the poset
    new_order = L.topological_sort()
    # Make sure that coatoms are in the end in proper order
    tail = [faces[atomes, frozenset([coatom])]
            for coatom, atomes in enumerate(coatom_to_atoms)]
    tail.append(faces[A, frozenset()])
    new_order = [n for n in new_order if n not in tail] + tail
    # Make sure that atoms are in the beginning in proper order
    head = [0] # We know that the empty face has index 0
    head.extend(faces[frozenset([atom]), coatoms]
                for atom, coatoms in enumerate(atom_to_coatoms)
                if required_atoms is None or atom in required_atoms)
    new_order = head + [n for n in new_order if n not in head]
    # "Invert" this list to a dictionary
    labels = {}
    for new, old in enumerate(new_order):
        labels[old] = new
    L.relabel(labels)
    # Construct the actual poset elements
    elements = [None] * next_index
    for face, index in faces.items():
        atoms, coatoms = face
        elements[labels[index]] = face_constructor(
                        tuple(sorted(atoms)), tuple(sorted(coatoms)), **kwds)
    D = {i: f for i, f in enumerate(elements)}
    L.relabel(D)
    return FiniteLatticePoset(L, elements, key=key)
