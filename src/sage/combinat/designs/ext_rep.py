r"""
External representations of block designs

This module is an API to the abstract tree represented by
an XML document containing the External Representation of a list of
block designs. The module also provides the related I/O operations for
reading/writing ``ext-rep`` files or data. The parsing is based on expat.

This is a modified form of the module ``ext_rep.py`` (version 0.8)
written by Peter Dobcsanyi [Do2009]_ peter@designtheory.org.

.. TODO::

    The XML data from the designtheory.org database contains a wealth of
    information about things like automorphism groups, transitivity, cycle type
    representatives, etc, but none of this data is made available through the
    current implementation.
"""

# ***********************************************************************
# This software is released under the terms of the GNU General Public
# License, version 2 or above (your choice). For details on licensing,
# see the accompanying documentation.
#
# This is a modified form of the module ext_rep.py (version 0.8)
# written by Peter Dobcsanyi peter@designtheory.org.
#
# Copyright 2004 Peter Dobcsanyi peter@designtheory.org
#           2009 Carlo Hamalainen carlo.hamalainen@gmail.com
# ***********************************************************************

import sys
import xml.parsers.expat
import re
import os.path
import gzip
import bz2

from urllib.request import urlopen

from sage.misc.temporary_file import tmp_filename


XML_NAMESPACE = 'http://designtheory.org/xml-namespace'
DTRS_PROTOCOL = '2.0'

# The following string is the file
# http://designtheory.org/database/v-b-k/v2-b2-k2.icgsa.txt.bz2
# We use this for doctests to make sure that the parsing works.

v2_b2_k2_icgsa = \
"""<?xml version="1.0"?>
<list_of_designs
 design_type="block_design"
 dtrs_protocol="2.0"
 no_designs="1"
 pairwise_nonisomorphic="true"
 xmlns="http://designtheory.org/xml-namespace">
 <info>
  <software>
   [ DESIGN-1.1, GRAPE-4.2, GAPDoc-0.9999, GAP-4.4.3 ]
  </software>
  <software>
   [ bdstat-0.8/280, numarray-1.1.1, pydesign-0.5/274, python-2.4.0.final.0 ]
  </software>
 </info>
 <designs>
  <block_design
   b="2"
   id="v2-b2-k2-0"
   v="2">
   <blocks ordered="true">
    <block>
     <z>0</z>
     <z>1</z>
    </block>
    <block>
     <z>0</z>
     <z>1</z>
    </block>
   </blocks>
   <indicators>
    <repeated_blocks flag="true"/>
    <resolvable flag="true"/>
    <affine_resolvable
     flag="true"
     mu="2"/>
    <equireplicate
     flag="true"
     r="2"/>
    <constant_blocksize
     flag="true"
     k="2"/>
    <t_design
     flag="true"
     maximum_t="2"/>
    <connected
     flag="true"
     no_components="1"/>
    <pairwise_balanced
     flag="true"
     lambda="2"/>
    <variance_balanced flag="true"/>
    <efficiency_balanced flag="true"/>
    <cyclic flag="true"/>
    <one_rotational flag="true"/>
   </indicators>
   <combinatorial_properties>
    <point_concurrences>
     <function_on_ksubsets_of_indices
      domain_base="points"
      k="1"
      n="2"
      ordered="true"
      title="replication_numbers">
      <map>
       <preimage>
        <entire_domain/>
       </preimage>
       <image>
        <z>2</z>
       </image>
      </map>
     </function_on_ksubsets_of_indices>
     <function_on_ksubsets_of_indices
      domain_base="points"
      k="2"
      n="2"
      ordered="true"
      title="pairwise_point_concurrences">
      <map>
       <preimage>
        <entire_domain/>
       </preimage>
       <image>
        <z>2</z>
       </image>
      </map>
     </function_on_ksubsets_of_indices>
    </point_concurrences>
    <block_concurrences>
     <function_on_ksubsets_of_indices
      domain_base="blocks"
      k="1"
      n="2"
      ordered="unknown"
      title="block_sizes">
      <map>
       <preimage_cardinality>
        <z>2</z>
       </preimage_cardinality>
       <image>
        <z>2</z>
       </image>
      </map>
     </function_on_ksubsets_of_indices>
     <function_on_ksubsets_of_indices
      domain_base="blocks"
      k="2"
      n="2"
      ordered="unknown"
      title="pairwise_block_intersection_sizes">
      <map>
       <preimage_cardinality>
        <z>1</z>
       </preimage_cardinality>
       <image>
        <z>2</z>
       </image>
      </map>
     </function_on_ksubsets_of_indices>
    </block_concurrences>
    <t_design_properties>
     <parameters
      b="2"
      k="2"
      lambda="2"
      r="2"
      t="2"
      v="2"/>
     <square flag="true"/>
     <projective_plane flag="false"/>
     <affine_plane flag="false"/>
     <steiner_system flag="false"/>
     <steiner_triple_system flag="false"/>
    </t_design_properties>
    <alpha_resolvable>
     <index_flag
      flag="true"
      index="2"/>
    </alpha_resolvable>
    <t_wise_balanced>
     <index_flag
      flag="true"
      index="1"/>
     <index_flag
      flag="true"
      index="2"/>
    </t_wise_balanced>
   </combinatorial_properties>
   <automorphism_group>
    <permutation_group
     degree="2"
     domain="points"
     order="2">
     <generators>
      <permutation>
       <z>1</z>
       <z>0</z>
      </permutation>
     </generators>
     <permutation_group_properties>
      <primitive flag="true"/>
      <generously_transitive flag="true"/>
      <multiplicity_free flag="true"/>
      <stratifiable flag="true"/>
      <no_orbits value="1"/>
      <degree_transitivity value="2"/>
      <rank value="2"/>
      <cycle_type_representatives>
       <cycle_type_representative>
        <permutation>
         <z>1</z>
         <z>0</z>
        </permutation>
        <cycle_type ordered="true">
         <z>2</z>
        </cycle_type>
        <no_having_cycle_type>
         <z>1</z>
        </no_having_cycle_type>
       </cycle_type_representative>
       <cycle_type_representative>
        <permutation>
         <z>0</z>
         <z>1</z>
        </permutation>
        <cycle_type ordered="true">
         <z>1</z>
         <z>1</z>
        </cycle_type>
        <no_having_cycle_type>
         <z>1</z>
        </no_having_cycle_type>
       </cycle_type_representative>
      </cycle_type_representatives>
     </permutation_group_properties>
    </permutation_group>
    <automorphism_group_properties>
     <block_primitive flag="not_applicable"/>
     <no_block_orbits value="not_applicable"/>
     <degree_block_transitivity value="not_applicable"/>
    </automorphism_group_properties>
   </automorphism_group>
   <resolutions
    all_classes_represented="true"
    pairwise_nonisomorphic="true">
    <resolution>
     <function_on_indices
      domain="blocks"
      n="2"
      ordered="true"
      title="resolution">
      <map>
       <preimage>
        <z>0</z>
       </preimage>
       <image>
        <z>0</z>
       </image>
      </map>
      <map>
       <preimage>
        <z>0</z>
       </preimage>
       <image>
        <z>1</z>
       </image>
      </map>
     </function_on_indices>
     <automorphism_group>
      <permutation_group
       degree="2"
       domain="points"
       order="2">
       <generators>
        <permutation>
         <z>1</z>
         <z>0</z>
        </permutation>
       </generators>
      </permutation_group>
     </automorphism_group>
    </resolution>
   </resolutions>
   <statistical_properties precision="9">
    <canonical_variances
     no_distinct="1"
     ordered="true">
     <value multiplicity="1">
      <d>0.5</d>
     </value>
    </canonical_variances>
    <pairwise_variances>
     <function_on_ksubsets_of_indices
      domain_base="points"
      k="2"
      n="2"
      ordered="true">
      <map>
       <preimage>
        <entire_domain/>
       </preimage>
       <image>
        <d>1.0</d>
       </image>
      </map>
     </function_on_ksubsets_of_indices>
    </pairwise_variances>
    <optimality_criteria>
     <phi_0>
      <value>
       <d>-0.693147181</d>
      </value>
      <absolute_efficiency>
       <z>1</z>
      </absolute_efficiency>
      <calculated_efficiency>
       <z>1</z>
      </calculated_efficiency>
     </phi_0>
     <phi_1>
      <value>
       <d>0.5</d>
      </value>
      <absolute_efficiency>
       <z>1</z>
      </absolute_efficiency>
      <calculated_efficiency>
       <z>1</z>
      </calculated_efficiency>
     </phi_1>
     <phi_2>
      <value>
       <d>0.25</d>
      </value>
      <absolute_efficiency>
       <z>1</z>
      </absolute_efficiency>
      <calculated_efficiency>
       <z>1</z>
      </calculated_efficiency>
     </phi_2>
     <maximum_pairwise_variances>
      <value>
       <d>1.0</d>
      </value>
      <absolute_efficiency>
       <z>1</z>
      </absolute_efficiency>
      <calculated_efficiency>
       <z>1</z>
      </calculated_efficiency>
     </maximum_pairwise_variances>
     <E_criteria>
      <E_value index="1">
       <value>
        <d>0.5</d>
       </value>
       <absolute_efficiency>
        <z>1</z>
       </absolute_efficiency>
       <calculated_efficiency>
        <z>1</z>
       </calculated_efficiency>
      </E_value>
     </E_criteria>
    </optimality_criteria>
    <other_ordering_criteria>
     <trace_of_square_of_C>
      <value>
       <d>4.0</d>
      </value>
      <absolute_comparison>
       <z>1</z>
      </absolute_comparison>
      <calculated_comparison>
       <z>1</z>
      </calculated_comparison>
     </trace_of_square_of_C>
     <max_min_ratio_canonical_variances>
      <value>
       <d>1.0</d>
      </value>
      <absolute_comparison>
       <z>1</z>
      </absolute_comparison>
      <calculated_comparison>
       <z>1</z>
      </calculated_comparison>
     </max_min_ratio_canonical_variances>
     <max_min_ratio_pairwise_variances>
      <value>
       <d>1.0</d>
      </value>
      <absolute_comparison>
       <z>1</z>
      </absolute_comparison>
      <calculated_comparison>
       <z>1</z>
      </calculated_comparison>
     </max_min_ratio_pairwise_variances>
     <no_distinct_canonical_variances>
      <value>
       <z>1</z>
      </value>
      <absolute_comparison>
       <z>1</z>
      </absolute_comparison>
      <calculated_comparison>
       <z>1</z>
      </calculated_comparison>
     </no_distinct_canonical_variances>
     <no_distinct_pairwise_variances>
      <value>
       <z>1</z>
      </value>
      <absolute_comparison>
       <z>1</z>
      </absolute_comparison>
      <calculated_comparison>
       <z>1</z>
      </calculated_comparison>
     </no_distinct_pairwise_variances>
    </other_ordering_criteria>
    <canonical_efficiency_factors
     no_distinct="1"
     ordered="true">
     <value multiplicity="1">
      <d>1.0</d>
     </value>
    </canonical_efficiency_factors>
    <functions_of_efficiency_factors>
     <harmonic_mean alias="A">
      <value>
       <d>1.0</d>
      </value>
     </harmonic_mean>
     <geometric_mean alias="D">
      <value>
       <d>1.0</d>
      </value>
     </geometric_mean>
     <minimum alias="E">
      <value>
       <d>1.0</d>
      </value>
     </minimum>
    </functions_of_efficiency_factors>
   </statistical_properties>
  </block_design>
 </designs>
</list_of_designs>
"""


def dump_to_tmpfile(s):
    """
    Utility function to dump a string to a temporary file.

    EXAMPLES::

        sage: from sage.combinat.designs import ext_rep
        sage: file_loc = ext_rep.dump_to_tmpfile("boo")
        sage: os.remove(file_loc)
    """

    file_loc = tmp_filename()
    f = open(file_loc,"w")
    f.write(v2_b2_k2_icgsa)
    f.close()
    return file_loc


def check_dtrs_protocols(input_name, input_pv):
    """
    Check that the XML data is in a valid format. We can currently
    handle version 2.0. For more information see
    http://designtheory.org/library/extrep/

    EXAMPLES::

        sage: from sage.combinat.designs import ext_rep
        sage: ext_rep.check_dtrs_protocols('source', '2.0')
        sage: ext_rep.check_dtrs_protocols('source', '3.0')
        Traceback (most recent call last):
        ...
        RuntimeError: Incompatible dtrs_protocols: program: 2.0 source: 3.0
    """

    program_pv = DTRS_PROTOCOL
    ppv_major, ppv_minor = program_pv.split('.')
    ipv_major, ipv_minor = input_pv.split('.')
    if ppv_major != ipv_major or int(ppv_minor) < int(ipv_minor):
        msg = ('''Incompatible dtrs_protocols: program: %s %s: %s''' % (program_pv, input_name, input_pv))
        raise RuntimeError(msg)


def open_extrep_file(fname):
    """
    Try to guess the compression type from extension
    and open the extrep file.

    EXAMPLES::

        sage: from sage.combinat.designs import ext_rep
        sage: file_loc = ext_rep.dump_to_tmpfile(ext_rep.v2_b2_k2_icgsa)
        sage: proc = ext_rep.XTreeProcessor()
        sage: f = ext_rep.open_extrep_file(file_loc)
        sage: proc.parse(f)
        sage: f.close()
        sage: os.remove(file_loc)
    """

    if fname == '-':
        f = sys.stdin
    else:
        root, ext = os.path.splitext(fname)
        if ext == '.gz':
            f = gzip.GzipFile(fname)
        elif ext == '.bz2':
            f = bz2.BZ2File(fname)
        else:
            f = open(fname, 'rb')
    return f


def open_extrep_url(url):
    """
    Try to guess the compression type from extension
    and open the extrep file pointed to by the url. This function
    (unlike open_extrep_file) returns the uncompressed text contained in
    the file.

    EXAMPLES::

        sage: from sage.combinat.designs import ext_rep
        sage: file_loc = ext_rep.dump_to_tmpfile(ext_rep.v2_b2_k2_icgsa)
        sage: proc = ext_rep.XTreeProcessor()
        sage: s = ext_rep.open_extrep_url("file://" + file_loc)
        sage: proc.parse(s)
        sage: os.remove(file_loc)

        sage: from sage.combinat.designs import ext_rep
        sage: s = ext_rep.designs_from_XML_url("http://designtheory.org/database/v-b-k/v3-b6-k2.icgsa.txt.bz2") # optional - internet
    """

    f = urlopen(url)

    root, ext = os.path.splitext(url)
    if ext == '.gz':
        raise NotImplementedError
    elif ext == '.bz2':
        return bz2.decompress(f.read())
    else:
        return f.read()


pattern_integer = re.compile(r'\d+$')
pattern_decimal = re.compile(r'-?\d+\.\d+$')
pattern_rational = re.compile(r'-?\d+/\d+$')


def _encode_attribute(string):
    """
    Convert numbers in attributes into binary format.
    Currently integer and floating point conversions are implemented.

    EXAMPLES::

        sage: from sage.combinat.designs.ext_rep import _encode_attribute
        sage: _encode_attribute('1')
        1
        sage: _encode_attribute('2')
        2
        sage: _encode_attribute('12')
        12
        sage: _encode_attribute('true')
        'true'
        sage: _encode_attribute('A')
        'A'
        sage: _encode_attribute('D')
        'D'
        sage: _encode_attribute('E')
        'E'
    """
    if pattern_integer.match(string):
        return int(string)
    elif pattern_decimal.match(string):
        return float(string)
    else:
        return string


class XTree:
    '''
    A lazy class to wrap a rooted tree representing an XML document.
    The tree's nodes are tuples of the structure:
    (name, {dictionary of attributes}, [list of children])

    Methods and services of an XTree object ``t``:

    - ``t.attribute`` -- attribute named
    - ``t.child`` -- first child named
    - ``t[i]`` -- `i`-th child
    - ``for child in t:`` -- iterate over ``t``'s children
    - ``len(t)`` -- number of ``t``'s children

    If child is not an empty subtree, return the subtree as an ``XTree``
    object. If child is an empty subtree, return ``_name`` of the subtree.
    Otherwise return the child itself.

    The lazy tree idea originated from a utility class of the
    pyRXP 0.9 package by Robin Becker at ReportLab.
    '''

    def __init__(self, node):
        """
        Initialisation method given a node in an XML document.

        EXAMPLES::

            sage: from sage.combinat.designs.ext_rep import *
            sage: xt = XTree(('blocks', {'ordered': 'true'}, [('block', {}, [[0, 1, 2]]), ('block', {}, [[0, 3, 4]]), ('block', {}, [[0, 5, 6]]), ('block', {}, [[0, 7, 8]]), ('block', {}, [[0, 9, 10]]), ('block', {}, [[0, 11, 12]]), ('block', {}, [[1, 3, 5]]), ('block', {}, [[1, 4, 6]]), ('block', {}, [[1, 7, 9]]), ('block', {}, [[1, 8, 11]]), ('block', {}, [[1, 10, 12]]), ('block', {}, [[2, 3, 7]]), ('block', {}, [[2, 4, 8]]), ('block', {}, [[2, 5, 10]]), ('block', {}, [[2, 6, 12]]), ('block', {}, [[2, 9, 11]]), ('block', {}, [[3, 6, 9]]), ('block', {}, [[3, 8, 12]]), ('block', {}, [[3, 10, 11]]), ('block', {}, [[4, 5, 11]]), ('block', {}, [[4, 7, 10]]), ('block', {}, [[4, 9, 12]]), ('block', {}, [[5, 7, 12]]), ('block', {}, [[5, 8, 9]]), ('block', {}, [[6, 7, 11]]), ('block', {}, [[6, 8, 10]])]))
            sage: xt.xt_children
            [('block', {}, [[0, 1, 2]]),
             ('block', {}, [[0, 3, 4]]),
             ('block', {}, [[0, 5, 6]]),
             ('block', {}, [[0, 7, 8]]),
             ('block', {}, [[0, 9, 10]]),
             ('block', {}, [[0, 11, 12]]),
             ('block', {}, [[1, 3, 5]]),
             ('block', {}, [[1, 4, 6]]),
             ('block', {}, [[1, 7, 9]]),
             ('block', {}, [[1, 8, 11]]),
             ('block', {}, [[1, 10, 12]]),
             ('block', {}, [[2, 3, 7]]),
             ('block', {}, [[2, 4, 8]]),
             ('block', {}, [[2, 5, 10]]),
             ('block', {}, [[2, 6, 12]]),
             ('block', {}, [[2, 9, 11]]),
             ('block', {}, [[3, 6, 9]]),
             ('block', {}, [[3, 8, 12]]),
             ('block', {}, [[3, 10, 11]]),
             ('block', {}, [[4, 5, 11]]),
             ('block', {}, [[4, 7, 10]]),
             ('block', {}, [[4, 9, 12]]),
             ('block', {}, [[5, 7, 12]]),
             ('block', {}, [[5, 8, 9]]),
             ('block', {}, [[6, 7, 11]]),
             ('block', {}, [[6, 8, 10]])]
        """

        if isinstance(node, str):
            node = (node, {}, [])
        name, attributes, children = node
        self.xt_node = node
        self.xt_name = name
        self.xt_attributes = attributes
        self.xt_children = children

    def __repr__(self):
        """
        String representation of an XTree object.

        EXAMPLES::

            sage: from sage.combinat.designs.ext_rep import *
            sage: xt = XTree(('blocks', {'ordered': 'true'}, [('block', {}, [[0, 1, 2]]), ('block', {}, [[0, 3, 4]]), ('block', {}, [[0, 5, 6]]), ('block', {}, [[0, 7, 8]]), ('block', {}, [[0, 9, 10]]), ('block', {}, [[0, 11, 12]]), ('block', {}, [[1, 3, 5]]), ('block', {}, [[1, 4, 6]]), ('block', {}, [[1, 7, 9]]), ('block', {}, [[1, 8, 11]]), ('block', {}, [[1, 10, 12]]), ('block', {}, [[2, 3, 7]]), ('block', {}, [[2, 4, 8]]), ('block', {}, [[2, 5, 10]]), ('block', {}, [[2, 6, 12]]), ('block', {}, [[2, 9, 11]]), ('block', {}, [[3, 6, 9]]), ('block', {}, [[3, 8, 12]]), ('block', {}, [[3, 10, 11]]), ('block', {}, [[4, 5, 11]]), ('block', {}, [[4, 7, 10]]), ('block', {}, [[4, 9, 12]]), ('block', {}, [[5, 7, 12]]), ('block', {}, [[5, 8, 9]]), ('block', {}, [[6, 7, 11]]), ('block', {}, [[6, 8, 10]])]))
            sage: xt.__repr__()
            'XTree<blocks>'
        """

        return 'XTree<%s>' % self.xt_name

    def __getattr__(self, attr):
        """
        Return the data for the first attribute with name attr.

        EXAMPLES::

            sage: from sage.combinat.designs.ext_rep import *
            sage: xt = XTree(('blocks', {'ordered': 'true'}, [('block', {}, [[0, 1, 2]]), ('block', {}, [[0, 3, 4]]), ('block', {}, [[0, 5, 6]]), ('block', {}, [[0, 7, 8]]), ('block', {}, [[0, 9, 10]]), ('block', {}, [[0, 11, 12]]), ('block', {}, [[1, 3, 5]]), ('block', {}, [[1, 4, 6]]), ('block', {}, [[1, 7, 9]]), ('block', {}, [[1, 8, 11]]), ('block', {}, [[1, 10, 12]]), ('block', {}, [[2, 3, 7]]), ('block', {}, [[2, 4, 8]]), ('block', {}, [[2, 5, 10]]), ('block', {}, [[2, 6, 12]]), ('block', {}, [[2, 9, 11]]), ('block', {}, [[3, 6, 9]]), ('block', {}, [[3, 8, 12]]), ('block', {}, [[3, 10, 11]]), ('block', {}, [[4, 5, 11]]), ('block', {}, [[4, 7, 10]]), ('block', {}, [[4, 9, 12]]), ('block', {}, [[5, 7, 12]]), ('block', {}, [[5, 8, 9]]), ('block', {}, [[6, 7, 11]]), ('block', {}, [[6, 8, 10]])]))
            sage: xt.__getattr__('block')
            [0, 1, 2]
        """

        if attr in self.xt_attributes:
            return self.xt_attributes[attr]
        else:
            for child in self.xt_children:
                name, attributes, children = child
                if name == attr:
                    if len(attributes) > 0:
                        return XTree(child)
                    else:
                        if len(children) == 0:
                            # need this to get an empty Xtree, for append
                            return XTree(child)
                        grandchild = children[0]
                        if isinstance(grandchild, tuple):
                            if len(grandchild[1]) == 0 and \
                                len(grandchild[2]) == 0:
                                return grandchild[0]
                            else:
                                return XTree(child)
                        else:
                            return grandchild
        msg = '"%s" is not found in attributes of %s or its children.' % \
              (attr, self)
        raise AttributeError(msg)

    def __getitem__(self, i):
        """
        Get the ``i``-th item in the current node.

        EXAMPLES::

            sage: from sage.combinat.designs.ext_rep import *
            sage: xt = XTree(('blocks', {'ordered': 'true'}, [('block', {}, [[0, 1, 2]]), ('block', {}, [[0, 3, 4]]), ('block', {}, [[0, 5, 6]]), ('block', {}, [[0, 7, 8]]), ('block', {}, [[0, 9, 10]]), ('block', {}, [[0, 11, 12]]), ('block', {}, [[1, 3, 5]]), ('block', {}, [[1, 4, 6]]), ('block', {}, [[1, 7, 9]]), ('block', {}, [[1, 8, 11]]), ('block', {}, [[1, 10, 12]]), ('block', {}, [[2, 3, 7]]), ('block', {}, [[2, 4, 8]]), ('block', {}, [[2, 5, 10]]), ('block', {}, [[2, 6, 12]]), ('block', {}, [[2, 9, 11]]), ('block', {}, [[3, 6, 9]]), ('block', {}, [[3, 8, 12]]), ('block', {}, [[3, 10, 11]]), ('block', {}, [[4, 5, 11]]), ('block', {}, [[4, 7, 10]]), ('block', {}, [[4, 9, 12]]), ('block', {}, [[5, 7, 12]]), ('block', {}, [[5, 8, 9]]), ('block', {}, [[6, 7, 11]]), ('block', {}, [[6, 8, 10]])]))
            sage: xt.__getitem__(0)
            [0, 1, 2]
            sage: xt.__getitem__(1)
            [0, 3, 4]

        TESTS::

            sage: xt.__getitem__(119)
            Traceback (most recent call last):
            ...
            IndexError: XTree<blocks> has no index 119
        """
        try:
            child = self.xt_children[i]
        except IndexError:
            raise IndexError('{!r} has no index {}'.format(self, i))
        if isinstance(child, tuple):
            name, attributes, children = child
            if len(attributes) > 0:
                return XTree(child)
            else:
                grandchild = children[0]
                if isinstance(grandchild, tuple):
                    if len(grandchild[1]) == 0 and len(grandchild[2]) == 0:
                        return grandchild[0]
                    else:
                        return XTree(child)
                else:
                    return grandchild
        else:
            return child

    def __len__(self):
        """
        Return the length of the current node.

        EXAMPLES::

            sage: from sage.combinat.designs.ext_rep import *
            sage: xt = XTree(('blocks', {'ordered': 'true'}, [('block', {}, [[0, 1, 2]]), ('block', {}, [[0, 3, 4]]), ('block', {}, [[0, 5, 6]]), ('block', {}, [[0, 7, 8]]), ('block', {}, [[0, 9, 10]]), ('block', {}, [[0, 11, 12]]), ('block', {}, [[1, 3, 5]]), ('block', {}, [[1, 4, 6]]), ('block', {}, [[1, 7, 9]]), ('block', {}, [[1, 8, 11]]), ('block', {}, [[1, 10, 12]]), ('block', {}, [[2, 3, 7]]), ('block', {}, [[2, 4, 8]]), ('block', {}, [[2, 5, 10]]), ('block', {}, [[2, 6, 12]]), ('block', {}, [[2, 9, 11]]), ('block', {}, [[3, 6, 9]]), ('block', {}, [[3, 8, 12]]), ('block', {}, [[3, 10, 11]]), ('block', {}, [[4, 5, 11]]), ('block', {}, [[4, 7, 10]]), ('block', {}, [[4, 9, 12]]), ('block', {}, [[5, 7, 12]]), ('block', {}, [[5, 8, 9]]), ('block', {}, [[6, 7, 11]]), ('block', {}, [[6, 8, 10]])]))
            sage: xt.__len__()
            26
        """

        return len(self.xt_children)


class XTreeProcessor:
    '''
    An incremental event-driven parser for ext-rep documents.
    The processing stages:

    - ``<list_of_designs ...>`` opening element.
      call-back: ``list_of_designs_proc``

    - ``<list_definition>`` subtree.
      call-back: ``list_definition_proc``

    - ``<info>`` subtree.
      call-back: ``info_proc``

    - iterating over ``<designs>`` processing each ``<block_design>``
      separately.
      call-back: ``block_design_proc``

    - finishing with closing ``</designs>`` and ``</list_of_designs>``.
    '''
    def _init(self):
        """
        Internal initialisation for the processor of XTrees.

        EXAMPLES::

            sage: from sage.combinat.designs import ext_rep
            sage: proc = ext_rep.XTreeProcessor()
            sage: proc._init()
            sage: proc.current_node
            ('root0', {}, [])
        """

        self.current_node = ('root0', {}, [])
        self.node_stack = [self.current_node]
        self.in_item = False

    def __init__(self):
        """
        Internal initialisation for the processor of XTrees.

        EXAMPLES::

            sage: from sage.combinat.designs.ext_rep import *
            sage: proc = XTreeProcessor()
            sage: proc.current_node
            ('root0', {}, [])
            sage: proc.node_stack
            [('root0', {}, [])]
            sage: proc.in_item
            False
        """

        self._init()
        self.outf = sys.stdout
        # call-back handlers
        self.list_of_designs_start_proc = None
        self.list_definition_proc = None
        self.info_proc = None
        self.designs_start_proc = None
        self.block_design_proc = None
        self.designs_end_proc = None
        self.list_of_designs_end_proc = None

        self.save_designs = False
        self.list_of_designs = []

    def _start_element(self, name, attrs):
        """
        Process the start of an element with certain name and
        attributes.

        EXAMPLES::

            sage: from sage.combinat.designs.ext_rep import *
            sage: name = "block_design"
            sage: attrs = {'b': '26', 'id': 't2-v13-b26-r6-k3-L1-0', 'v': '13'}
            sage: proc = XTreeProcessor()
            sage: proc._start_element(name, attrs)
            sage: proc.current_node
            ('block_design', {'b': 26, 'id': 't2-v13-b26-r6-k3-L1-0', 'v': 13}, [])
        """

        if name == 'block_design' or name == 'info' or name == 'list_definition':
            self.in_item = True
        elif name == 'list_of_designs':
            check_dtrs_protocols('source', attrs['dtrs_protocol'])
            if self.list_of_designs_start_proc:
                self.list_of_designs_start_proc(attrs)
            # self.outf.write('<%s' % name)
            # pp_attributes(self.outf, attrs, indent='', precision_stack=[])
            # self.outf.write('>\n')
        elif name == 'designs':
            pass  # self.outf.write(' <%s>\n' % name)
        if self.in_item:
            for k, v in attrs.items():
                attrs[k] = _encode_attribute(v)
            new_node = (name, attrs, [])
            self.node_stack.append(new_node)
            self.current_node[2].append(new_node)
            self.current_node = new_node

    def _end_element(self, name):
        """
        Finish processing the element name.

        EXAMPLES::

            sage: from sage.combinat.designs.ext_rep import *
            sage: name = "block_design"
            sage: attrs = {'b': '26', 'id': 't2-v13-b26-r6-k3-L1-0', 'v': '13'}
            sage: proc = XTreeProcessor()
            sage: proc._start_element(name, attrs)
            sage: proc._end_element(name)
            sage: proc.current_node
            ('root0', {}, [])
        """

        if self.in_item:
            children = self.current_node[2]
            if len(children) > 0 and isinstance(children[0], tuple):
                if children[0][0] == 'z' or children[0][0] == 'd' \
                   or children[0][0] == 'q':
                    if children[0][0] == 'z':
                        convert = int
                    elif children[0][0] == 'd':
                        convert = float
                    else:
                        raise NotImplementedError('rational numbers')
                    ps = []
                    for x in children:
                        ps.append(convert(''.join(x[2])))
                    del children[:]
                    if name == 'block' or name == 'permutation' \
                       or name == 'preimage' or name == 'ksubset' \
                       or name == 'cycle_type' or name == 'row':
                        # these enclose lists of numbers
                        children.append(ps)
                    else:
                        # the rest is a single number
                        children.append(ps[0])
            self.node_stack.pop()
            self.current_node = self.node_stack[-1]
            if name == 'block_design':
                if self.block_design_proc:
                    self.block_design_proc(self.current_node[2][0])
                if self.save_designs:
                    init_bd = XTree(self.current_node[2][0])
                    self.list_of_designs.append((init_bd.v, list(init_bd.blocks)))
                # print_subxt(self.current_node[2][0], level=2, outf=self.outf)
                self._init()
            elif name == 'info':
                if self.info_proc:
                    self.info_proc(self.current_node[2][0])
                # print_subxt(self.current_node[2][0], level=1, outf=self.outf)
                self._init()
        else:
            if name == 'designs':
                if self.designs_end_proc:
                    self.designs_end_proc()
                # self.outf.write(' ')
            # self.outf.write('</%s>\n' % name)

    def _char_data(self, data):
        """
        Internal function to tidy up character data.

        EXAMPLES::

            sage: from sage.combinat.designs.ext_rep import *
            sage: name = "block_design"
            sage: attrs = {'b': '26', 'id': 't2-v13-b26-r6-k3-L1-0', 'v': '13'}

            sage: proc = XTreeProcessor()
            sage: proc._start_element(name, attrs)

            sage: proc._char_data(r"   [ DESIGN-1.1, GRAPE-4.2, GAPDoc-0.9999, GAP-4.4.3]")
            sage: proc.current_node
            ('block_design',
             {'b': 26, 'id': 't2-v13-b26-r6-k3-L1-0', 'v': 13},
             ['[ DESIGN-1.1, GRAPE-4.2, GAPDoc-0.9999, GAP-4.4.3]'])
        """
        if self.in_item:
            # @ this stripping may distort char data in the <info> subtree
            # if they are not bracketed in some way.
            data = data.strip()
            if data:
                # we use the xtree's children list here to collect char data
                # since only leaves have char data.
                self.current_node[2].append(data)

    def parse(self, xml_source):
        """
        The main parsing function. Given an XML source (either a file
        handle or a string), parse the entire XML source.

        EXAMPLES::

            sage: from sage.combinat.designs import ext_rep
            sage: file_loc = ext_rep.dump_to_tmpfile(ext_rep.v2_b2_k2_icgsa)
            sage: proc = ext_rep.XTreeProcessor()
            sage: proc.save_designs = True
            sage: f = ext_rep.open_extrep_file(file_loc)
            sage: proc.parse(f)
            sage: f.close()
            sage: os.remove(file_loc)
            sage: proc.list_of_designs[0]
            (2, [[0, 1], [0, 1]])
        """

        p = xml.parsers.expat.ParserCreate()
        p.StartElementHandler = self._start_element
        p.EndElementHandler = self._end_element
        p.CharacterDataHandler = self._char_data

        if isinstance(xml_source, (str, bytes)):
            p.Parse(xml_source)
        else:
            p.ParseFile(xml_source)


def designs_from_XML(fname):
    """
    Return a list of designs contained in an XML file fname. The list
    contains tuples of the form (v, bs) where v is the number of points of
    the design and bs is the list of blocks.

    EXAMPLES::

        sage: from sage.combinat.designs import ext_rep
        sage: file_loc = ext_rep.dump_to_tmpfile(ext_rep.v2_b2_k2_icgsa)
        sage: ext_rep.designs_from_XML(file_loc)[0]
        (2, [[0, 1], [0, 1]])
        sage: os.remove(file_loc)

        sage: from sage.combinat.designs import ext_rep
        sage: from sage.combinat.designs.block_design import BlockDesign
        sage: file_loc = ext_rep.dump_to_tmpfile(ext_rep.v2_b2_k2_icgsa)
        sage: v, blocks = ext_rep.designs_from_XML(file_loc)[0]
        sage: d = BlockDesign(v, blocks)
        sage: d.blocks()
        [[0, 1], [0, 1]]
        sage: d.is_t_design(t=2)
        True
        sage: d.is_t_design(return_parameters=True)
        (True, (2, 2, 2, 2))
    """

    proc = XTreeProcessor()
    proc.save_designs = True
    f = open_extrep_file(fname)
    proc.parse(f)
    f.close()

    return proc.list_of_designs


def designs_from_XML_url(url):
    """
    Return a list of designs contained in an XML file named by a URL.
    The list contains tuples of the form (v, bs) where v is the number
    of points of the design and bs is the list of blocks.

    EXAMPLES::

        sage: from sage.combinat.designs import ext_rep
        sage: file_loc = ext_rep.dump_to_tmpfile(ext_rep.v2_b2_k2_icgsa)
        sage: ext_rep.designs_from_XML_url("file://" + file_loc)[0]
        (2, [[0, 1], [0, 1]])
        sage: os.remove(file_loc)

        sage: from sage.combinat.designs import ext_rep
        sage: ext_rep.designs_from_XML_url("http://designtheory.org/database/v-b-k/v3-b6-k2.icgsa.txt.bz2") # optional - internet
        [(3, [[0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 2]]),
         (3, [[0, 1], [0, 1], [0, 1], [0, 1], [0, 2], [0, 2]]),
         (3, [[0, 1], [0, 1], [0, 1], [0, 1], [0, 2], [1, 2]]),
         (3, [[0, 1], [0, 1], [0, 1], [0, 2], [0, 2], [0, 2]]),
         (3, [[0, 1], [0, 1], [0, 1], [0, 2], [0, 2], [1, 2]]),
         (3, [[0, 1], [0, 1], [0, 2], [0, 2], [1, 2], [1, 2]])]
    """
    proc = XTreeProcessor()
    proc.save_designs = True
    s = open_extrep_url(url)
    proc.parse(s)

    return proc.list_of_designs
