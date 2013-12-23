lib_LTLIBRARIES       = libTOPCOM.la

libTOPCOM_la_SOURCES = \
Admissibles.cc \
CheckTriang.cc \
Circuits.cc \
Cocircuits.cc \
CommandlineOptions.cc \
CompressedIntegerSet.cc \
ComputeTriangs.cc \
Facets.cc \
FastSimplicialComplex.cc \
Field.cc \
FineTriang.cc \
Flip.cc \
IndexTableArray.cc \
IntegerSet.cc \
InteriorFacets.cc \
MarkedFlips.cc \
Matrix.cc \
Permutation.cc \
PlacingTriang.cc \
PointConfiguration.cc \
RealChiro.cc \
SimplicialComplex.cc \
SimplicialComplexTemplate.cc \
SparseIntegerSet.cc \
SparseSimplicialComplex.cc \
StairCaseMatrix.cc \
SymmetricBFS.cc \
Symmetry.cc \
TriangFlips.cc \
TriangNode.cc \
Vector.cc \
VertexFacetTable.cc \
VertexFacetTableArray.cc \
VirtualChiro.cc

include_HEADERS = \
Admissibles.hh \
Array.hh \
CheckTriang.hh \
Chirotope.hh \
Circuits.hh \
Cocircuits.hh \
CommandlineOptions.hh \
CompressedIntegerSet.hh \
ComputeTriangs.hh \
Facets.hh \
FastSimplicialComplex.hh \
Field.hh \
FineTriang.hh \
Flip.hh \
Global.hh \
HashKey.hh \
HashMap.hh \
HashSet.hh \
IndexTableArray.hh \
IntegerSet.hh \
InteriorFacets.hh \
MarkedFlips.hh \
Matrix.hh \
Pair.hh \
PartialTriang.hh \
Permutation.hh \
PlacingTriang.hh \
PlainArray.hh \
PlainHashIndexTable.hh \
PlainHashMap.hh \
PlainHashSet.hh \
PlainHashTable.hh \
PointConfiguration.hh \
RealChiro.hh \
RefCount.hh \
SimplicialComplex.hh \
SimplicialComplexTemplate.hh \
SimplicialComplexTemplate_Declarations.hh \
SimplicialComplexTemplate_Definitions.hh \
SparseIntegerSet.hh \
SparseSimplicialComplex.hh \
StairCaseMatrix.hh \
SymmetricBFS.hh \
Symmetry.hh \
TriangFlips.hh \
TriangNode.hh \
Vector.hh \
VertexFacetTable.hh \
VertexFacetTableArray.hh \
VirtualChiro.hh


libTOPCOM_la_CPPFLAGS = \
	-DVERBOSE -DGMPRATIONAL -ftemplate-depth-35 \
	$(GMP_CFLAGS) $(CDD_CFLAGS) \
	-I../lib-src-reg -I../wrap-gmp-gmpxx 

libTOPCOM_la_LIBADD = $(GMP_LIBS) $(CDD_LIBS)

