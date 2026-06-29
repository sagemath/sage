/*
 * Include the planarity public header in a way that works across the header
 * layouts shipped by different planarity versions and distributions:
 *
 * - planarity 5 (and the SageMath spkg) install the umbrella header as
 *   <planarity/graphLib.h>;
 * - some packages put their include path inside the planarity directory, so
 *   the umbrella header is reached as <graphLib.h>;
 * - older planarity (e.g. the libplanarity-dev currently shipped by Debian and
 *   Ubuntu) provides <planarity/graph.h> but no graphLib.h.
 *
 * Sage used <planarity/graph.h> until gh-42405; switching unconditionally to
 * <planarity/graphLib.h> there broke the build on the last layout above.
 */
#if !defined(__has_include)
#  include <planarity/graphLib.h>
#elif __has_include(<planarity/graphLib.h>)
#  include <planarity/graphLib.h>
#elif __has_include(<graphLib.h>)
#  include <graphLib.h>
#else
#  include <planarity/graph.h>
#endif

#if GP_PROJECTVERSION_MAJOR < 5
#define gp_EnsureVertexCapacity gp_InitGraph
#endif
