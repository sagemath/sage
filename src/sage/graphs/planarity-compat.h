#ifndef SAGE_GRAPHS_PLANARITY_COMPAT_H
#define SAGE_GRAPHS_PLANARITY_COMPAT_H

/*
 * Compatibility shim across releases of the edge-addition-planarity-suite.
 *
 * Header location:
 *   - planarity 3 installs only <planarity/graph.h>;
 *   - planarity 4 installs <planarity/graphLib.h> and a <planarity/graph.h> shim;
 *   - planarity 5 installs only <planarity/graphLib.h>.
 * graphLib.h is therefore preferred whenever it is present.
 *
 * Function name:
 *   - up to planarity 4 the vertex-capacity routine is gp_InitGraph();
 *   - planarity 5 renamed it to gp_EnsureVertexCapacity().
 * planarity.pyx always calls gp_EnsureVertexCapacity(); for older releases it is
 * aliased back to gp_InitGraph() below.
 *
 * meson.build probes for the headers and defines SAGE_PLANARITY_HAVE_GRAPHLIB_H
 * or SAGE_PLANARITY_HAVE_GRAPH_H; the __has_include path is only a fallback for
 * compiling this header outside the meson build.
 */
#if defined(SAGE_PLANARITY_HAVE_GRAPHLIB_H)
#include <planarity/graphLib.h>
#elif defined(SAGE_PLANARITY_HAVE_GRAPH_H)
#include <planarity/graph.h>
#define SAGE_PLANARITY_NEEDS_INITGRAPH_ALIAS
#elif defined(__has_include)
#if __has_include(<planarity/graphLib.h>)
#include <planarity/graphLib.h>
#elif __has_include(<planarity/graph.h>)
#include <planarity/graph.h>
#define SAGE_PLANARITY_NEEDS_INITGRAPH_ALIAS
#else
#error "neither planarity/graphLib.h nor planarity/graph.h is available"
#endif
#else
#include <planarity/graphLib.h>
#endif

/*
 * Alias gp_EnsureVertexCapacity() to gp_InitGraph() for releases predating the
 * rename: the graph.h path (planarity 3, which does not define
 * GP_PROJECTVERSION_MAJOR) and planarity 4 (GP_PROJECTVERSION_MAJOR == 4).
 */
#if defined(SAGE_PLANARITY_NEEDS_INITGRAPH_ALIAS) || \
    (defined(GP_PROJECTVERSION_MAJOR) && GP_PROJECTVERSION_MAJOR < 5)
#define gp_EnsureVertexCapacity gp_InitGraph
#endif

#endif
