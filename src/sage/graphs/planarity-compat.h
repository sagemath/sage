#include <planarity/graphLib.h>

#if GP_PROJECTVERSION_MAJOR < 5
#define gp_EnsureVertexCapacity gp_InitGraph
#endif
