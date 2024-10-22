#include <bliss/graph.hh>
#include <bliss/digraph.hh>

inline void bliss_find_automorphisms(bliss::Graph *graph, void (*hook)(void *user_param, unsigned int n, const unsigned int *aut), void *hook_user_param, bliss::Stats s)
{
  auto report_aut = [&](unsigned int n, const unsigned int *aut) -> void {
    if(hook)
      (*hook)(hook_user_param, n, aut);
  };

  graph->find_automorphisms(s, report_aut);
}

inline void bliss_find_automorphisms(bliss::Digraph *graph, void (*hook)(void *user_param, unsigned int n, const unsigned int *aut), void *hook_user_param, bliss::Stats s)
{
  auto report_aut = [&](unsigned int n, const unsigned int *aut) -> void {
    if(hook)
      (*hook)(hook_user_param, n, aut);
  };

  graph->find_automorphisms(s, report_aut);
}
