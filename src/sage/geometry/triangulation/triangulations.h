#ifndef TRIANGULATION__H
#define TRIANGULATION__H

#include "data.h"
#include <Python.h>


class triangulations: public std::vector<compact_simplices>
{
private:
  hash_value hash_max;
  compact_simplices no_triangulation_instance;
  compact_simplices::const_iterator no_triangulation;
  std::vector<size_t> hash_list;
  flips bistellar_flips;
  int position;
  int star;
  bool fine;
  mutable bool need_resize;
protected:
  void find_hash_position(const compact_simplices&, hash_value&, bool&) const;
  void add_triangulation(const compact_simplices &);
public:
  triangulations(const flips&);

  void add_triang_if_new(const compact_simplices &);
  void add_neighbours(const simplices &);

  void require_star_triangulation(const int s=-1) { star=s; };
  void require_fine_triangulation(const bool f=true) { fine=f; };

  bool have_more_triangulations();
  const compact_simplices& next_triangulation();
};


typedef triangulations* triangulations_ptr;

triangulations_ptr init_triangulations
(int n, int d, int star, bool fine, PyObject* py_seed, PyObject* py_flips);

PyObject* next_triangulation(triangulations_ptr);

void delete_triangulations(triangulations_ptr);


#endif
