#include <algorithm>
#include "triangulations.h"


//------------ traverse the edges of the secondary polytope ------------
triangulations::triangulations(const flips& all_flips)
  : hash_max(10000),
    hash_list(hash_max,hash_max),
    bistellar_flips(all_flips),
    position(0),
    star(-1),
    fine(false),
    need_resize(false)
{}


// find the correct position for t in the hash
void triangulations::find_hash_position(const compact_simplices& t,
                                        hash_value& pos, bool& is_new) const
{
  const hash_value initial_guess = t.hash_function() % hash_max;

  for (hash_value i=0; i<hash_max; ++i) {
    pos = (initial_guess+i) % hash_max;
    if (hash_list[pos]==hash_max) {
      // found empty place in hash
      is_new=true;
      if (i>5)
        need_resize = true;
      return;
    }
    else if ((*this)[hash_list[pos]]==t) {
      // found ourselves in hash
      is_new = false;
      return;
    }
    // hash collision: continue
  }
  assert(false);  // the need_resize was not honored, bug?
}


void triangulations::add_triang_if_new(const compact_simplices & new_triang)
{
  hash_value pos;
  bool is_new;
  find_hash_position(new_triang, pos, is_new);
  if (!is_new) return;

  while (need_resize) {
    hash_max = 2*hash_max+1;
    hash_list = std::vector<size_t>(hash_max, hash_max);
    for (size_t i=0; i<size(); i++) {
      find_hash_position( (*this)[i], pos, is_new);
      assert(is_new);
      hash_list[pos] = i;
    }
    need_resize = false;
    find_hash_position(new_triang, pos, is_new);
  }

  push_back(new_triang);
  hash_list[pos] = size()-1;
}


void triangulations::add_neighbours(const simplices & s)
{
  for (flips::const_iterator
         f=bistellar_flips.begin(); f!=bistellar_flips.end(); ++f) {
    goodcircuit goody(s,*f);
    if (goody.is_good()) {
      goody.do_flip(s,*f);
      compact_simplices new_triang=goody.get_neighbor();
      add_triang_if_new(new_triang);
    }
  }
}


bool triangulations::have_more_triangulations()
{
  while (position != this->size()) {
     // eat all non-star triangulations
    simplices triangulation((*this)[position]);

    if ((star>=0) && !(triangulation.starshaped(star))) {
      next_triangulation();
      continue;
    }
    // eat non-fine triangulations
    if (fine && !(triangulation.fine())) {
      next_triangulation();
      continue;
    }
    return true;
  }
  return false;
}


const compact_simplices& triangulations::next_triangulation()
{
  add_neighbours((*this)[position]);
  return (*this)[position++];
}




//------------ to be called from Python -------------------------
triangulations_ptr init_triangulations
(int n, int d, int star, bool fine, PyObject* py_seed, PyObject* py_flips)
{
  vertices().set_dimensions(n,d);

  compact_simplices seed;
  for (int i=0; i<PySequence_Size(py_seed); i++) {
    PyObject* simplex = PySequence_GetItem(py_seed,i);
    seed.push_back(PyLong_AsLong(simplex));
    Py_DECREF(simplex);
  }

  flips all_flips;
  for (int i=0; i<PySequence_Size(py_flips); i++) {
    PyObject* py_flip = PySequence_GetItem(py_flips,i);
    PyObject* py_flip_pos = PySequence_GetItem(py_flip,0);
    PyObject* py_flip_neg = PySequence_GetItem(py_flip,1);

    std::vector<vertices> pos;
    for (int j=0; j<PySequence_Size(py_flip_pos); j++) {
      PyObject* py_simplex = PySequence_GetItem(py_flip_pos,j);
      vertices simplex;
      for (int k=0; k<PySequence_Size(py_simplex); k++) {
        PyObject* py_vertex = PySequence_GetItem(py_simplex,k);
        simplex.insert(simplex.begin(), PyLong_AsLong(py_vertex));
        Py_DECREF(py_vertex);
      }
      pos.push_back(simplex);
      Py_DECREF(py_simplex);
    }

    std::vector<vertices> neg;
    for (int j=0; j<PySequence_Size(py_flip_neg); j++) {
      PyObject* py_simplex = PySequence_GetItem(py_flip_neg,j);
      vertices simplex;
      for (int k=0; k<PySequence_Size(py_simplex); k++) {
        PyObject* py_vertex = PySequence_GetItem(py_simplex,k);
        simplex.insert(simplex.begin(), PyLong_AsLong(py_vertex));
        Py_DECREF(py_vertex);
      }
      neg.push_back(simplex);
      Py_DECREF(py_simplex);
    }

    all_flips.push_back(flip(pos,neg));
    all_flips.push_back(flip(neg,pos));
    Py_DECREF(py_flip_pos);
    Py_DECREF(py_flip_neg);
    Py_DECREF(py_flip);
  }

  // print data so we can see that it worked
  /*
  std::cout << "\n" << "Seed triangulation: " << seed << "\n";
  {
    std::cout << "Vertices of seed triangulation: ";
    simplices s(seed);
    std::cout << s << "\n";
  }
  std::cout << "All flips:\n";
  for (flips::const_iterator i=all_flips.begin(); i!=all_flips.end(); ++i)
    std::cout << *i << std::endl;
  */

  triangulations_ptr t = new triangulations(all_flips);

  if (star>=0)
    t->require_star_triangulation(star);

  if (fine)
    t->require_fine_triangulation(fine);

  t->add_triang_if_new(seed);

  return t;
}


PyObject* next_triangulation(triangulations_ptr t)
{
  if (!t->have_more_triangulations())
    return PyTuple_New(0);

  const compact_simplices& triang = t->next_triangulation();
  PyObject* py_triang = PyTuple_New(triang.size());
  for (size_t i=0; i<triang.size(); i++)
    PyTuple_SET_ITEM(py_triang, i, PyLong_FromLong(triang[i]));

  return py_triang;
}


void delete_triangulations(triangulations_ptr t)
{
  delete(t);
}

