#ifndef __TRIANGULATION_DATA_H__
#define __TRIANGULATION_DATA_H__

#include <cstdlib>
#include <time.h>
#include <iostream>
#include <set>
#include <functional>
#include <vector>


// may be insufficient if you have 16 bit int's
typedef unsigned int hash_value;



// typedef unsigned long int simplex;
typedef int vertex;
typedef int simplex;

class flip;
class flips;
class vertices_lookup;


// a set of vertices; usually vertices of a simplex
class vertices: public std::set<vertex,std::less<vertex> >
{
 private:
  static int n;
  static int d;
  static vertices_lookup lookup;
 protected:
 public:
  vertices();
  vertices(const simplex&);
  const simplex vertices_to_simplex() const;
  const vertices & simplex_to_vertices(const simplex &); // uses lookup-table
  void set_dimensions(int N, int D);
  friend std::ostream & operator << (std::ostream &, const vertices &);
  friend bool operator==(const vertices &, const vertices &);
  bool full_set() const { return this->size() == n; }
};

// tables used by vertices to make computation faster
class vertices_lookup
{
 private:
  int n,d;
  std::vector<vertices> SimplexToVertices;
  std::vector<std::vector<int> > fast_binomial;
 protected:
  vertices manual_vertices_to_simplex(const simplex &) const;
 public:
  void generate_tables(int N, int D);
  const vertices & simplex_to_vertices(const simplex &) const;
  int get_binomial(int i, int j) const;
};


// total order on the vertices
struct vertices_order
{
  bool operator() (const vertices &, const vertices &) const;
};


// simplicial complex, data compressed by vertices_to_simplex
class compact_simplices: public std::vector<simplex>
{
 private:
 protected:
 public:
  compact_simplices();
  virtual ~compact_simplices();
  hash_value hash_function() const;
  friend std::ostream & operator <<
    (std::ostream &, const compact_simplices &);
  friend const bool operator==
    (const compact_simplices &, const compact_simplices &);
};

// same but with complete vertex information
class simplices: public compact_simplices
{
 private:
  std::vector<vertices> v;
 protected:
  void compress();
  void decompress();
 public:
  simplices() { };
  simplices(const std::set<vertices,vertices_order> &);
  simplices(const compact_simplices &);
  bool starshaped(const vertex) const;
  bool fine() const;
  virtual ~simplices();
  const std::vector<vertices> & get_vertices() const { return(v); }
  friend std::ostream & operator << (std::ostream & out, const simplices & s);
};

class flip
{
 private:
  std::vector<vertices> deltaplus, deltaminus;
 protected:
 public:
  flip();
  flip(const std::vector<vertices>& pos, const std::vector<vertices>& neg);
  flip(const flip &);
  virtual ~flip();
  flip & operator =(const flip &);
  friend std::ostream & operator <<(std::ostream &, const flip &);
  const std::vector<vertices> & get_deltaplus()  const { return deltaplus; };
  const std::vector<vertices> & get_deltaminus() const { return deltaminus;};
  void mirror();
};

class flips: public std::vector<flip>
{
 protected:
  void sort_flips_containing();
 public:
  flips();
  virtual ~flips();
};


class goodcircuit
{
 private:
  std::vector<std::vector<vertices> > supported;
  flip supporter;
  std::vector<std::set<vertices,vertices_order> > link;
  std::set<vertices,vertices_order> bistellarneighbor;
  bool good;
 public:
  goodcircuit(const simplices &, const flip &);
  bool is_good() const { return(good); }
  void do_flip(const simplices &, const flip &);
  simplices get_neighbor() const { return(bistellarneighbor); }
};


#endif
