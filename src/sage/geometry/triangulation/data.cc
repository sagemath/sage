#include <iostream>
#include <fstream>
#include <set>
#include <functional>
#include <vector>
#include <algorithm>
#include "functions.h"
#include "data.h"




// ------- class vertices ----------------------------

int vertices::n=0;

int vertices::d=0;

vertices_lookup vertices::lookup=vertices_lookup();


void vertices::set_dimensions(int N, int D)
{
  n=N; d=D;
  lookup.generate_tables(n,d);
}


vertices::vertices()
{
}


vertices::vertices(const simplex& S)
{
  simplex_to_vertices(S);
}


const simplex vertices::vertices_to_simplex() const
{
  simplex s=1;
  vertex k=1, l;
  const_iterator l_iterator=begin();
  for (vertex i=1; i<=d; ++i) {
    l = (*l_iterator)+1;
    for (vertex j=k; j<=l-1; ++j)
      s=s+lookup.get_binomial(n-j,d-i);
    k=l+1;
    ++l_iterator;
  }
  return(s);
}


const vertices & vertices::simplex_to_vertices(const simplex & S)
{
  (*this) = lookup.simplex_to_vertices(S);
  return(*this);
}


std::ostream & operator << (std::ostream & out, const vertices & v)
{
  out << *( v.begin() );
  vertices::const_iterator i=v.begin();   i++;
  for (; i!=v.end(); ++i)
    out << "_" << *i;
  return(out);
}


bool vertices_order::operator()(const vertices & a, const vertices & b) const
{
  size_t n1=a.size(), n2=b.size();
  if (n1<n2) return(true);
  if (n1>n2) return(false);
  // a,b have same size
  for (vertices::iterator i=a.begin(), j=b.begin();
       i!=a.end() && j!=b.end(); ++i, ++j) {
    if (*i < *j) return(true);
    if (*i > *j) return(false);
  }
  // a == b
  return(false);
}


bool operator==(const vertices & a, const vertices & b)
{
  return( std::set<vertex,std::less<vertex> >(a) ==
          std::set<vertex,std::less<vertex> >(b) );
}


// ------- class vertices_lookup --------------------

const vertices & vertices_lookup::simplex_to_vertices(const simplex & S) const
{
  return( SimplexToVertices[S-1] );
}


// this one is messy, but it only reverses vertices_to_simplex()
vertices vertices_lookup::manual_vertices_to_simplex(const simplex & S) const
{
  vertices result;
  simplex s=S;
  simplex b;
  vertex i,j,l=0,k;
  for (k=1; k<d; k++) {
    l++;  i=l;   j=1;
    b=binomial(n-l,d-k);
    while (s>b && b>0) {
      j++;  l++;
      s=s-b;
      b=binomial(n-l,d-k);
    }
    result.insert(result.begin(),l-1);
  }
  result.insert(result.begin(),s+l-1);
  return(result);
}

// tweaked for special (e.g. negative) values
int vertices_lookup::get_binomial(int i, int j) const
{
  if (i>=0 && i<=n && j>=0 && j<=d && j<=i)
    return(fast_binomial[i][j]);
  else
    return(1);
}


void vertices_lookup::generate_tables(int N, int D)
{
  n=N;  d=D;
  fast_binomial.erase(fast_binomial.begin(),fast_binomial.end());
  for(int i=0; i<=n; ++i) {
    std::vector<int> row;
    for(int j=0; j<=i && j<=d; j++)
      row.push_back( binomial(i,j) );
    fast_binomial.push_back( row );
  }
  SimplexToVertices.erase(SimplexToVertices.begin(),SimplexToVertices.end());
  for(int s=1; s<=binomial(n,d); ++s)
    SimplexToVertices.push_back(manual_vertices_to_simplex(s));
}




// ------- class compact_simplices -------------------

compact_simplices::compact_simplices()
{
  reserve(50);
}

compact_simplices::~compact_simplices()
{
}

std::ostream & operator << (std::ostream & out, const compact_simplices & t)
{
  out << "[" << *( t.begin() );
  for (compact_simplices::const_iterator i=t.begin()+1; i!=t.end(); ++i)
    out << "," << *i;
  out << "]";
  return(out);
}

const bool operator==(const compact_simplices & a, const compact_simplices & b)
{
  // friendship is not inherited
  return( std::vector<simplex>(a) == std::vector<simplex>(b)  );
}

hash_value compact_simplices::hash_function() const
{
  hash_value result=0;
  int n=101;
  for (const_iterator i=begin(); i!=end(); ++i, n+=37)
    result += n* (*i) + n*n;
  return(result);
}





// ------- class simplices ---------------------------

simplices::simplices(const compact_simplices & s)
  :compact_simplices(s)
{
  decompress();    // de-compress vertex information
}

simplices::simplices(const std::set<vertices,vertices_order> & s)
{
  v.erase(v.begin(), v.end());
  copy(s.begin(), s.end(), back_inserter(v) );
  compress();
}


simplices::~simplices()
{
}

void simplices::compress()
{
  erase(begin(),end());
  for (std::vector<vertices>::const_iterator i=v.begin(); i!=v.end(); i++)
    push_back( (*i).vertices_to_simplex() );
  std::sort(begin(),end());
}


void simplices::decompress()
{
  v.erase(v.begin(),v.end());
  for (const_iterator i=begin(); i!=end(); i++)
    v.push_back( vertices().simplex_to_vertices(*i) );
}


bool simplices::starshaped(const vertex origin) const
{
  bool result = true;
  for (std::vector<vertices>::const_iterator
         vi=v.begin(); vi!=v.end(); vi++) {
    result=result && (find(vi->begin(),vi->end(),origin)!=vi->end());
  }
  return(result);
}


bool simplices::fine() const
{
  vertices support;
  for (std::vector<vertices>::const_iterator
         vi=v.begin(); vi!=v.end(); vi++) {
    support.insert(vi->begin(), vi->end());
  }
  return support.full_set();
}


std::ostream & operator << (std::ostream & out, const simplices & s)
{
  out << "[" << s.v.front();
  for (std::vector<vertices>::const_iterator
         si=s.v.begin()+1; si!=s.v.end(); si++)
    out << ", " << *si;
  out << "]";
  return(out);
}

// ------- class flip --------------------------------

flip::flip()
{
  deltaplus.reserve(10);
  deltaminus.reserve(10);
}

flip::flip(const std::vector<vertices>& pos, const std::vector<vertices>& neg)
  : deltaplus(pos), deltaminus(neg)
{}

flip::flip(const flip & copyfrom)
{
  deltaplus.reserve(10);
  deltaminus.reserve(10);
  (*this)=copyfrom;
}

flip::~flip()
{
}

flip & flip::operator =(const flip & copyfrom)
{
  deltaplus  = copyfrom.get_deltaplus();
  deltaminus = copyfrom.get_deltaminus();
  return(*this);
}

std::ostream & operator << (std::ostream & out, const flip & f)
{
  out << "< ";
  for (std::vector<vertices>::const_iterator i=f.get_deltaplus().begin();
       i!=f.get_deltaplus().end(); ++i) {
    std::cout << *i << " ";
  }
  out << "| ";
  for (std::vector<vertices>::const_iterator i=f.get_deltaminus().begin();
       i!=f.get_deltaminus().end(); ++i) {
    std::cout << *i << " ";
  }
  out << ">";
  return(out);
}

void flip::mirror()
{
  swap(deltaplus,deltaminus);
}




// ------- class flips -------------------------------

flips::flips()
{
}

flips::~flips()
{
}


// ------- class goodcircuit -------------------------


goodcircuit::goodcircuit(const simplices & s, const flip & f)
  :supporter(f)
{
  const std::vector<vertices> & deltaplus = f.get_deltaplus();
  const std::vector<vertices> & v = s.get_vertices();
  supported.reserve(deltaplus.size());
  good=true;


  // find simplices that include members of deltaplus
  std::vector<vertices> empty_result;
  for (size_t n=0; n<deltaplus.size(); n++) {
    supported.push_back(empty_result);  supported[n].reserve(10);
    bool found_match=false;
    for (std::vector<vertices>::const_iterator i=v.begin(); i!=v.end(); ++i) {
      if (includes( (*i).begin(), (*i).end(),
                    deltaplus[n].begin(), deltaplus[n].end() )) {
        found_match=true;
        supported[n].push_back(*i);
      }
    }
    good=good && found_match;
    if (!good) return;
  }

  // check that all links are equal
  for (size_t n=0; n<deltaplus.size(); n++) {
    std::set<vertices,vertices_order> l;
    for (std::vector<vertices>::const_iterator i=supported[n].begin();
             i!=supported[n].end(); ++i) {
      vertices one_simplex;
      std::insert_iterator<vertices>
        ins(one_simplex,one_simplex.begin());
      set_difference( (*i).begin(), (*i).end(),
                       deltaplus[n].begin(), deltaplus[n].end(),
                       ins );
      l.insert(l.begin(),one_simplex);
    }
    link.push_back(l);
    if (n>0 && link[0]!=link[n]) {
      good=false;
      return;
    }
  }
}


void goodcircuit::do_flip(const simplices & s, const flip & f)
{
  bistellarneighbor.erase(bistellarneighbor.begin(), bistellarneighbor.end());

  const std::set<vertices,vertices_order> & use_link = link[0];
  const std::vector<vertices> & v = s.get_vertices();
  const std::vector<vertices> & deltaplus  = f.get_deltaplus();
  const std::vector<vertices> & deltaminus = f.get_deltaminus();

  // start with all simplices
  std::set<vertices,vertices_order> tri;
  copy(v.begin(), v.end(), inserter(tri,tri.begin()) );

  // remove all supported simplices
  std::set<vertices,vertices_order> to_remove;

  for (std::set<vertices,vertices_order>::const_iterator i=use_link.begin();
       i!=use_link.end(); ++i)
    for (std::vector<vertices>::const_iterator j=deltaplus.begin();
         j!=deltaplus.end(); ++j) {
      vertices one_simplex;
      set_union( (*j).begin(), (*j).end(),
                 (*i).begin(), (*i).end(),
                 inserter(one_simplex,one_simplex.begin()) );
      to_remove.insert(to_remove.begin(),one_simplex);
    }
  set_difference( tri.begin(), tri.end(),
                  to_remove.begin(), to_remove.end(),
                  inserter(bistellarneighbor,bistellarneighbor.begin()),
                  vertices_order() );

  // now add the flip of the supported simplices
  for (std::set<vertices,vertices_order>::const_iterator i=use_link.begin();
       i!=use_link.end(); ++i)
    for (std::vector<vertices>::const_iterator j=deltaminus.begin();
         j!=deltaminus.end(); j++) {
      vertices one_simplex;
      set_union( (*j).begin(), (*j).end(),
                 (*i).begin(), (*i).end(),
                 inserter(one_simplex,one_simplex.begin()) );
      bistellarneighbor.insert( bistellarneighbor.begin(), one_simplex );
    }
}













