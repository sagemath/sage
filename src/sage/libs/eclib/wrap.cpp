/* sage_setup: distribution = sagemath-eclib
 */
#include <iostream>
#include <sstream>
using namespace std;

/**************** Miscellaneous functions ****************/

/* NB In eclib versions before v20190226, precision was decimal
   precision by default, while now it is bit precision. */

long mwrank_get_precision()
{
  return bit_precision();
}

void mwrank_set_precision(long n)
{
  set_bit_precision(n);
}

void mwrank_initprimes(char* pfilename, int verb)
{
  initprimes((const char*)pfilename, verb);
}


char* stringstream_to_char(ostringstream& instore) {
  int n = strlen(instore.str().data());
  char* buf = (char*)malloc(n+1);
  strcpy(buf, instore.str().data());
  return buf;
}

//////// bigint //////////

bigint* new_bigint() {
  return new bigint();
}

void del_bigint(bigint* n) {
  delete n;
}

bigint* str_to_bigint(char* s) {
  istringstream *out = new istringstream(s);
  bigint* y = new bigint();
  *out >> *y;
  delete out;
  return y;
}

char* bigint_to_str(bigint* x)
{
  ostringstream instore;
  instore << (*x);
  return stringstream_to_char(instore);
}


//////// Curvedata //////////

char* Curvedata_repr(struct Curvedata* curve)
{
  ostringstream instore;
  instore << (*curve);
  return stringstream_to_char(instore);
}

double Curvedata_silverman_bound(const Curvedata* curve)
{
  return silverman_bound(*curve);
}

double Curvedata_cps_bound(const Curvedata* curve)
{
  return cps_bound(*curve);
}

double Curvedata_height_constant(const Curvedata* curve)
{
  return height_constant(*curve);
}

char* Curvedata_getdiscr(struct Curvedata* curve)
{
  ostringstream instore;
  instore << getdiscr(*curve);
  return stringstream_to_char(instore);
}

char* Curvedata_conductor(struct Curvedata* curve)
{
  CurveRed E(*curve);
  ostringstream instore;
  instore << getconductor(E);
  return stringstream_to_char(instore);
}

char* Curvedata_isogeny_class(struct Curvedata* E, int verbose)
{
  //copied from allisog.cc
  ostringstream instore;
  CurveRed C;
  C = CurveRed(*E);
  IsogenyClass cl(C, verbose);
  cl.grow();
  vector<CurveRed> crs=cl.getcurves();
  vector<Curve> cs;
  for(unsigned int i=0; i<crs.size(); i++) cs.push_back((Curve)(crs[i]));

  //copied from point_vector_to_str
  // TODO: Would just use instore << cs, but it's buggy and printouts to stdout!
  instore << "([";
  for (unsigned int i=0; i<cs.size(); i++) {
    instore << cs[i];
    if (i+1 < cs.size())
      instore << ", ";
  }
  instore << "], ";

  instore << cl.getmatrix() << ")";
  return stringstream_to_char(instore);
}

//////// mw //////////


int mw_process(struct Curvedata* curve, struct mw* m,
                      const bigint* x, const bigint* y,
                      const bigint* z, int sat)
{
  Point P(*curve, *x, *y, *z);
  if (!P.isvalid())
    return 1;
  m->process(P, sat);
  cout<<flush;
  return 0;
}

char* point_vector_to_str(const vector<Point>& v)
{
  ostringstream instore;
  instore << "[";
  for (unsigned int i=0; i<v.size(); i++) {
    instore << v[i];
    if (i+1 < v.size())
      instore << ", ";
  }
  instore << "]";
  return stringstream_to_char(instore);
}

char* p2point_vector_to_str(const vector<P2Point>& v)
{
  ostringstream instore;
  instore << "[";
  for (unsigned int i=0; i<v.size(); i++) {
    instore << v[i];
    if (i+1 < v.size())
      instore << ", ";
  }
  instore << "]";
  return stringstream_to_char(instore);
}

char* mw_getbasis(struct mw* m)
{
  return point_vector_to_str(m->getbasis());
}

double mw_regulator(struct mw* m)
{
  return to_double(m->regulator());
}

int mw_rank(struct mw* m)
{
  return m->getrank();
}

/* Returns index and unsat long array, which user must deallocate */
int mw_saturate(struct mw* m, long* index, char** unsat,
                       long sat_bd, long sat_low_bd)
{
  vector<long> v;
  int s = m->saturate(*index, v, sat_bd, sat_low_bd);
  ostringstream instore;
  instore << v;
  *unsat  = stringstream_to_char(instore);
  return s;
}

bigfloat str_to_bigfloat(char* s) {
  istringstream *out = new istringstream(s);
  bigfloat y;
  *out >> y;
  delete out;
  return y;
}

void mw_search(struct mw* m, char* h_lim, int moduli_option, int verb)
{
  m->search(str_to_bigfloat(h_lim), moduli_option, verb);
  cout<<flush;
}


//////// two_descent //////////

long two_descent_get_rank(struct two_descent* t)
{
  return t->getrank();
}

long two_descent_get_rank_bound(struct two_descent* t)
{
  return t->getrankbound();
}

long two_descent_get_selmer_rank(struct two_descent* t)
{
  return t->getselmer();
}

char* two_descent_get_basis(struct two_descent* t)
{
  return p2point_vector_to_str(t->getbasis());
}

int two_descent_ok(const two_descent* t)
{
  return t->ok();
}

long two_descent_get_certain(const two_descent* t)
{
  return t->getcertain();
}

void two_descent_saturate(struct two_descent* t, long sat_bd, long sat_low_bd)
{
  t->saturate(sat_bd, sat_low_bd);
}

double two_descent_regulator(struct two_descent* t)
{
  return to_double(t->regulator());
}
