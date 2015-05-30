
#include <polybori/BoolePolynomial.h>
#include <polybori/BoolePolyRing.h>
#include <polybori/factories/VariableBlock.h>
#include <polybori/factories/VariableFactory.h>
#include <polybori/groebner/add_up.h>
#include <polybori/groebner/contained_variables.h>
#include <polybori/groebner/FGLMStrategy.h>
#include <polybori/groebner/groebner_alg.h>
#include <polybori/groebner/interpolate.h>
#include <polybori/groebner/linear_algebra_step.h>
#include <polybori/groebner/LiteralFactorization.h>
#include <polybori/groebner/ll_red_nf.h>
#include <polybori/groebner/nf.h>
#include <polybori/groebner/red_tail.h>
#include <polybori/groebner/minimal_elements.h>
#include <polybori/groebner/randomset.h>
#include <polybori.h>
#include <polybori/iterators/COrderedIter.h>
#include <polybori/orderings/pbori_order.h>
#include <polybori/pbori_defs.h>
#include "ccobject.h"

// M4RI
#ifdef __cplusplus
extern "C" {
#endif //__cplusplus
#define PACKED 1
#include "m4ri/m4ri.h"
#ifdef __cplusplus
}
#endif //__cplusplus

#include <sstream>
#include <vector>

USING_NAMESPACE_PBORI
USING_NAMESPACE_PBORIGB

//#define PBORI_PREALLOCATED_DEBUG

/* Allocate and Construct */
template <class T, class Arg>
T* New_p(const Arg& arg){
  return new T(arg);
}

BoolePolynomial get_ith_gen(const GroebnerStrategy& strat, int i){
    return strat.generators[i].p;
}

static int pairs_top_sugar(const GroebnerStrategy& strat){
    if (strat.pairs.pairSetEmpty())
        return -1;
    else
        return (strat.pairs.queue.top().sugar);
}

static std::vector<BoolePolynomial> someNextDegreeSpolys(GroebnerStrategy& strat, int n){
    std::vector<BoolePolynomial> res;
    assert(!(strat.pairs.pairSetEmpty()));
    strat.pairs.cleanTopByChainCriterion();
    deg_type deg=strat.pairs.queue.top().sugar;

    while((!(strat.pairs.pairSetEmpty())) && \
                (strat.pairs.queue.top().sugar<=deg) && (res.size()<n)){
        assert(strat.pairs.queue.top().sugar==deg);
        res.push_back(strat.nextSpoly());
        strat.pairs.cleanTopByChainCriterion();
    }
    return res;
}

static std::vector<Polynomial> nextDegreeSpolys(GroebnerStrategy& strat){
    std::vector<Polynomial> res;
    assert(!(strat.pairs.pairSetEmpty()));
    strat.pairs.cleanTopByChainCriterion();
    deg_type deg=strat.pairs.queue.top().sugar;

    while((!(strat.pairs.pairSetEmpty())) &&
            (strat.pairs.queue.top().sugar<=deg)){

        assert(strat.pairs.queue.top().sugar==deg);
        res.push_back(strat.nextSpoly());
        strat.pairs.cleanTopByChainCriterion();
    }
    return res;
}

static std::vector<Polynomial> small_next_degree_spolys(GroebnerStrategy& strat,
        double f, int n){
    std::vector<Polynomial> res;
    assert(!(strat.pairs.pairSetEmpty()));
    strat.pairs.cleanTopByChainCriterion();
    deg_type deg=strat.pairs.queue.top().sugar;
    wlen_type wlen=strat.pairs.queue.top().wlen;
    while((!(strat.pairs.pairSetEmpty())) &&
            (strat.pairs.queue.top().sugar<=deg) &&
            (strat.pairs.queue.top().wlen<=wlen*f+2)&& (res.size()<n)){

        assert(strat.pairs.queue.top().sugar==deg);
        res.push_back(strat.nextSpoly());
        strat.pairs.cleanTopByChainCriterion();
    }
    return res;
}

static void implications(GroebnerStrategy& strat, int i){
    strat.addNonTrivialImplicationsDelayed(strat.generators[i]);
}

inline BooleSet::const_iterator*
construct_bset_begin(void* mem, const BooleSet& bset) {
  return new(mem)  BooleSet::const_iterator(bset.begin());
}

inline BooleSet::const_iterator*
construct_bset_end(void* mem, const BooleSet& bset) {
  return new(mem)  BooleSet::const_iterator(bset.end());
}

#define PBPolyVector_set(v,i,p) v[i] = p


class ring_singleton{
 public:
  static BoolePolyRing instance() {
    static BoolePolyRing ring(1);
    return ring;
  }
};



template <class ValueType>
class DefaultRinged:
  public ValueType {
  typedef DefaultRinged self;

 public:
  typedef ValueType value_type;
  typedef value_type base;

  // Default constructor allocated memory
  DefaultRinged();

  DefaultRinged(const value_type& rhs): base(rhs) {}

  DefaultRinged(const self& rhs): base(rhs) {}

  ~DefaultRinged() { }

  self& operator=(const self& rhs) {
    return operator=(static_cast<const value_type&>(rhs));
  }

  self& operator=(const value_type& rhs) {
    base::operator=(rhs);
    return *this;
  }
};

template <class ValueType>
DefaultRinged<ValueType>::DefaultRinged():
  base(ring_singleton::instance()) { }

template <>
DefaultRinged<FGLMStrategy>::DefaultRinged():
base(ring_singleton::instance(), ring_singleton::instance(),
     PolynomialVector()) { }


template <class T>
PyObject* preallocated_to_PyString(const DefaultRinged<T> *wrapped) {
  return _to_PyString<T>(wrapped) ;
}


template <class Type>
class WrappedPtr:
  public boost::shared_ptr<Type> {
    typedef WrappedPtr self;
    typedef boost::shared_ptr<Type> base;

 public:
    WrappedPtr(): base() {}
    WrappedPtr(const self& rhs): base(rhs) {}

    template <class T1>
    WrappedPtr(const T1& arg): base(new Type(arg)) {}

    template <class T1, class T2>
    WrappedPtr(const T1& arg1, const T2& arg2): base(new Type(arg1, arg2)) {}

      template <class T1, class T2, class T3>
        WrappedPtr(const T1& arg1, const T2& arg2, const T3& arg3):
        base(new Type(arg1, arg2, arg3)) {}

    operator Type&() { return base::operator*();}
    operator const Type&() const { return base::operator*();}
};


class PBRefCounter {
 public:

  PBRefCounter(): p_count(new long(0)) {}
  PBRefCounter(const PBRefCounter& rhs):
    p_count(rhs.p_count) { ++(*p_count); }

  ~PBRefCounter() {  decrease(); }

  PBRefCounter&
  operator=(const PBRefCounter& rhs) {
    decrease();
    p_count = rhs.p_count;
    ++(*p_count);
    return *this;
  }

  bool released() {
    assert(p_count);
    if(!*p_count) {
      kill();
      return true;
    }
    return false;
  }

 private:
  void decrease() {  if(p_count) --(*p_count); }
  void kill() {
    delete p_count;
    p_count = NULL;
  }
  long* p_count;
};
