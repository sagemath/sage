/* sage_setup: distribution = sagemath-brial
 */
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

#include <sstream>
#include <vector>

USING_NAMESPACE_PBORI
USING_NAMESPACE_PBORIGB


static int pairs_top_sugar(const GroebnerStrategy& strat){
    if (strat.pairs.pairSetEmpty())
        return -1;
    else
        return (strat.pairs.queue.top().sugar);
}

static std::vector<BoolePolynomial> someNextDegreeSpolys(GroebnerStrategy& strat, size_t n){
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
        double f, size_t n){
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
