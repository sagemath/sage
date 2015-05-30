/**
 * Cython does not properly handly const constructions, hence we use
 * this rather ugly workaround to call get_sorted_learnts
 */

static CYTHON_INLINE void *sage_malloc(size_t);

uint32_t ** get_sorted_learnts_helper(CMSat::Solver* solver, uint32_t *num) {
  const CMSat::vec<CMSat::Clause *>& learnt = solver->get_sorted_learnts();
  *num = learnt.size();
  uint32_t **ret = (uint32_t**)sage_malloc(sizeof(uint32_t*)* learnt.size());
  for(size_t i=0; i<learnt.size(); i++) {
    CMSat::Clause *clause = learnt[i];
    ret[i] = (uint32_t*)sage_malloc(sizeof(uint32_t)*(clause->size()+1));
    ret[i][0] = clause->size();
    for(size_t j=0; j<clause->size(); j++) {
      ret[i][j+1] = (*clause)[j].toInt();
    }
  }
  return ret;
}


uint32_t * get_unitary_learnts_helper(CMSat::Solver* solver, uint32_t *num) {
  const CMSat::vec<CMSat::Lit> learnt = solver->get_unitary_learnts();
  *num = learnt.size();
  uint32_t *ret = (uint32_t*)sage_malloc(sizeof(uint32_t) * learnt.size());
  for(size_t i=0; i<learnt.size(); i++) {
    ret[i] = learnt[i].toInt();
  }
  return ret;
}

