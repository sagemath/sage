// cmatcher.cpp: commutative matching
// Distributed under GPL2+
// Author: 2018: Ralf Stephan <ralf@ark.in-berlin.de>
//
// This differs from GiNaC's original basic::match() because commutative
// structures need a special algorithm. We follow the outline in
// Manuel Krebber's Master Thesis, "Non-linear Associative-Commutative
// Many-to-One Pattern Matching with Sequence Variables", section 3.2
// https://arxiv.org/abs/1705.00907
//
// Features:
//       - commutative matching (sums and products, backtracking also in
//         powers and functions)
//       - more than two args with noncommutative functions
//       - one "global wildcard" (x^+) per sum or product
//
// TODO:
//       - "global wildcard"s (x^+)
//       - "zero wildcards" (matching superfluous wilds with 0(+) or 1(*)
//       - constant wildcards (those lacking specified symbols)
//       - commutative functions

#include "cmatcher.h"
#include "expairseq.h"
#include "symbol.h"
#include "wildcard.h"
#include "add.h"
#include "mul.h"
#include "power.h"
#include "function.h"
#include "fderivative.h"
#include "operators.h"
#include "utils.h"

#include <unistd.h>
#include <iostream>
#include <algorithm>


namespace GiNaC {

int CMatcher::level = 0;

static bool trivial_match(const ex& s, const ex& pattern, exmap& map)
{
	if (is_exactly_a<wildcard>(pattern)) {

		// Wildcard matches anything, but check whether we already have found
		// a match for that wildcard first (if so, the earlier match must be
		// the same expression)
                const auto& it = map.find(pattern);
                if (it != map.end())
		        return s.is_equal(ex_to<basic>(it->second));
		map[pattern] = s;
		return true;
	} 

        // Expression must be of the same type as the pattern
        if (s.bp->tinfo() != ex_to<basic>(pattern).tinfo())
                return false;

        // No subexpressions (such expressions are handled in cmatch()
        // So just compare the objects (there can't be
        // wildcards in the pattern)
        return s.bp->is_equal_same_type(ex_to<basic>(pattern));
}

inline bool is_ncfunc(const ex& e)
{
        return is_exactly_a<power>(e)
            or is_exactly_a<function>(e)
            or is_exactly_a<fderivative>(e)
            or is_exactly_a<exprseq>(e)
            or is_exactly_a<lst>(e);
}

inline bool is_func(const ex& e)
{
        return is_ncfunc(e)
            or is_a<expairseq>(e);
}

std::vector<exmap> ex::all_matches(const ex & pattern) const
{
        exmap map;
        std::vector<exmap> vec;
        if (not is_func(*this)) {
                bool ret = bp->match(pattern, map);
                if (ret)
                        vec.push_back(map);
                return vec;
        }
        CMatcher::level=0;
        CMatcher cm(*this, pattern, map);
        do {
                opt_exmap m = cm.get();
                if (m) {
                        vec.push_back(m.value());
                        cm.clear_ret();
                        if (cm.finished)
                                break;
                }
                else
                        break;
        }
        while (1);
        return vec;
}

inline bool CMatcher::get_alt(size_t i)
{
        CMatcher& cm = cms[i].value();
        if (cm.ret_val) {
                bool ret = cm.ret_val.value();
                if (ret) {
                        ret_map = cm.ret_map.value();
                        cm.ret_map.reset();
                }
                cm.ret_val.reset();
                return ret;
        }
        if (cm.finished)
                return false;
        bool ret;
        cm.map = map_repo[i];
        const opt_exmap& opm = cm.get();
        ret_map = opm;
        ret_val = ret = opm? true : false;
        cm.clear_ret();
        if (not cm.finished)
                finished = false;
        return ret;
}

CMatcher::~CMatcher()
{
}

opt_bool CMatcher::init()
{
        bool global_wild = false;
        const size_t snops = source.nops(), pnops = pattern.nops();
	if (is_exactly_a<wildcard>(pattern)) {
                const auto& it = map.find(pattern);
                finished = true;
                if (it != map.end()) {
		        if (not source.is_equal(ex_to<basic>(it->second)))
                                return false;
                        ret_map = map;
                        return true;
                }
		map[pattern] = source;
                ret_map = map;
		return true;
	} 
	if (ex_to<basic>(source).tinfo() != ex_to<basic>(pattern).tinfo()) {
                return false;
        }
        if (is_exactly_a<function>(source)
           and ex_to<function>(source).get_serial()
            != ex_to<function>(pattern).get_serial()) {
                return false;
        }
        if (snops < pnops) {
                return false;
        }
        symbolset oset = source.symbols();
        symbolset pset = substitute(pattern.wilds(), map);
        if (not subset_of(pset, oset)) {
                return false;
        }

        // Chop into terms
        for (size_t i=0; i<snops; i++) {
                ops.push_back(source.op(i));
        }
        for (size_t i=0; i<pnops; i++)
                pat.push_back(pattern.op(i).subs(map,
                                        subs_options::no_pattern));
        for (size_t i=0; i<pat.size(); ++i)
                if (is_exactly_a<wildcard>(pat[i]))
                        wild_ind.push_back(i);

        if (not is_ncfunc(source)) { // source is commutative
                if (snops > pnops) {
                        if (wild_ind.empty())
                                return false;
                        // we have a global wildcard, i.e., one that matches
                        // more than one term (x^+ in the paper) in a sum or
                        // product
                        global_wild = true;
                }
                // Check that all "constants" in the pattern are matched
                // Presets are already substituted now
                for (auto it1 = pat.begin(); it1 != pat.end(); ) {
                        if (haswild(*it1)) {
                                ++it1;
                                continue;
                        }
                        const auto& it2 = std::find_if(ops.begin(), ops.end(),
                                              [it1](const ex& e) {
                                                      return e.is_equal(*it1);
                                              } );
                        if (it2 != ops.end()) {
                                (void)ops.erase(it2);
                                it1 = pat.erase(it1);
                        }
                        else {
                                return false;
                        }
                }
        }
        else
        {
                // Check that all "constants" in the pattern are matched
                // Presets are already substituted now
                for (auto it1 = pat.begin(), it2 = ops.begin();
                                it1 != pat.end(); ) {
                        if (haswild(*it1)) {
                                ++it1;
                                ++it2;
                                continue;
                        }
                        if (it2->is_equal(*it1)) {
                                it2 = ops.erase(it2);
                                it1 = pat.erase(it1);
                        }
                        else {
                                return false;
                        }
                }
                if (pat.empty()) {
                        finished = true;
                        ret_map = map;
                        return true;
                }
        }
        if (wild_ind.empty() and ops.empty() and pat.empty()) {
                finished = true;
                ret_map = map;
                return true;
        }
        if (wild_ind.empty() and (ops.empty() or pat.empty()))
                throw std::runtime_error("matching gotcha");

        N = ops.size();
        P = pat.size();
        m_cmatch.assign(N, false);
        std::transform(pat.begin(), pat.end(), m_cmatch.begin(),
                        [](const ex& e) {
                        return is_func(e); } );
        if (P == 1) {
            if (not global_wild) {
                if (not m_cmatch[0]) {
                        ret_val = trivial_match(ops[0], pat[0], map);
                        if (ret_val.value())
                                ret_map = map;
                        finished = true;
                        return ret_val.value();
                }
            }
            else {
                // only the global wildcard left
                ret_val = true;
                ret_map = map;
                ex r;
                if (is_exactly_a<add>(source)) {
                        r = _ex0;
                        for (size_t i=0; i<ops.size(); ++i)
                                r += ops[i];
                }
                if (is_exactly_a<mul>(source)) {
                        r = _ex1;
                        for (size_t i=0; i<ops.size(); ++i)
                                r *= ops[i];
                }
                ret_map.value()[pat[0]] = r;
                finished = true;
                return true;
            }
        }
        
        // completely initialize cmatcher
        for (size_t i=0; i<N; ++i)
                cms.emplace_back();
        map_repo = std::vector<exmap>(N);

        finished = false;
        perm.reserve(P);        
        if (is_ncfunc(source)) {
                type = Type::noncomm;
                // perm is not used but needs to be filled
                for (size_t i=0; i<P; ++i)
                        perm.push_back(i);
                if (std::all_of(m_cmatch.begin(), m_cmatch.end(),
                                [](bool b) { return not b; } )) {
                        finished = true;
                        for (size_t i=0; i<P; ++i)
                                if (not trivial_match(ops[i], pat[i], map))
                                        return false;
                        ret_map = map;
                        return true;
                }
                return nullopt;
        }

        if (global_wild) {
                type = Type::comm_plus;
                --P;
                gws.reserve(P);
                next_combination(comb, P, N);
        }
        if (not global_wild) {
                type = Type::comm;
                for (size_t i=0; i<P; ++i)
                        perm.push_back(i);
        }
        return nullopt;
}

void CMatcher::run()
{
        clear_ret();
        if (finished or
            (type == Type::comm and perm.empty())) {
                ret_val = false;
                return;
        }
        switch (type) {
                case Type::noncomm:
                       noncomm_run();
                       return;
                case Type::comm:
                       no_global_wild();
                       return;
                case Type::comm_plus:
                       with_global_wild();
                       return;
                default:
                       throw std::runtime_error("can't happen");
        }
}

void CMatcher::noncomm_run()
{
        int ii = P;
        while (--ii >= 0) {
                if (m_cmatch[ii]
                    and cms[ii]
                    and not cms[ii].value().finished)
                        break;
        }
        size_t index;
        if (ii >= 0) {
                index = static_cast<size_t>(ii);
                while (unsigned(++ii) < P)
                        cms[ii].reset();
        }
        else // no cmatcher or all uninitialized
                index = 0;
        // The second loop increases index to get a new ops term
        do {
                const ex& e = ops[index];
                if (index == 0)
                        map_repo[0] = map;
                else
                        map_repo[index] = map_repo[index-1];
                // At this point we try matching p to e 
                const ex& p = pat[index];
                if (not m_cmatch[index]) {
                        // normal matching attempt
                        exmap m = map_repo[index];
                        bool ret = trivial_match(e, p, m);
                        if (ret) {
                                map_repo[index] = m;
                                continue;
                        }
                }
                else {
                        if (not cms[index]) {
                                exmap m = map_repo[index];
                                cms[index].emplace(e, p, m);
                        }
                        else {
                                cms[index].value().ret_val.reset();
                                cms[index].value().finished = false;
                        }
                        if (get_alt(index)) {
                                map_repo[index] = ret_map.value();
                                continue;
                        }
                        cms[index].reset();
                }
                bool alt_solution_found = false;
                int i = static_cast<int>(index);
                while (--i >= 0) {
                        if (cms[i]) {
                                if (i == 0)
                                        map_repo[0] = map;
                                else
                                        map_repo[i] = map_repo[i-1];
                                auto& cm = cms[i].value();
                                cm.ret_val.reset();
                                cm.map = map_repo[i];
                                if (get_alt(i)) {
                                        map_repo[i] = ret_map.value();
                                        index = i;
                                        alt_solution_found = true;
                                        break;
                                }
                                else
                                        cms[i].reset();
                        }
                }
                if (not alt_solution_found)
                        break;
        }
        while (++index < P);

        if (index >= P) {
                ret_val = true;
                ret_map = map_repo[P-1];
                finished = std::all_of(cms.begin(), cms.end(),
                      [](const opt_CMatcher& cm) { return not cm or cm.value().finished; } );
                return;
        }
        finished = true;
        ret_map.reset();
        ret_val = false;
}

void CMatcher::no_global_wild()
{

        perm_run(ops, pat);
}

void CMatcher::perm_run(const exvector& sterms, const exvector& pterms)
{
        // The outer loop goes though permutations of perm
        while (true) {
                int ii = P;
                while (--ii >= 0) {
                        if (cms[ii]
                            and not cms[ii].value().finished)
                                break;
                }
                size_t index;
                if (ii >= 0)
                        index = static_cast<size_t>(ii);
                else // no cmatcher or all uninitialized or finished
                        index = 0;
                while (unsigned(++ii) < P)
                        if (cms[ii])
                                cms[ii].reset();
                finished = true;
                // The second loop increases index to get a new ops term
                do {
                        const ex& e = sterms[index];
                        if (index == 0)
                                map_repo[0] = map;
                        else
                                map_repo[index] = map_repo[index-1];
                        // At this point we try matching p to e 
                        const ex& p = pterms[perm[index]];
                        if (not m_cmatch[perm[index]]) {
                                // normal matching attempt
                                exmap m = map_repo[index];
                                bool ret = trivial_match(e, p, m);
                                if (ret) {
                                        map_repo[index] = m;
                                        continue;
                                }
                        }
                        else {
                                if (not cms[index]) {
                                        exmap m = map_repo[index];
                                        cms[index].emplace(e, p, m);
                                }
                                else {
                                        cms[index].value().ret_val.reset();
                                        cms[index].value().finished = false;
                                }
                                if (get_alt(index)) {
                                        map_repo[index] = ret_map.value();
                                        continue;
                                }
                                cms[index].reset();
                        }
                        // unfinished cmatchers in this permutation?
                        bool alt_solution_found = false;
                        int i = static_cast<int>(index);
                        while (--i >= 0) {
                                if (cms[i]) {
                                        if (i == 0)
                                                map_repo[0] = map;
                                        else
                                                map_repo[i] = map_repo[i-1];
                                        auto& cm = cms[i].value();
                                        cm.ret_val.reset();
                                        cm.map = map_repo[i];
                                        if (get_alt(i)) {
                                                map_repo[i] = ret_map.value();
                                                index = i;
                                                alt_solution_found = true;
                                                break;
                                        }
                                        cms[i].reset();
                                }
                        }
                        if (not alt_solution_found)
                                break;
                }
                while (++index < P);

                if (index >= P) {
                        // give back one solution of this cmatch call
                        // possibly still permutations left
                        if (std::all_of(cms.begin(), cms.end(),
                              [](const opt_CMatcher& cm) { return not cm or cm.value().finished; } )) {
                                finished = not std::next_permutation(perm.begin(),
                                                perm.end());
                        }
                        else
                                finished = false;

                        ret_val = true;
                        ret_map = map_repo[P-1];
                        return;
                }
                else {
                        // no cmatch calls have alternative solutions
                        // to their cmatch, so get the next permutation
                        // that changes current index position
                        size_t old = perm[index];
                        while (perm[index] == old) {
                                bool more = std::next_permutation(perm.begin(),
                                                perm.end());
                                if (not more) {
                                        ret_val = false;
                                        ret_map.reset();
                                        finished = true;
                                        return;
                                }
                        }
                }
        }
        ret_val = false;
}

// The case with one global wildcard term in a sum or product.
//
// For all wildcard terms in the pattern, we remove it and define it
// as the global (so P is now one less); for all combinations of P-1
// source terms (out of N) the permutation vector is filled with 0,...,P-1
// Like before, all permutations of those P-1 source terms are matched
// against the reduced patterns. When matched, those source terms not
// chosen by the combination combined (+/*) are the solution of the
// global wildcard.
//
// We call comb_run() with the following arguments: a source term vector
// which contains those P-1 terms the current combination has picked out;
// and the P-1 pattern terms vector (missing the current global wildcard).
// comb_run() is nearly the same as perm_run(); we just need to adapt the
// indices pointing to the respective entries of m_cmatch, cms, and map_repo.
void CMatcher::with_global_wild()
{
        // NOTE: P was decremented earlier
        do {
                // loop through all global wildcards
                // caveat: is it used elsewhere?
                size_t wwi = wild_ind[wi];
                const wildcard& gw = ex_to<wildcard>(pat[wwi]);
                size_t ii;
                for (ii=0; ii<P+1; ++ii)
                        if (ii != wwi and haswild(pat[ii], gw))
                                break;
                if (ii < P+1)
                        continue;

                do {
                        // loop through all combinations
                        bool comb_finished = finished = false;
                        if (perm.empty()) { // new combination
                                for (size_t i=0; i<P; ++i)
                                        perm.push_back(i);
                                gwp = pat;
                                gwp.erase(gwp.begin() + wild_ind[wi]);
                                mcm = m_cmatch;
                                mcm.erase(mcm.begin() + wild_ind[wi]);
                                for (size_t i=0; i<P; ++i)
                                        gws.push_back(ops[comb[i]]);
                        }

                        comb_run(gws, gwp, mcm);

                        if (finished) {
                                finished = false;
                                comb_finished = true;
                                perm.clear();
                                gwp.clear();
                                gws.clear();
                        }
                        if (ret_val and ret_val.value())
                        {
                                ex gwe;
                                if (is_exactly_a<add>(source))
                                        gwe = _ex0;
                                if (is_exactly_a<mul>(source))
                                        gwe = _ex1;
                                std::vector<bool> t(P);
                                t.assign(ops.size(), false);
                                for (size_t i=0; i<P; ++i)
                                       t[comb[i]] = true;
                                if (is_exactly_a<add>(source))
                                        for (size_t i=0; i<ops.size(); ++i)
                                                if (not t[i])
                                                        gwe += ops[i];
                                if (is_exactly_a<mul>(source))
                                        for (size_t i=0; i<ops.size(); ++i)
                                                if (not t[i])
                                                        gwe *= ops[i];
                                ret_map.value()[pat[wild_ind[wi]]] = gwe;
                                finished = false;
                                if (comb_finished
                                    and not next_combination(comb, P, N)) {
                                        comb.clear();
                                        next_combination(comb, P, N);
                                        if (not (++wi < wild_ind.size()))
                                                finished = true;
                                }
                                return;
                        }
                        if (not comb_finished)
                                continue;
                }
                while (next_combination(comb, P, N));
                comb.clear();
                next_combination(comb, P, N);
        }
        while (++wi < wild_ind.size());
        finished = true;
        ret_val = false;
}

// sterms[index] == ops[comb[index]]
void CMatcher::comb_run(const exvector& sterms, const exvector& pterms,
                const std::vector<bool>& cmneeded)
{
        // The outer loop goes though permutations of perm
        while (true) {
                int ii = P;
                while (--ii >= 0) {
                        if (cms[comb[ii]]
                            and not cms[comb[ii]].value().finished)
                                break;
                }
                size_t index;
                if (ii >= 0) {
                        index = static_cast<size_t>(ii);
                }
                else // no cmatcher or all uninitialized
                        index = 0;
                while (unsigned(++ii) < P)
                        cms[comb[ii]].reset();
                finished = true;
                // The second loop increases index to get a new ops term
                do {
                        const ex& e = sterms[index];
                        if (index == 0)
                                map_repo[comb[0]] = map;
                        else
                                map_repo[comb[index]] = map_repo[comb[index-1]];
                        // At this point we try matching p to e 
                        const ex& p = pterms[perm[index]];
                        if (not cmneeded[perm[index]]) {
                                // normal matching attempt
                                exmap m = map_repo[comb[index]];
                                bool ret = trivial_match(e, p, m);
                                if (ret) {
                                        map_repo[comb[index]] = m;
                                        continue;
                                }
                        }
                        else {
                                if (not cms[comb[index]]) {
                                        exmap m = map_repo[comb[index]];
                                        cms[comb[index]].emplace(e, p, m);
                                }
                                else {
                                        cms[comb[index]].value().ret_val.reset();
                                        cms[comb[index]].value().finished = false;
                                }
                                if (get_alt(comb[index])) {
                                        map_repo[comb[index]] = ret_map.value();
                                        continue;
                                }
                                else
                                        cms[comb[index]].reset();
                        }
                        // unfinished cmatchers in this permutation?
                        bool alt_solution_found = false;
                        int i = static_cast<int>(index);
                        while (--i >= 0) {
                                if (cms[comb[i]]) {
                                        if (i == 0)
                                                map_repo[comb[0]] = map;
                                        else
                                                map_repo[comb[i]] = map_repo[comb[i-1]];
                                        auto& cm = cms[comb[i]].value();
                                        cm.ret_val.reset();
                                        cm.map = map_repo[comb[i]];
                                        if (get_alt(comb[i])) {
                                                map_repo[comb[i]] = ret_map.value();
                                                index = i;
                                                alt_solution_found = true;
                                                break;
                                        }
                                        cms[comb[i]].reset();
                                }
                        }
                        if (not alt_solution_found)
                                break;
                }
                while (++index < P);

                if (index >= P) {
                        // give back one solution of this cmatch call
                        // possibly still permutations left
                        size_t i;
                        for (i=0; i<P; ++i)
                                if (cms[comb[i]]
                                    and not cms[comb[i]].value().finished)
                                        break;
                        if (i >= P)
                                finished =  not std::next_permutation(perm.begin(),
                                                perm.end());
                        else
                                finished = false;
                        ret_val = true;
                        ret_map = map_repo[comb[P-1]];
                        return;
                }
                else {
                        // no cmatch calls have alternative solutions
                        // to their cmatch, so get the next permutation
                        // that changes current index position
                        size_t old = perm[index];
                        while (perm[index] == old) {
                                bool more = std::next_permutation(perm.begin(),
                                                perm.end());
                                if (not more) {
                                        ret_val = false;
                                        ret_map.reset();
                                        finished = true;
                                        return;
                                }
                        }
                }
        }
        ret_val = false;
}

}
