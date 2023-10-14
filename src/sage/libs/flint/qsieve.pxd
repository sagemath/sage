# distutils: libraries = flint
# distutils: depends = flint/qsieve.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    mp_limb_t qsieve_knuth_schroeppel(qs_t qs_inf)
    # Return the Knuth-Schroeppel multiplier for the `n`, integer to be factored
    # based upon the Knuth-Schroeppel function.

    mp_limb_t qsieve_primes_init(qs_t qs_inf)
    # Compute the factor base prime along with there inverse for `kn`, where `k`
    # is Knuth-Schroeppel multiplier and `n` is the integer to be factored. It
    # also computes the square root of `kn` modulo factor base primes.

    mp_limb_t qsieve_primes_increment(qs_t qs_inf, mp_limb_t delta)
    # It increase the number of factor base primes by amount 'delta' and
    # calculate inverse of those primes along with the square root of `kn` modulo
    # those primes.

    void qsieve_init_A0(qs_t qs_inf)
    # First it chooses the possible range of factor of `A _0`, based on the number
    # of bits in optimal value of `A _0`. It tries to select range such that we have
    # plenty of primes to choose from as well as number of factor in `A _0` are
    # sufficient. For input of size less than 130 bit, this selection method doesn't
    # work therefore we randomly generate 2 or 3-subset of all the factor base prime
    # as the factor of `A _0`.
    # Otherwise, if we have to select `s` factor for `A _0`, we generate `s - 1`-
    # subset from odd indices of the possible range of factor and then search last
    # factor using binary search from the even indices of possible range of factor
    # such that value of `A _0` is close to it's optimal value.

    void qsieve_next_A0(qs_t qs_inf)
    # Find next candidate for `A _0` as follows:
    # generate next lexicographic `s - 1`-subset from the odd indices of possible
    # range of factor base and choose the last factor from even indices using binary
    # search so that value `A _0` is close to it's optimal value.

    void qsieve_compute_pre_data(qs_t qs_inf)
    # Precompute all the data associated with factor's of `A _0`, since `A _0` is going
    # to be fixed for several `A`.

    void qsieve_init_poly_first(qs_t qs_inf)
    # Initializes the value of `A = q _0 * A _0`, where `q _0` is non-factor base prime.
    # precompute the data necessary for generating different `B` value using grey code
    # formula. Combine the data calculated for the factor of `A _0` along with the
    # parameter `q _0` to obtain data as for factor of `A`. It also calculates the sieve
    # offset for all the factor base prime, for first polynomial.

    void qsieve_init_poly_next(qs_t qs_inf, slong i)
    # Generate next polynomial or next `B` value for particular `A` and also updates the
    # sieve offsets for all the factor base prime, for this `B` value.

    void qsieve_compute_C(fmpz_t C, qs_t qs_inf, qs_poly_t poly)
    # Given `A` and `B`, calculate `C = (B ^2 - A) / N`.

    void qsieve_do_sieving(qs_t qs_inf, unsigned char * sieve, qs_poly_t poly)
    # First initialize the sieve array to zero, then for each `p \in` ``factor base``, add
    # `\log_2(p)` to the locations `\operatorname{soln1} _p + i * p` and `\operatorname{soln2} _p + i * p` for
    # `i = 0, 1, 2,\dots`, where `\operatorname{soln1} _p` and `\operatorname{soln2} _p` are the sieve offsets calculated
    # for `p`.

    void qsieve_do_sieving2(qs_t qs_inf, unsigned char * seive, qs_poly_t poly)
    # Perform the same task as above but instead of sieving over whole array at once divide
    # the array in blocks and then sieve over each block for all the primes in factor base.

    slong qsieve_evaluate_candidate(qs_t qs_inf, ulong i, unsigned char * sieve, qs_poly_t poly)
    # For location `i` in sieve array value at which, is greater than sieve threshold, check
    # the value of `Q(x)` at position `i` for smoothness. If value is found to be smooth then
    # store it for later processing, else check the residue for the partial if it is found to
    # be partial then store it for late processing.

    slong qsieve_evaluate_sieve(qs_t qs_inf, unsigned char * sieve, qs_poly_t poly)
    # Scan the sieve array for location at, which accumulated value is greater than sieve
    # threshold.

    slong qsieve_collect_relations(qs_t qs_inf, unsigned char * sieve)
    # Call for initialization of polynomial, sieving, and scanning of sieve
    # for all the possible polynomials for particular hypercube i.e. `A`.

    void qsieve_write_to_file(qs_t qs_inf, mp_limb_t prime, fmpz_t Y, qs_poly_t poly)
    # Write a relation to the file. Format is as follows,
    # first write large prime, in case of full relation it is 1, then write exponent
    # of small primes, then write number of factor followed by offset of factor in
    # factor base and their exponent and at last value of `Q(x)` for particular relation.
    # each relation is written in new line.

    hash_t * qsieve_get_table_entry(qs_t qs_inf, mp_limb_t prime)
    # Return the pointer to the location of 'prime' is hash table if it exist, else
    # create and entry for it in hash table and return pointer to that.

    void qsieve_add_to_hashtable(qs_t qs_inf, mp_limb_t prime)
    # Add 'prime' to the hast table.

    relation_t qsieve_parse_relation(qs_t qs_inf, char * str)
    # Given a string representation of relation from the file, parse it to obtain
    # all the parameters of relation.

    relation_t qsieve_merge_relation(qs_t qs_inf, relation_t  a, relation_t  b)
    # Given two partial relation having same large prime, merge them to obtain a full
    # relation.

    int qsieve_compare_relation(const void * a, const void * b)
    # Compare two relation based on, first large prime, then number of factor and then
    # offsets of factor in factor base.

    int qsieve_remove_duplicates(relation_t * rel_list, slong num_relations)
    # Remove duplicate from given list of relations by sorting relations in the list.

    void qsieve_insert_relation2(qs_t qs_inf, relation_t * rel_list, slong num_relations)
    # Given a list of relations, insert each relation from the list into the matrix for
    # further processing.

    int qsieve_process_relation(qs_t qs_inf)
    # After we have accumulated required number of relations, first process the file by
    # reading all the relations, removes singleton. Then merge all the possible partial
    # to obtain full relations.

    void qsieve_factor(fmpz_factor_t factors, const fmpz_t n)
    # Factor `n` using the quadratic sieve method. It is required that `n` is not a
    # prime and not a perfect power. There is no guarantee that the factors found will
    # be prime, or distinct.
