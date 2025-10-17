#ifndef UTIL_H
#define UTIL_H

#include "common.h"

void fq_nmod_mat_randtest_not_zero(fq_nmod_mat_t v, const fq_nmod_ctx_t ctx,
                                   flint_rand_t state);

slong fq_nmod_mat_first_pos_entry(fq_nmod_mat_t v, const fq_nmod_ctx_t ctx);

void fq_nmod_mat_set_minor(fq_nmod_mat_t mat, const ulong skip_i,
                           const ulong skip_j, const fq_nmod_mat_t src,
                           const fq_nmod_ctx_t ctx);

void fq_nmod_mat_scalar_mul(fq_nmod_mat_t B, const fq_nmod_mat_t A,
                            const fq_nmod_t c, const fq_nmod_ctx_t ctx);

void fmpz_randlimb_m(fmpz_t f, flint_rand_t state, const fmpz_t m);

#endif // UTIL_H
