#ifndef SSPACE_H
#define SSPACE_H

#include "common.h"

void unrank_subspace(fq_nmod_mat_t v, const fq_nmod_mat_t B, const fmpz_t r,
                     const fq_nmod_ctx_t ctx);

void rank_subspace(fmpz_t r, const fq_nmod_mat_t v, const fq_nmod_mat_t B,
                   const fq_nmod_ctx_t ctx);

#endif // SSPACE_H
