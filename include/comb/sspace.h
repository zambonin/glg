#ifndef SSPACE_H
#define SSPACE_H

#include "common.h"

void unrank_subspace(fq_nmod_mat_t v, const fq_nmod_mat_t B, const fmpz_t r,
                     const fq_nmod_ctx_t ctx);

void rank_subspace(fmpz_t r, const fq_nmod_mat_t v, const fq_nmod_mat_t B,
                   const fq_nmod_ctx_t ctx);

void unrank_subspace_pivots(fq_nmod_mat_t v, const fq_nmod_mat_t B,
                            const slong *J, const slong n_avail, const fmpz_t r,
                            const fq_nmod_ctx_t ctx, slong *m_last_out);

void rank_subspace_pivots(fmpz_t r, const fq_nmod_mat_t v,
                          const fq_nmod_mat_t B, const slong *J,
                          const slong n_avail, const fq_nmod_ctx_t ctx,
                          slong *m_last_out);

#endif // SSPACE_H
