#ifndef VSPACE_H
#define VSPACE_H

#include "common.h"

void unrank_vspace(fq_nmod_mat_t v, const fmpz_t r, const fq_nmod_ctx_t ctx);

void rank_vspace(fmpz_t r, const fq_nmod_mat_t v, const fq_nmod_ctx_t ctx);

#endif // VSPACE_H
