#ifndef POLY_H
#define POLY_H

#include "common.h"

void rank_poly(fmpz_t r, const fq_nmod_t poly, const fq_nmod_ctx_t ctx);

void unrank_poly(fq_nmod_t poly, const fmpz_t r, const fq_nmod_ctx_t ctx);

#endif // POLY_H
