#ifndef GL_H
#define GL_H

#include "common.h"

void unrank_gl(fq_nmod_mat_t M, const fmpz_t r, const fq_nmod_ctx_t ctx);

void rank_gl(fmpz_t r, const fq_nmod_mat_t M, const fq_nmod_ctx_t ctx);

#endif // GL_H
