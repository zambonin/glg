#ifndef MATH_H
#define MATH_H

#include "common.h"

void gl_order(fmpz_t rop, const fq_nmod_ctx_t ctx, const slong n);

void complete_basis(fq_nmod_mat_t B, const fq_nmod_mat_t S, const slong n,
                    const fq_nmod_ctx_t ctx);

#endif // MATH_H
