#include "lu.h"

void lu(fq_nmod_mat_t M, const fq_nmod_ctx_t ctx, flint_rand_t state) {
  fq_nmod_mat_t L;
  fq_nmod_mat_init(L, M->r, M->c, ctx);

  fq_nmod_mat_t U;
  fq_nmod_mat_init(U, M->r, M->c, ctx);

  fq_nmod_mat_randtril(L, state, 1, ctx);
  fq_nmod_mat_randtriu(U, state, 0, ctx);

  fq_nmod_mat_mul(M, L, U, ctx);

  fq_nmod_mat_clear(U, ctx);
  fq_nmod_mat_clear(L, ctx);
}
