#include "rej.h"

void rej(fq_nmod_mat_t M, const fq_nmod_ctx_t ctx, flint_rand_t state) {
  fq_nmod_mat_t inv;
  fq_nmod_mat_init(inv, M->r, M->c, ctx);

  do {
    fq_nmod_mat_randtest(M, state, ctx);
  } while (!fq_nmod_mat_inv(inv, M, ctx));

  fq_nmod_mat_clear(inv, ctx);
}
