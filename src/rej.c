#include "rej.h"
#include "util/fq_nmod_mat_extra.h"
#include "util/rand.h"

void rej(fq_nmod_mat_t M, const fq_nmod_ctx_t ctx, flint_rand_t state) {
  bit_buffer_t buf;
  bit_buffer_init(&buf, state);

  fq_nmod_mat_t inv;
  fq_nmod_mat_init(inv, M->r, M->c, ctx);

  do {
    fq_nmod_mat_randtest_buf(M, ctx, &buf);
  } while (!fq_nmod_mat_inv(inv, M, ctx));

  fq_nmod_mat_clear(inv, ctx);
  bit_buffer_clear(&buf);
}
