#include "gl.h"
#include "util/math.h"
#include "unrank.h"

void unrank(fq_nmod_mat_t M, const fq_nmod_ctx_t ctx, flint_rand_t state) {
  fmpz_t count;
  fmpz_init(count);

  fmpz_t rank;
  fmpz_init_set_ui(rank, 0);

  gl_order(count, ctx, M->r);
  if (fmpz_cmp_ui(count, 0) > 0) {
    fmpz_randm(rank, state, count);
  }

  unrank_gl(M, rank, ctx);

  fmpz_clear(rank);
  fmpz_clear(count);
}
