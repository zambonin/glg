#include "unrank.h"
#include "comb/gl.h"
#include "util/fq_nmod_mat_extra.h"
#include "util/math.h"
#include "util/rand.h"

void unrank(fq_nmod_mat_t M, const fq_nmod_ctx_t ctx, flint_rand_t state) {
  bit_buffer_t buf;
  bit_buffer_init(&buf, state);

  fmpz_t count;
  fmpz_init(count);

  fmpz_t rank;
  fmpz_init_set_ui(rank, 0);

  gl_order(count, ctx, M->r);

  if (fmpz_cmp_ui(count, 0) > 0) {
    fast_dice_roller(rank, count, &buf);
  }

  unrank_gl(M, rank, ctx);

  fmpz_clear(rank);
  fmpz_clear(count);
  bit_buffer_clear(&buf);
}
