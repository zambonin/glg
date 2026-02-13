#include "nonsing.h"
#include "util/fq_nmod_mat_extra.h"
#include "util/rand.h"

void inner_nonsing(fq_nmod_mat_t A, fq_nmod_mat_t T, const fq_nmod_ctx_t ctx,
                   bit_buffer_t *buf) {
  const slong n = fq_nmod_mat_nrows(A, ctx);
  if (n == 1) {
    fq_nmod_mat_one(A, ctx);
    fq_nmod_mat_entry_rand_not_zero_buf(T, 0, 0, ctx, buf);
    return;
  }

  fq_nmod_mat_t A_sub;
  fq_nmod_mat_init(A_sub, n - 1, n - 1, ctx);

  fq_nmod_mat_t T_sub;
  fq_nmod_mat_init(T_sub, n - 1, n - 1, ctx);

  inner_nonsing(A_sub, T_sub, ctx, buf);

  fq_nmod_mat_t v;
  fq_nmod_mat_init(v, n, 1, ctx);
  fq_nmod_mat_randtest_not_zero_buf(v, ctx, buf);

  const slong r = fq_nmod_mat_first_pos_entry(v, ctx);

  fq_nmod_t one;
  fq_nmod_init(one, ctx);
  fq_nmod_one(one, ctx);

  fq_nmod_mat_zero(A, ctx);
  fq_nmod_mat_entry_set(A, 0, r, one, ctx);

  fq_nmod_mat_set_minor(A, 0, r, A_sub, ctx);

  fq_nmod_mat_t window;
  fq_nmod_mat_window_init(window, A, 1, r, n, r + 1, ctx);
  fq_nmod_mat_randtest_buf(window, ctx, buf);
  fq_nmod_mat_window_clear(window, ctx);

  fq_nmod_mat_zero(T, ctx);
  fq_nmod_mat_set_minor(T, r, r, T_sub, ctx);

  fq_nmod_mat_t vT;
  fq_nmod_mat_init(vT, 1, n, ctx);
  fq_nmod_mat_transpose(vT, v, ctx);

  fq_nmod_mat_window_init(window, T, r, 0, r + 1, n, ctx);
  fq_nmod_mat_set(window, vT, ctx);
  fq_nmod_mat_window_clear(window, ctx);

  fq_nmod_mat_entry_set(T, r, r, one, ctx);

  fq_nmod_mat_clear(vT, ctx);
  fq_nmod_mat_clear(A_sub, ctx);
  fq_nmod_mat_clear(T_sub, ctx);
  fq_nmod_mat_clear(v, ctx);
  fq_nmod_clear(one, ctx);
}

void nonsing(fq_nmod_mat_t M, const fq_nmod_ctx_t ctx, flint_rand_t state) {
  bit_buffer_t buf;
  bit_buffer_init(&buf, state);

  const slong n = fq_nmod_mat_nrows(M, ctx);
  fq_nmod_mat_zero(M, ctx);

  fq_nmod_mat_t A;
  fq_nmod_mat_init(A, n, n, ctx);

  fq_nmod_mat_t T;
  fq_nmod_mat_init(T, n, n, ctx);

  inner_nonsing(A, T, ctx, &buf);

  fq_nmod_mat_mul(M, A, T, ctx);

  fq_nmod_mat_clear(T, ctx);
  fq_nmod_mat_clear(A, ctx);
  bit_buffer_clear(&buf);
}
