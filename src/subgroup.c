#include "subgroup.h"
#include "util/fq_nmod_mat_extra.h"
#include "util/rand.h"

void inner_subgroup(fq_nmod_mat_t M, const fq_nmod_ctx_t ctx,
                    bit_buffer_t *buf) {
  const slong n = fq_nmod_mat_nrows(M, ctx);

  if (n == 1) {
    fq_nmod_mat_entry_rand_not_zero_buf(M, 0, 0, ctx, buf);
    return;
  }

  fq_nmod_mat_t Mp;
  fq_nmod_mat_init(Mp, n, n, ctx);
  fq_nmod_mat_zero(Mp, ctx);

  fq_nmod_mat_t A;
  fq_nmod_mat_init(A, n, n, ctx);
  fq_nmod_mat_one(A, ctx);

  fq_nmod_mat_t X;
  fq_nmod_mat_init(X, n, n, ctx);
  fq_nmod_mat_zero(X, ctx);

  fq_nmod_mat_t tmp;
  fq_nmod_mat_init(tmp, n, n, ctx);

  fq_nmod_mat_t M_sub;
  fq_nmod_mat_init(M_sub, n - 1, n - 1, ctx);

  fq_nmod_mat_t v;
  fq_nmod_mat_init(v, n, 1, ctx);

  fq_nmod_mat_t view;

  fq_nmod_t v_k_inv;
  fq_nmod_init(v_k_inv, ctx);

  inner_subgroup(M_sub, ctx, buf);

  fq_nmod_one(fq_nmod_mat_entry(Mp, 0, 0), ctx);

  fq_nmod_mat_window_init(view, Mp, 1, 1, n, n, ctx);
  fq_nmod_mat_set(view, M_sub, ctx);
  fq_nmod_mat_window_clear(view, ctx);

  fq_nmod_mat_entry_rand_not_zero_buf(A, 0, 0, ctx, buf);

  fq_nmod_mat_window_init(view, A, 0, 1, 1, n, ctx);
  fq_nmod_mat_randtest_buf(view, ctx, buf);
  fq_nmod_mat_window_clear(view, ctx);

  fq_nmod_mat_randtest_not_zero_buf(v, ctx, buf);

  const slong k = fq_nmod_mat_first_pos_entry(v, ctx);
  fq_nmod_inv(v_k_inv, fq_nmod_mat_entry(v, k, 0), ctx);

  fq_nmod_mat_scalar_mul(v, v, v_k_inv, ctx);

  fq_nmod_mat_window_init(view, X, 0, 0, n, 1, ctx);
  fq_nmod_mat_set(view, v, ctx);
  fq_nmod_mat_window_clear(view, ctx);

  slong col = 1;
  for (slong basis = 0; basis < n; ++basis) {
    if (basis == k) {
      continue;
    }

    fq_nmod_one(fq_nmod_mat_entry(X, basis, col), ctx);
    ++col;
  }

  fq_nmod_mat_mul(tmp, A, Mp, ctx);
  fq_nmod_mat_mul(M, X, tmp, ctx);

  fq_nmod_clear(v_k_inv, ctx);
  fq_nmod_mat_clear(v, ctx);
  fq_nmod_mat_clear(M_sub, ctx);
  fq_nmod_mat_clear(tmp, ctx);
  fq_nmod_mat_clear(X, ctx);
  fq_nmod_mat_clear(A, ctx);
  fq_nmod_mat_clear(Mp, ctx);
}

void subgroup(fq_nmod_mat_t M, const fq_nmod_ctx_t ctx, flint_rand_t state) {
  bit_buffer_t buf;
  bit_buffer_init(&buf, state);
  inner_subgroup(M, ctx, &buf);
  bit_buffer_clear(&buf);
}
