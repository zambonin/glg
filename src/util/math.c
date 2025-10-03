#include "util/math.h"

void gl_order(fmpz_t rop, const fq_nmod_ctx_t ctx, const slong n) {
  fmpz_set_ui(rop, 1);

  fmpz_t q_;
  fmpz_init(q_);
  fq_nmod_ctx_order(q_, ctx);

  fmpz_t q_pow_n;
  fmpz_init(q_pow_n);
  fmpz_pow_ui(q_pow_n, q_, n);

  fmpz_t q_pow_i;
  fmpz_init_set_ui(q_pow_i, 1);

  fmpz_t count;
  fmpz_init(count);

  for (slong i = 0; i < n; ++i) {
    fmpz_sub(count, q_pow_n, q_pow_i);
    fmpz_mul(rop, rop, count);
    fmpz_mul(q_pow_i, q_pow_i, q_);
  }

  fmpz_clear(count);
  fmpz_clear(q_pow_i);
  fmpz_clear(q_pow_n);
  fmpz_clear(q_);
}

void complete_basis(fq_nmod_mat_t B, const fq_nmod_mat_t S, const slong n,
                    const fq_nmod_ctx_t ctx) {
  const slong i = fq_nmod_mat_ncols(S, ctx);
  fq_nmod_mat_zero(B, ctx);

  if (i >= n) {
    fq_nmod_mat_set(B, S, ctx);
    return;
  }

  fq_nmod_mat_t M;
  fq_nmod_mat_init(M, n, i + n, ctx);
  fq_nmod_mat_zero(M, ctx);
  for (slong r = 0; r < n; ++r) {
    fq_nmod_one(fq_nmod_mat_entry(M, r, i + r), ctx);
  }

  fq_nmod_mat_t view;
  fq_nmod_mat_window_init(view, M, 0, 0, n, i, ctx);
  fq_nmod_mat_set(view, S, ctx);
  fq_nmod_mat_window_clear(view, ctx);

  const slong rank = fq_nmod_mat_rref(M, M, ctx);

  slong l = 0;
  slong s = -1;

  for (slong r = 0; r < rank && l < n; ++r) {
    for (slong c = s + 1; c < i + n; ++c) {
      if (fq_nmod_is_zero(fq_nmod_mat_entry(M, r, c), ctx)) {
        continue;
      }

      s = c;
      if (c < i) {
        for (slong idx = 0; idx < n; ++idx) {
          fq_nmod_set(fq_nmod_mat_entry(B, idx, l),
                      fq_nmod_mat_entry(S, idx, c), ctx);
        }
      } else {
        for (slong idx = 0; idx < n; ++idx) {
          fq_nmod_set_ui(fq_nmod_mat_entry(B, idx, l), idx == (c - i), ctx);
        }
      }
      l++;
      break;
    }
  }

  fq_nmod_mat_clear(M, ctx);
}
