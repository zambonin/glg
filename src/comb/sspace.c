#include "comb/sspace.h"
#include "comb/vspace.h"

void _greedy_complement(slong *J, const slong n, const fq_nmod_mat_t B,
                        const fq_nmod_ctx_t ctx) {
  if (B == NULL) {
    for (slong i = 0; i < n; ++i) {
      J[i] = i;
    }
    return;
  }

  const slong n_i = fq_nmod_mat_ncols(B, ctx);

  fq_nmod_mat_t rref;
  fq_nmod_mat_init(rref, n_i, n, ctx);
  fq_nmod_mat_transpose(rref, B, ctx);
  fq_nmod_mat_rref(rref, rref, ctx);

  slong count = 0;
  slong col = 0;
  for (slong j = 0; j < n_i; ++j) {
    while (col < n && fq_nmod_is_zero(fq_nmod_mat_entry(rref, j, col), ctx)) {
      J[count++] = col++;
    }
    ++col;
  }

  while (count < n - n_i) {
    J[count++] = col++;
  }

  fq_nmod_mat_clear(rref, ctx);
}

void unrank_subspace(fq_nmod_mat_t v, const fq_nmod_mat_t B, const fmpz_t r,
                     const fq_nmod_ctx_t ctx) {
  const slong n = fq_nmod_mat_nrows(v, ctx);
  const slong n_i = (B == NULL) ? 0 : fq_nmod_mat_ncols(B, ctx);

  slong *J = flint_malloc((n - n_i) * sizeof(slong));

  _greedy_complement(J, n, B, ctx);
  unrank_subspace_pivots(v, B, J, n - n_i, r, ctx, NULL);

  flint_free(J);
}

void rank_subspace(fmpz_t r, const fq_nmod_mat_t v, const fq_nmod_mat_t B,
                   const fq_nmod_ctx_t ctx) {
  const slong n = fq_nmod_mat_nrows(v, ctx);
  const slong n_i = (B == NULL) ? 0 : fq_nmod_mat_ncols(B, ctx);

  slong *J = flint_malloc((n - n_i) * sizeof(slong));

  _greedy_complement(J, n, B, ctx);
  rank_subspace_pivots(r, v, B, J, n - n_i, ctx, NULL);

  flint_free(J);
}

void unrank_subspace_pivots(fq_nmod_mat_t v, const fq_nmod_mat_t B,
                            const slong *J, const slong n_avail, const fmpz_t r,
                            const fq_nmod_ctx_t ctx, slong *m_last_out) {
  const slong n = fq_nmod_mat_nrows(v, ctx);
  const slong n_i = (B == NULL) ? 0 : fq_nmod_mat_ncols(B, ctx);

  fq_nmod_mat_zero(v, ctx);

  if (m_last_out) {
    *m_last_out = -1;
  }

  if (n_avail == 0 && n_i == 0) {
    return;
  }

  fmpz_t q;
  fmpz_init(q);
  fq_nmod_ctx_order(q, ctx);

  fmpz_t q_pow_i;
  fmpz_init(q_pow_i);
  fmpz_pow_ui(q_pow_i, q, n_i);

  fmpz_t r_prime;
  fmpz_init(r_prime);
  fmpz_add(r_prime, r, q_pow_i);

  fq_nmod_t tmp;
  fq_nmod_init(tmp, ctx);

  fq_nmod_mat_t c;
  fq_nmod_mat_init(c, n, 1, ctx);
  unrank_vspace(c, r_prime, ctx);

  for (slong j = 0; j < n_i; ++j) {
    if (fq_nmod_is_zero(fq_nmod_mat_entry(c, j, 0), ctx)) {
      continue;
    }
    for (slong row = 0; row < n; ++row) {
      fq_nmod_mul(tmp, fq_nmod_mat_entry(B, row, j),
                  fq_nmod_mat_entry(c, j, 0), ctx);
      fq_nmod_add(fq_nmod_mat_entry(v, row, 0), fq_nmod_mat_entry(v, row, 0),
                  tmp, ctx);
    }
  }

  for (slong m = 0; m < n_avail; ++m) {
    if (m_last_out && *m_last_out == -1 &&
        !fq_nmod_is_zero(fq_nmod_mat_entry(c, n_i + m, 0), ctx)) {
      *m_last_out = m;
    }
    fq_nmod_add(fq_nmod_mat_entry(v, J[m], 0), fq_nmod_mat_entry(v, J[m], 0),
                fq_nmod_mat_entry(c, n_i + m, 0), ctx);
  }

  fq_nmod_mat_clear(c, ctx);
  fq_nmod_clear(tmp, ctx);
  fmpz_clear(r_prime);
  fmpz_clear(q_pow_i);
  fmpz_clear(q);
}

void rank_subspace_pivots(fmpz_t r, const fq_nmod_mat_t v,
                          const fq_nmod_mat_t B, const slong *J,
                          const slong n_avail, const fq_nmod_ctx_t ctx,
                          slong *m_last_out) {
  const slong n = fq_nmod_mat_nrows(v, ctx);
  const slong n_i = (B == NULL) ? 0 : fq_nmod_mat_ncols(B, ctx);

  fmpz_zero(r);

  if (m_last_out) {
    *m_last_out = -1;
  }

  if (n_avail == 0 && n_i == 0) {
    return;
  }

  fmpz_t q;
  fmpz_init(q);
  fq_nmod_ctx_order(q, ctx);

  fmpz_t q_pow_i;
  fmpz_init(q_pow_i);
  fmpz_pow_ui(q_pow_i, q, n_i);

  fq_nmod_mat_t F;
  fq_nmod_mat_init(F, n, n, ctx);

  for (slong j = 0; j < n_i; ++j) {
    for (slong row = 0; row < n; ++row) {
      fq_nmod_mat_entry_set(F, row, j, fq_nmod_mat_entry(B, row, j), ctx);
    }
  }

  for (slong m = 0; m < n_avail; ++m) {
    fq_nmod_one(fq_nmod_mat_entry(F, J[m], n_i + m), ctx);
  }

  fq_nmod_mat_t c;
  fq_nmod_mat_init(c, n, 1, ctx);
  (void)fq_nmod_mat_solve(c, F, v, ctx);

  if (m_last_out) {
    for (slong m = 0; m < n_avail; ++m) {
      if (!fq_nmod_is_zero(fq_nmod_mat_entry(c, n_i + m, 0), ctx)) {
        *m_last_out = m;
        break;
      }
    }
  }

  fmpz_t r_prime;
  fmpz_init(r_prime);
  rank_vspace(r_prime, c, ctx);

  fmpz_sub(r, r_prime, q_pow_i);

  fq_nmod_mat_clear(c, ctx);
  fq_nmod_mat_clear(F, ctx);
  fmpz_clear(r_prime);
  fmpz_clear(q_pow_i);
  fmpz_clear(q);
}
