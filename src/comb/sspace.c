#include "comb/sspace.h"
#include "comb/vspace.h"
#include "util/math.h"

void unrank_subspace(fq_nmod_mat_t v, const fq_nmod_mat_t B, const fmpz_t r,
                     const fq_nmod_ctx_t ctx) {
  const slong n = fq_nmod_mat_nrows(B, ctx);
  const slong n_i = fq_nmod_mat_ncols(B, ctx);
  const slong n_o = n - n_i;

  fq_nmod_mat_zero(v, ctx);
  if (n_o == 0) {
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

  fq_nmod_mat_t c;
  fq_nmod_mat_init(c, n, 1, ctx);
  unrank_vspace(c, r_prime, ctx);

  fq_nmod_mat_t F;
  fq_nmod_mat_init(F, n, n, ctx);
  complete_basis(F, B, n, ctx);

  fq_nmod_mat_mul(v, F, c, ctx);

  fq_nmod_mat_clear(F, ctx);
  fq_nmod_mat_clear(c, ctx);
  fmpz_clear(r_prime);
  fmpz_clear(q_pow_i);
  fmpz_clear(q);
}

void rank_subspace(fmpz_t r, const fq_nmod_mat_t v, const fq_nmod_mat_t B,
                   const fq_nmod_ctx_t ctx) {
  const slong n = fq_nmod_mat_nrows(B, ctx);
  const slong n_i = fq_nmod_mat_ncols(B, ctx);
  const slong n_o = n - n_i;

  fmpz_zero(r);
  if (n_o == 0) {
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
  complete_basis(F, B, n, ctx);

  fq_nmod_mat_t c;
  fq_nmod_mat_init(c, n, 1, ctx);
  (void)fq_nmod_mat_solve(c, F, v, ctx);

  fmpz_t r_prime;
  fmpz_init(r_prime);
  rank_vspace(r_prime, c, ctx);

  fmpz_sub(r, r_prime, q_pow_i);

  fq_nmod_mat_clear(F, ctx);
  fq_nmod_mat_clear(c, ctx);
  fmpz_clear(r_prime);
  fmpz_clear(q_pow_i);
  fmpz_clear(q);
}
