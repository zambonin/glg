#include <string.h>

#include <flint/perm.h>

#include "schubert.h"

void schubert(fq_nmod_mat_t M, const fq_nmod_ctx_t ctx, flint_rand_t state) {
  const slong n = fq_nmod_mat_nrows(M, ctx);
  if (n == 1) {
    fq_nmod_rand_not_zero(fq_nmod_mat_entry(M, 0, 0), state, ctx);
    return;
  }

  slong *sigma = _perm_init(n);
  sigma[0] = 0;

  slong *sigma_inv = _perm_init(n);

  fmpz_t q;
  fmpz_init(q);
  fq_nmod_ctx_order(q, ctx);

  fmpz_t q_pow_i;
  fmpz_init_set(q_pow_i, q);

  fmpz_t total;
  fmpz_init_set_ui(total, 1);

  fmpz_t acc_weight;
  fmpz_init(acc_weight);

  fmpz_t weight;
  fmpz_init(weight);

  fmpz_t r;
  fmpz_init(r);

  fq_nmod_mat_t Mp;
  fq_nmod_mat_init(Mp, n, n, ctx);

  fq_nmod_mat_t B;
  fq_nmod_mat_init(B, n, n, ctx);

  for (slong i = 2; i <= n; ++i) {
    fmpz_add(total, total, q_pow_i);
    fmpz_mul(q_pow_i, q_pow_i, q);

    fmpz_randm(r, state, total);
    fmpz_zero(acc_weight);

    slong j = -1;
    for (slong k = 1; k <= i; ++k) {
      fmpz_pow_ui(weight, q, i - k);
      fmpz_add(acc_weight, acc_weight, weight);

      if (fmpz_cmp(r, acc_weight) < 0) {
        j = k;
        break;
      }
    }

    memmove(&sigma[i - j + 1], &sigma[i - j], (j - 1) * sizeof(slong));
    sigma[i - j] = i - 1;
  }

  _perm_inv(sigma_inv, sigma, n);

  for (slong i = 0; i < n; ++i) {
    for (slong j = 0; j < n; ++j) {
      if (j == sigma[i]) {
        fq_nmod_one(fq_nmod_mat_entry(Mp, i, j), ctx);
      } else if (j < sigma[i] || sigma_inv[j] < i) {
        fq_nmod_zero(fq_nmod_mat_entry(Mp, i, j), ctx);
      } else {
        fq_nmod_rand(fq_nmod_mat_entry(Mp, i, j), state, ctx);
      }
    }
  }

  fq_nmod_mat_randtriu(B, state, 0, ctx);
  fq_nmod_mat_transpose(B, B, ctx);

  fq_nmod_mat_mul(M, B, Mp, ctx);

  fq_nmod_mat_clear(B, ctx);
  fq_nmod_mat_clear(Mp, ctx);
  fmpz_clear(r);
  fmpz_clear(weight);
  fmpz_clear(acc_weight);
  fmpz_clear(total);
  fmpz_clear(q_pow_i);
  fmpz_clear(q);
  _perm_clear(sigma_inv);
  _perm_clear(sigma);
}
