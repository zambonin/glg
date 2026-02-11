#include <string.h>

#include <flint/perm.h>

#include "bruhat.h"
#include "util/fq_nmod_mat_extra.h"

void bruhat(fq_nmod_mat_t M, const fq_nmod_ctx_t ctx, flint_rand_t state) {
  const slong n = fq_nmod_mat_nrows(M, ctx);
  if (n == 1) {
    fq_nmod_rand_not_zero(fq_nmod_mat_entry(M, 0, 0), state, ctx);
    return;
  }

  slong *sigma = _perm_init(n);
  sigma[0] = 0;

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

  fq_nmod_mat_t D;
  fq_nmod_mat_init(D, n, n, ctx);

  fq_nmod_mat_t B1;
  fq_nmod_mat_init(B1, n, n, ctx);

  fq_nmod_mat_t B2;
  fq_nmod_mat_init(B2, n, n, ctx);

  fq_nmod_mat_t Mp;
  fq_nmod_mat_init(Mp, n, n, ctx);

  for (slong i = 2; i <= n; ++i) {
    fmpz_add(total, total, q_pow_i);
    fmpz_mul(q_pow_i, q_pow_i, q);

    fmpz_randlimb_m(r, state, total);
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

  fq_nmod_mat_zero(D, ctx);
  for (slong i = 0; i < n; i++) {
    fq_nmod_one(fq_nmod_mat_entry(D, i, sigma[i]), ctx);
  }

  fq_nmod_mat_randtriu(B1, state, 0, ctx);
  fq_nmod_mat_transpose(B1, B1, ctx);

  fq_nmod_mat_randtriu(B2, state, 0, ctx);

  fq_nmod_mat_mul(Mp, B1, D, ctx);

  fq_nmod_mat_mul(M, Mp, B2, ctx);

  fq_nmod_mat_clear(Mp, ctx);
  fq_nmod_mat_clear(B2, ctx);
  fq_nmod_mat_clear(B1, ctx);
  fq_nmod_mat_clear(D, ctx);
  fmpz_clear(r);
  fmpz_clear(weight);
  fmpz_clear(acc_weight);
  fmpz_clear(total);
  fmpz_clear(q_pow_i);
  fmpz_clear(q);
  _perm_clear(sigma);
}
