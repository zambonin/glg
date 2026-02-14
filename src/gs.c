#include "gs.h"
#include "comb/poly.h"
#include "util/fq_nmod_mat_extra.h"
#include "util/math.h"
#include "util/rand.h"

void gram_schmidt(fq_nmod_mat_t M, const fmpz_t r, const fq_nmod_ctx_t ctx) {
  const slong n = fq_nmod_mat_nrows(M, ctx);

  fmpz_t rp;
  fmpz_init_set(rp, r);

  fmpz_t q;
  fmpz_init(q);
  fq_nmod_ctx_order(q, ctx);

  fmpz_t q_pow_n;
  fmpz_init(q_pow_n);
  fmpz_pow_ui(q_pow_n, q, n);

  fmpz_t q_pow_i;
  fmpz_init_set_ui(q_pow_i, 1);

  fmpz_t k;
  fmpz_init(k);

  fmpz_t count;
  fmpz_init(count);

  fmpz_t a;
  fmpz_init(a);

  fmpz_t b;
  fmpz_init(b);

  fmpz_t rem;
  fmpz_init(rem);

  fmpz_t temp_rank;
  fmpz_init(temp_rank);

  fq_nmod_t coeff;
  fq_nmod_init(coeff, ctx);

  fq_nmod_t tmp;
  fq_nmod_init(tmp, ctx);

  slong *J = flint_malloc(n * sizeof(slong));
  for (slong i = 0; i < n; ++i) {
    J[i] = i;
  }

  fq_nmod_mat_zero(M, ctx);

  for (slong i = 0; i < n; ++i) {
    fmpz_sub(count, q_pow_n, q_pow_i);
    fmpz_fdiv_qr(rp, k, rp, count);

    fmpz_fdiv_qr(b, a, k, q_pow_i);
    fmpz_add_ui(b, b, 1);

    fmpz_set(temp_rank, a);
    for (slong j = 0; j < i; ++j) {
      fmpz_fdiv_qr(temp_rank, rem, temp_rank, q);
      unrank_poly(coeff, rem, ctx);
      if (fq_nmod_is_zero(coeff, ctx)) {
        continue;
      }

      for (slong r = 0; r < n; ++r) {
        fq_nmod_mul(tmp, fq_nmod_mat_entry(M, r, j), coeff, ctx);
        fq_nmod_add(fq_nmod_mat_entry(M, r, i), fq_nmod_mat_entry(M, r, i), tmp,
                    ctx);
      }
    }

    fmpz_set(temp_rank, b);
    slong last = -1;

    for (slong j = 0; j < n - i; ++j) {
      fmpz_fdiv_qr(temp_rank, rem, temp_rank, q);
      unrank_poly(coeff, rem, ctx);
      if (fq_nmod_is_zero(coeff, ctx)) {
        continue;
      }

      if (last == -1) {
        last = j;
      }
      fq_nmod_add(fq_nmod_mat_entry(M, J[j], i), fq_nmod_mat_entry(M, J[j], i),
                  coeff, ctx);
    }

    if (last != -1) {
      for (slong k = last; k < n - i - 1; ++k) {
        J[k] = J[k + 1];
      }
    }

    fmpz_mul(q_pow_i, q_pow_i, q);
  }

  flint_free(J);
  fq_nmod_clear(tmp, ctx);
  fq_nmod_clear(coeff, ctx);
  fmpz_clear(temp_rank);
  fmpz_clear(rem);
  fmpz_clear(b);
  fmpz_clear(a);
  fmpz_clear(count);
  fmpz_clear(k);
  fmpz_clear(q_pow_i);
  fmpz_clear(q_pow_n);
  fmpz_clear(q);
  fmpz_clear(rp);
}

void gs(fq_nmod_mat_t M, const fq_nmod_ctx_t ctx, flint_rand_t state) {
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

  gram_schmidt(M, rank, ctx);

  fmpz_clear(rank);
  fmpz_clear(count);
  bit_buffer_clear(&buf);
}
