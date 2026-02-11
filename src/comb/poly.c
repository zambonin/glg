#include "comb/poly.h"

void unrank_poly(fq_nmod_t poly, const fmpz_t r, const fq_nmod_ctx_t ctx) {
  const ulong p = fq_nmod_ctx_prime(ctx);
  const slong d = fq_nmod_ctx_degree(ctx);
  fq_nmod_zero(poly, ctx);

  fmpz_t rp;
  fmpz_init_set(rp, r);

  for (slong j = 0; j < d; ++j) {
    if (fmpz_is_zero(rp)) {
      break;
    }
    nmod_poly_set_coeff_ui(poly, j, fmpz_fdiv_ui(rp, p));
    fmpz_fdiv_q_ui(rp, rp, p);
  }

  fmpz_clear(rp);
}

void rank_poly(fmpz_t r, const fq_nmod_t poly, const fq_nmod_ctx_t ctx) {
  const ulong p = fq_nmod_ctx_prime(ctx);
  const slong d = fq_nmod_ctx_degree(ctx);
  fmpz_zero(r);

  for (slong j = d - 1; j >= 0; --j) {
    fmpz_mul_ui(r, r, p);
    fmpz_add_ui(r, r, nmod_poly_get_coeff_ui(poly, j));
  }
}
