#include "comb/poly.h"
#include "comb/vspace.h"

void unrank_vspace(fq_nmod_mat_t v, const fmpz_t r, const fq_nmod_ctx_t ctx) {
  fq_nmod_mat_zero(v, ctx);

  fmpz_t q;
  fmpz_init(q);
  fq_nmod_ctx_order(q, ctx);

  fmpz_t rp;
  fmpz_init_set(rp, r);

  fmpz_t rem;
  fmpz_init(rem);

  for (slong i = v->r - 1; i >= 0; --i) {
    fmpz_mod(rem, rp, q);
    unrank_poly(fq_nmod_mat_entry(v, i, 0), rem, ctx);
    fmpz_fdiv_q(rp, rp, q);
  }

  fmpz_clear(rem);
  fmpz_clear(rp);
  fmpz_clear(q);
}

void rank_vspace(fmpz_t r, const fq_nmod_mat_t v, const fq_nmod_ctx_t ctx) {
  fmpz_zero(r);

  fmpz_t q;
  fmpz_init(q);
  fq_nmod_ctx_order(q, ctx);

  fmpz_t tmp;
  fmpz_init(tmp);

  for (slong i = 0; i < v->r; ++i) {
    fmpz_mul(r, r, q);
    rank_poly(tmp, fq_nmod_mat_entry(v, i, 0), ctx);
    fmpz_add(r, r, tmp);
  }

  fmpz_clear(tmp);
  fmpz_clear(q);
}
