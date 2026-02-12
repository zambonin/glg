#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <flint/ulong_extras.h>

#include "comb/gl.h"
#include "comb/poly.h"
#include "comb/sspace.h"
#include "comb/vspace.h"
#include "common.h"
#include "util/math.h"

#define EXHAUSTIVE_LIMIT 8192
#define RANDOM_SAMPLES 1024

void test_poly(const fq_nmod_ctx_t ctx) {
  fmpz_t q;
  fmpz_init(q);
  fq_nmod_ctx_order(q, ctx);

  fmpz_t i;
  fmpz_init(i);

  fmpz_t r;
  fmpz_init(r);

  fq_nmod_t poly;
  fq_nmod_init(poly, ctx);

  for (fmpz_zero(i); fmpz_cmp(i, q) < 0; fmpz_add_ui(i, i, 1)) {
    unrank_poly(poly, i, ctx);
    rank_poly(r, poly, ctx);
    assert(fmpz_equal(i, r));
  }

  fq_nmod_clear(poly, ctx);
  fmpz_clear(r);
  fmpz_clear(i);
  fmpz_clear(q);
}

void test_vspace(const fq_nmod_ctx_t ctx, const slong n, flint_rand_t state) {
  fmpz_t q;
  fmpz_init(q);
  fq_nmod_ctx_order(q, ctx);

  fmpz_t count;
  fmpz_init(count);
  fmpz_pow_ui(count, q, n);

  fmpz_t i;
  fmpz_init(i);

  fmpz_t r;
  fmpz_init(r);

  fq_nmod_mat_t v;
  fq_nmod_mat_init(v, n, 1, ctx);

  if (fmpz_cmp_ui(count, EXHAUSTIVE_LIMIT) <= 0) {
    for (fmpz_zero(i); fmpz_cmp(i, count) < 0; fmpz_add_ui(i, i, 1)) {
      unrank_vspace(v, i, ctx);
      rank_vspace(r, v, ctx);
      assert(fmpz_equal(i, r));
    }
  } else {
    for (slong s = 0; s < RANDOM_SAMPLES; ++s) {
      fmpz_randm(i, state, count);
      unrank_vspace(v, i, ctx);
      rank_vspace(r, v, ctx);
      assert(fmpz_equal(i, r));
    }
  }

  fq_nmod_mat_clear(v, ctx);
  fmpz_clear(r);
  fmpz_clear(i);
  fmpz_clear(count);
  fmpz_clear(q);
}

void test_subspace(const fq_nmod_ctx_t ctx, const slong n, flint_rand_t state) {
  fmpz_t q;
  fmpz_init(q);
  fq_nmod_ctx_order(q, ctx);

  fmpz_t count;
  fmpz_init(count);
  fmpz_pow_ui(count, q, n);
  fmpz_sub_ui(count, count, 1);

  fmpz_t i;
  fmpz_init(i);

  fmpz_t r;
  fmpz_init(r);

  fq_nmod_mat_t v;
  fq_nmod_mat_init(v, n, 1, ctx);

  if (fmpz_cmp_ui(count, EXHAUSTIVE_LIMIT) <= 0) {
    for (fmpz_zero(i); fmpz_cmp(i, count) < 0; fmpz_add_ui(i, i, 1)) {
      unrank_subspace(v, NULL, i, ctx);
      rank_subspace(r, v, NULL, ctx);
      assert(fmpz_equal(i, r));
    }
  } else {
    for (slong s = 0; s < RANDOM_SAMPLES; ++s) {
      fmpz_randm(i, state, count);
      unrank_subspace(v, NULL, i, ctx);
      rank_subspace(r, v, NULL, ctx);
      assert(fmpz_equal(i, r));
    }
  }

  if (n >= 2) {
    fq_nmod_mat_t B;
    fq_nmod_mat_init(B, n, 1, ctx);
    fq_nmod_mat_randtest(B, state, ctx);

    if (fq_nmod_mat_is_zero(B, ctx)) {
      fq_nmod_one(fq_nmod_mat_entry(B, 0, 0), ctx);
    }

    fmpz_pow_ui(count, q, n);
    fmpz_sub(count, count, q);

    if (fmpz_cmp_ui(count, EXHAUSTIVE_LIMIT) <= 0) {
      for (fmpz_zero(i); fmpz_cmp(i, count) < 0; fmpz_add_ui(i, i, 1)) {
        unrank_subspace(v, B, i, ctx);
        rank_subspace(r, v, B, ctx);
        assert(fmpz_equal(i, r));
      }
    } else {
      for (slong s = 0; s < RANDOM_SAMPLES; ++s) {
        fmpz_randm(i, state, count);
        unrank_subspace(v, B, i, ctx);
        rank_subspace(r, v, B, ctx);
        assert(fmpz_equal(i, r));
      }
    }

    fq_nmod_mat_clear(B, ctx);
  }

  fq_nmod_mat_clear(v, ctx);
  fmpz_clear(r);
  fmpz_clear(i);
  fmpz_clear(count);
  fmpz_clear(q);
}

void test_gl(const fq_nmod_ctx_t ctx, const slong n, flint_rand_t state) {
  fmpz_t count;
  fmpz_init(count);
  gl_order(count, ctx, n);

  fmpz_t i;
  fmpz_init(i);

  fmpz_t r;
  fmpz_init(r);

  fq_nmod_mat_t M;
  fq_nmod_mat_init(M, n, n, ctx);

  if (fmpz_cmp_ui(count, EXHAUSTIVE_LIMIT) <= 0) {
    for (fmpz_zero(i); fmpz_cmp(i, count) < 0; fmpz_add_ui(i, i, 1)) {
      unrank_gl(M, i, ctx);
      rank_gl(r, M, ctx);
      assert(fmpz_equal(i, r));
    }
  } else {
    for (slong s = 0; s < RANDOM_SAMPLES; ++s) {
      fmpz_randm(i, state, count);
      unrank_gl(M, i, ctx);
      rank_gl(r, M, ctx);
      assert(fmpz_equal(i, r));
    }
  }

  fq_nmod_mat_clear(M, ctx);
  fmpz_clear(r);
  fmpz_clear(i);
  fmpz_clear(count);
}

void test_complete_basis(const fq_nmod_ctx_t ctx, const slong n,
                         flint_rand_t state) {
  for (slong k = 1; k < n; ++k) {
    fq_nmod_mat_t S;
    fq_nmod_mat_init(S, n, k, ctx);

    fq_nmod_mat_t tmp;
    fq_nmod_mat_init(tmp, n, n, ctx);
    fq_nmod_mat_randrank(tmp, state, n, ctx);

    fq_nmod_mat_t view;
    fq_nmod_mat_window_init(view, tmp, 0, 0, n, k, ctx);
    fq_nmod_mat_set(S, view, ctx);
    fq_nmod_mat_window_clear(view, ctx);
    fq_nmod_mat_clear(tmp, ctx);

    fq_nmod_mat_t B;
    fq_nmod_mat_init(B, n, n, ctx);
    complete_basis(B, S, n, ctx);

    fq_nmod_mat_t inv;
    fq_nmod_mat_init(inv, n, n, ctx);
    assert(fq_nmod_mat_inv(inv, B, ctx));

    fq_nmod_mat_clear(inv, ctx);
    fq_nmod_mat_clear(B, ctx);
    fq_nmod_mat_clear(S, ctx);
  }
}

int32_t main(void) {
  flint_rand_t state;
  flint_rand_init(state);
  state->__randval = (uint32_t)time(NULL);

  const ulong configs[][3] = {
      {2, 1, 2},
      {2, 1, 3},
      {2, 1, 4},
      {3, 1, 2},
      {3, 1, 3},
      {2, 2, 2},
      {2, 1, 8},
      {3, 1, 6},
      {5, 2, 4},
      {7, 1, 5},
      {11, 1, 4},
  };

  const slong num_configs = sizeof(configs) / sizeof(configs[0]);

  for (slong c = 0; c < num_configs; ++c) {
    const ulong p = configs[c][0];
    const ulong d = configs[c][1];
    const ulong n = configs[c][2];

    fq_nmod_ctx_t ctx;
    fq_nmod_ctx_init_ui(ctx, p, d, "x");

    printf("p = %3lu, d = %2lu, n = %2lu\n", p, d, n);

    test_poly(ctx);
    test_vspace(ctx, n, state);
    test_subspace(ctx, n, state);
    test_gl(ctx, n, state);
    test_complete_basis(ctx, n, state);

    fq_nmod_ctx_clear(ctx);
  }

  flint_rand_clear(state);
  flint_cleanup_master();

  return 0;
}
