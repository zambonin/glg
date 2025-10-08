#include <assert.h>
#include <getopt.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <flint/ulong_extras.h>

#include "test/stub.h"

static const struct option long_options[] = {
    {"iterations", required_argument, 0, 'i'},
    {"seed", required_argument, 0, 's'},
    {0, 0, 0, 0}};

static const char *help_text =
    "Usage: %s [OPTIONS]\n"
    "  -i, --iterations=<uint32_t>\n"
    "         Set the number of random tests per configuration.\n"
    "\n"
    "  -s, --seed=<uint32_t>\n"
    "         Set the seed for the random number generator.\n";

void run_suite(uint32_t iterations, flint_rand_t state, sample_gl_func sample) {
  for (uint32_t i = 0; i < iterations; ++i) {
    fq_nmod_ctx_t ctx;
    fq_nmod_ctx_init_randtest(ctx, state, 3);

    const slong p = fq_nmod_ctx_prime(ctx);
    const ulong d = fq_nmod_ctx_degree(ctx);
    const ulong n = 1 + n_urandint(state, 1 << 5);

    printf("p = %3ld, d = %3ld, n = %3ld\n", p, d, n);

    fq_nmod_mat_t M;
    fq_nmod_mat_init(M, n, n, ctx);

    fq_nmod_mat_t inv;
    fq_nmod_mat_init(inv, n, n, ctx);

    for (uint32_t j = 0; j < iterations; ++j) {
      fq_nmod_mat_zero(M, ctx);

      sample(M, ctx, state);

      assert(fq_nmod_mat_inv(inv, M, ctx));
    }

    fq_nmod_mat_clear(inv, ctx);
    fq_nmod_mat_clear(M, ctx);
    fq_nmod_ctx_clear(ctx);
  }
}

int32_t main_stub(int32_t argc, char **argv, sample_gl_func f) {
  uint32_t iterations = 16;
  uint32_t seed = (uint32_t)time(NULL);

  for (;;) {
    int c = getopt_long(argc, argv, "i:s:", long_options, NULL);
    if (c == -1) {
      break;
    }

    switch (c) {
    case 'i':
      iterations = strtoul(optarg, NULL, 10);
      break;
    case 's':
      seed = strtoul(optarg, NULL, 10);
      break;
    default:
      fprintf(stderr, help_text, argv[0]);
      return 1;
    }
  }

  flint_rand_t state;
  flint_rand_init(state);

  state->__randval = seed;

  run_suite(iterations, state, f);

  flint_rand_clear(state);

  flint_cleanup_master();

  return 0;
}
