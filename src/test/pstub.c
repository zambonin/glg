#include <assert.h>
#include <getopt.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "common.h"
#include "test/stub.h"

static const uint32_t NS_TO_SEC = 1000000000;

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

#define PERF(total_time, total_cycles, logic, var)                             \
  struct timespec var##_tstart;                                                \
  struct timespec var##_tstop;                                                 \
                                                                               \
  clock_gettime(CLOCK_MONOTONIC_RAW, &var##_tstart);                           \
  uint64_t var##_cstart = cycles();                                            \
                                                                               \
  logic;                                                                       \
                                                                               \
  uint64_t var##_cstop = cycles();                                             \
  clock_gettime(CLOCK_MONOTONIC_RAW, &var##_tstop);                            \
                                                                               \
  total_time +=                                                                \
      (long double)(var##_tstop.tv_sec - var##_tstart.tv_sec) * NS_TO_SEC +    \
      (long double)(var##_tstop.tv_nsec - var##_tstart.tv_nsec);               \
  total_cycles += var##_cstop - var##_cstart;

uint64_t cycles(void) {
  uint64_t result = 0;
  __asm volatile(".byte 15;.byte 49;shlq $32,%%rdx;orq %%rdx,%%rax"
                 : "=a"(result)::"%rdx");
  return result;
}

typedef struct {
  const fq_nmod_ctx_t *ctx;
  const ulong n;
  const ulong q;
} nq_cfg;

static fq_nmod_ctx_t snova_ctx;
static fq_nmod_ctx_t atfe_1_ctx;
static fq_nmod_ctx_t atfe_2_ctx;
static fq_nmod_ctx_t atfe_3_ctx;
static fq_nmod_ctx_t meds_1_ctx;
static fq_nmod_ctx_t meds_2_ctx;

static const nq_cfg USE_CASES[] = {
    {.ctx = &snova_ctx, .n = 2, .q = 16},
    {.ctx = &snova_ctx, .n = 3, .q = 16},
    {.ctx = &snova_ctx, .n = 4, .q = 16},
    {.ctx = &snova_ctx, .n = 5, .q = 16},
    {.ctx = &atfe_1_ctx, .n = 9, .q = 524287},
    {.ctx = &atfe_2_ctx, .n = 10, .q = 131071},
    {.ctx = &atfe_3_ctx, .n = 11, .q = 65521},
    {.ctx = &meds_1_ctx, .n = 14, .q = 4093},
    {.ctx = &meds_1_ctx, .n = 22, .q = 4093},
    {.ctx = &meds_2_ctx, .n = 30, .q = 2039},
};

void init_all_fq_nmod_ctx_t() {
  fq_nmod_ctx_init_ui(snova_ctx, 2, 4, "a");
  fq_nmod_ctx_init_ui(atfe_1_ctx, 524287, 1, "b");
  fq_nmod_ctx_init_ui(atfe_2_ctx, 131071, 1, "c");
  fq_nmod_ctx_init_ui(atfe_3_ctx, 65521, 1, "d");
  fq_nmod_ctx_init_ui(meds_1_ctx, 4093, 1, "e");
  fq_nmod_ctx_init_ui(meds_2_ctx, 2039, 1, "f");
}

void clear_all_fq_nmod_ctx_t() {
  fq_nmod_ctx_clear(snova_ctx);
  fq_nmod_ctx_clear(atfe_1_ctx);
  fq_nmod_ctx_clear(atfe_2_ctx);
  fq_nmod_ctx_clear(atfe_3_ctx);
  fq_nmod_ctx_clear(meds_1_ctx);
  fq_nmod_ctx_clear(meds_2_ctx);
}

void print_perf_stats(const fq_nmod_ctx_t ctx, const ulong n, const uint32_t it,
                      const long double time, const long double cycles) {
  const slong p = fq_nmod_ctx_prime(ctx);
  const ulong d = fq_nmod_ctx_degree(ctx);
  printf(
      "p = %8lu, d = %2lu, n = %3lu, i = %5u, avg = %14.2Lf ns, %14.2Lf cyc.\n",
      p, d, n, it, time / it, cycles / it);
}

void run_suite(uint32_t iterations, flint_rand_t state, sample_gl_func sample) {
  const slong len = sizeof(USE_CASES) / sizeof(USE_CASES[0]);
  for (uint32_t i = 0; i < len; ++i) {
    const nq_cfg *actual = &USE_CASES[i];

    const fq_nmod_ctx_t *ctx = actual->ctx;
    const ulong n = actual->n;

    long double total_time = 0;
    long double total_cycles = 0;

    fq_nmod_mat_t M;
    fq_nmod_mat_init(M, n, n, *ctx);

    fq_nmod_mat_t inv;
    fq_nmod_mat_init(inv, n, n, *ctx);

    for (uint32_t warmup = 0; warmup < iterations / 10; ++warmup) {
      fq_nmod_mat_zero(M, *ctx);

      sample(M, *ctx, state);

      assert(fq_nmod_mat_inv(inv, M, *ctx));
    }

    for (uint32_t j = 0; j < iterations; ++j) {
      fq_nmod_mat_zero(M, *ctx);

      PERF(total_time, total_cycles, sample(M, *ctx, state), p);

      assert(fq_nmod_mat_inv(inv, M, *ctx));
    }

    print_perf_stats(*ctx, n, iterations, total_time, total_cycles);

    fq_nmod_mat_clear(inv, *ctx);
    fq_nmod_mat_clear(M, *ctx);
  }
}

int32_t main_stub(int32_t argc, char **argv, sample_gl_func f) {
  uint32_t iterations = 1024;
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

  init_all_fq_nmod_ctx_t();

  run_suite(iterations, state, f);

  clear_all_fq_nmod_ctx_t();

  flint_rand_clear(state);

  flint_cleanup_master();

  return 0;
}
