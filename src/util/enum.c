#include <getopt.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "comb/gl.h"
#include "common.h"
#include "util/math.h"
#include "comb/poly.h"

static const struct option long_options[] = {
    {"characteristic", required_argument, 0, 'p'},
    {"dimension", required_argument, 0, 'd'},
    {"degree", required_argument, 0, 'n'},
    {0, 0, 0, 0}};

static const char *help_text =
    "Usage: %s [OPTIONS]\n"
    "  -p, --characteristic=<uint8_t>\n"
    "         Set the order of the base finite field.\n"
    "\n"
    "  -d, --dimension=<uint8_t>\n"
    "         Set the dimension of the extension field.\n"
    "\n"
    "  -n, --degree=<uint8_t>\n"
    "         Set the degree of the general linear group.\n";

static void matrix_print_oneline(const fq_nmod_mat_t M,
                                 const fq_nmod_ctx_t ctx) {
  fmpz_t r;
  fmpz_init(r);
  for (ulong i = 0; i < M->r; ++i) {
    for (ulong j = 0; j < M->c; ++j) {
      rank_poly(r, fq_nmod_mat_entry(M, i, j), ctx);
      fmpz_print(r);
      printf(" ");
    }
  }
  fmpz_clear(r);
  printf("\n");
}

int32_t main(int32_t argc, char **argv) {
  uint8_t p = 0;
  uint8_t d = 0;
  uint8_t n = 0;

  for (;;) {
    int c = getopt_long(argc, argv, "p:d:n:", long_options, NULL);
    if (c == -1) {
      break;
    }

    switch (c) {
    case 'p':
      p = strtol(optarg, NULL, 10);
      break;
    case 'd':
      d = strtol(optarg, NULL, 10);
      break;
    case 'n':
      n = strtol(optarg, NULL, 10);
      break;
    default:
      fprintf(stderr, help_text, argv[0]);
      return 1;
    }
  }

  if (p <= 1 || d <= 0 || n <= 0) {
    fprintf(stderr, help_text, argv[0]);
    return 1;
  }

  fq_nmod_ctx_t ctx;
  fq_nmod_ctx_init_ui(ctx, p, d, "x");

  fq_nmod_mat_t M;
  fq_nmod_mat_init(M, n, n, ctx);

  fmpz_t order;
  fmpz_init(order);
  gl_order(order, ctx, n);

  fmpz_t i;
  fmpz_init(i);

  for (fmpz_zero(i); fmpz_cmp(i, order) < 0; fmpz_add_ui(i, i, 1)) {
    unrank_gl(M, i, ctx);
    matrix_print_oneline(M, ctx);
  }

  fmpz_clear(i);
  fmpz_clear(order);
  fq_nmod_mat_clear(M, ctx);
  fq_nmod_ctx_clear(ctx);

  return 0;
}
