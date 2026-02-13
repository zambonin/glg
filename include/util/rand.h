#ifndef RAND_H
#define RAND_H

#include "common.h"
#include <flint/flint.h>

struct bit_buffer_struct {
  fmpz_t buffer;
  slong num_bits;
  slong bits_consumed;
  flint_rand_struct *state;
};

typedef struct bit_buffer_struct bit_buffer_t;

void bit_buffer_init(bit_buffer_t *buf, flint_rand_t state);

void bit_buffer_clear(bit_buffer_t *buf);

void bit_buffer_get_bits(fmpz_t result, bit_buffer_t *buf, const slong k);

void fast_dice_roller(fmpz_t result, const fmpz_t n, bit_buffer_t *buf);

void fq_nmod_rand_buf(fq_nmod_t out, const fq_nmod_ctx_t ctx,
                      bit_buffer_t *buf);

void fq_nmod_rand_not_zero_buf(fq_nmod_t out, const fq_nmod_ctx_t ctx,
                               bit_buffer_t *buf);

#endif // RAND_H
