#include "util/rand.h"
#include "comb/poly.h"

#ifndef CHUNK_SIZE
#define CHUNK_SIZE FLINT_BITS
#endif

_Static_assert(CHUNK_SIZE >= 0 && CHUNK_SIZE <= FLINT_BITS,
               "CHUNK_SIZE must be between 0 and FLINT_BITS");

int fmpz_is_pow2(const fmpz_t n) {
  if (fmpz_sgn(n) <= 0) {
    return 0;
  }
  return fmpz_sizeinbase(n, 2) == fmpz_val2(n) + 1;
}

void bit_buffer_init(bit_buffer_t *buf, flint_rand_t state) {
  fmpz_init(buf->buffer);
  buf->num_bits = 0;
  buf->bits_consumed = 0;
  buf->state = state;
}

void bit_buffer_clear(bit_buffer_t *buf) { fmpz_clear(buf->buffer); }

void bit_buffer_get_bits(fmpz_t result, bit_buffer_t *buf, const slong k) {
  if (k == 0) {
    fmpz_zero(result);
    return;
  }

  while (buf->num_bits < k) {
    const ulong limb = n_randlimb(buf->state);
    fmpz_mul_2exp(buf->buffer, buf->buffer, FLINT_BITS);
    fmpz_add_ui(buf->buffer, buf->buffer, limb);
    buf->num_bits += FLINT_BITS;
  }

  fmpz_tdiv_r_2exp(result, buf->buffer, k);
  fmpz_tdiv_q_2exp(buf->buffer, buf->buffer, k);
  buf->num_bits -= k;
  buf->bits_consumed += k;
}

void bit_buffer_discard(bit_buffer_t *buf) {
  slong to_discard = CHUNK_SIZE;
  if (to_discard > buf->num_bits) {
    to_discard = buf->num_bits;
  }

  if (to_discard > 0) {
    fmpz_tdiv_q_2exp(buf->buffer, buf->buffer, to_discard);
    buf->num_bits -= to_discard;
  }
}

void fast_dice_roller(fmpz_t result, const fmpz_t n, bit_buffer_t *buf) {
  if (fmpz_cmp_ui(n, 1) <= 0) {
    fmpz_zero(result);
    return;
  }

  if (fmpz_is_pow2(n)) {
    const slong k = fmpz_sizeinbase(n, 2) - 1;
    bit_buffer_get_bits(result, buf, k);
  } else {
    fmpz_t v, c, t;
    fmpz_init(v);
    fmpz_init_set_ui(c, 1);
    fmpz_init(t);

    fmpz_zero(v);

    for (;;) {
      fmpz_mul_2exp(v, v, 1);
      bit_buffer_get_bits(t, buf, 1);
      fmpz_add(v, v, t);

      fmpz_mul_2exp(c, c, 1);

      if (fmpz_cmp(c, n) >= 0) {
        if (fmpz_cmp(v, n) < 0) {
          fmpz_set(result, v);
          break;
        } else {
          fmpz_sub(v, v, n);
          fmpz_sub(c, c, n);
        }
      }
    }

    fmpz_clear(t);
    fmpz_clear(c);
    fmpz_clear(v);
  }

  bit_buffer_discard(buf);
}

void fq_nmod_rand_buf(fq_nmod_t out, const fq_nmod_ctx_t ctx,
                      bit_buffer_t *buf) {
  fmpz_t q;
  fmpz_init(q);
  fq_nmod_ctx_order(q, ctx);

  fmpz_t sample;
  fmpz_init(sample);

  fast_dice_roller(sample, q, buf);
  unrank_poly(out, sample, ctx);

  fmpz_clear(sample);
  fmpz_clear(q);
}

void fq_nmod_rand_not_zero_buf(fq_nmod_t out, const fq_nmod_ctx_t ctx,
                               bit_buffer_t *buf) {
  fmpz_t q;
  fmpz_init(q);
  fq_nmod_ctx_order(q, ctx);
  fmpz_sub_ui(q, q, 1);

  fmpz_t sample;
  fmpz_init(sample);

  fast_dice_roller(sample, q, buf);
  fmpz_add_ui(sample, sample, 1);
  unrank_poly(out, sample, ctx);

  fmpz_clear(sample);
  fmpz_clear(q);
}
