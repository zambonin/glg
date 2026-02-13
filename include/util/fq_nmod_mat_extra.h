#ifndef FQ_NMOD_MAT_EXTRA_H
#define FQ_NMOD_MAT_EXTRA_H

#include "common.h"

struct bit_buffer_struct;
typedef struct bit_buffer_struct bit_buffer_t;

void fq_nmod_mat_randtest_not_zero_buf(fq_nmod_mat_t M, const fq_nmod_ctx_t ctx,
                                       bit_buffer_t *buf);

slong fq_nmod_mat_first_pos_entry(fq_nmod_mat_t v, const fq_nmod_ctx_t ctx);

void fq_nmod_mat_set_minor(fq_nmod_mat_t mat, const ulong skip_i,
                           const ulong skip_j, const fq_nmod_mat_t src,
                           const fq_nmod_ctx_t ctx);

void fq_nmod_mat_scalar_mul(fq_nmod_mat_t B, const fq_nmod_mat_t A,
                            const fq_nmod_t c, const fq_nmod_ctx_t ctx);

void fq_nmod_mat_entry_rand_buf(fq_nmod_mat_t M, const slong i, const slong j,
                                const fq_nmod_ctx_t ctx, bit_buffer_t *buf);

void fq_nmod_mat_entry_rand_not_zero_buf(fq_nmod_mat_t M, const slong i,
                                         const slong j,
                                         const fq_nmod_ctx_t ctx,
                                         bit_buffer_t *buf);

void fq_nmod_mat_randtest_buf(fq_nmod_mat_t M, const fq_nmod_ctx_t ctx,
                              bit_buffer_t *buf);

void fq_nmod_mat_randtriu_buf(fq_nmod_mat_t M, const fq_nmod_ctx_t ctx,
                              bit_buffer_t *buf, const int unit);

void fq_nmod_mat_randtril_buf(fq_nmod_mat_t M, const fq_nmod_ctx_t ctx,
                              bit_buffer_t *buf, const int unit);

#endif // FQ_NMOD_MAT_EXTRA_H
