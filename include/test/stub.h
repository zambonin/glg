#ifndef STUB_H
#define STUB_H

#include <stdint.h>

#include "common.h"

typedef void (*sample_gl_func)(fq_nmod_mat_t, const fq_nmod_ctx_t,
                               flint_rand_t);

int32_t main_stub(int32_t argc, char **argv, sample_gl_func f);

#endif // STUB_H
