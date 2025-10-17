#include <dlfcn.h>
#include <stdint.h>

#include <flint/flint.h>
#include <flint/fq_nmod.h>
#include <gmp.h>

static mp_limb_t (*n_randlimb_func)(flint_rand_t) = NULL;

long double n_randlimb_count = 0;

mp_limb_t n_randlimb(flint_rand_t state) {
  if (!n_randlimb_func) {
    n_randlimb_func = dlsym(RTLD_NEXT, "n_randlimb");
  }

  ++n_randlimb_count;
  return n_randlimb_func(state);
}
