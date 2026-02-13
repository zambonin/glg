#include <dlfcn.h>
#include <stdint.h>
#include <stdio.h>

static const char fmt[] = ", limb = %8.2Lf, chunk = %3d";
static long double *count;

void reset_cnt(void) {
  if (count) {
    *count = 0;
  }
}

void norm_cnt(uint32_t i) {
  if (count) {
    *count /= (double)i;
  }
}

void print_cnt(void) {
  if (count) {
    printf(fmt, *count, CHUNK_SIZE);
  }
}

void find_cnt(void) {
  void *handle = dlopen(NULL, RTLD_LAZY);
  count = (long double *)dlsym(handle, "n_randlimb_count");
  dlclose(handle);
}
