#include <time.h>

#define STOP 0
#define START 1

void chrono (int kind, double *time) {
  static clock_t counts;
  if (kind == START) {
    *time = 0.0;
    counts = clock();
    return;
  }
  if (kind == STOP) {
    *time = ((double)(clock()-counts))/((double)CLOCKS_PER_SEC);
  }
}

