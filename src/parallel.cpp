#include "parallel.h"

#include <Rcpp.h>
#include <sccore_par.hpp>

namespace cacoa {
void parallelFor(int start, int end, const std::function<void(int)> &task, int n_cores, bool verbose) {
  if (n_cores <= 1) {
    for (int i = start; i < end; ++i) task(i);
    return;
  }
  sccore::runTaskParallelFor(start, end, task, n_cores, verbose);
}
}
