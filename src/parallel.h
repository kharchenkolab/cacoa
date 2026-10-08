// Thin wrapper around sccore's std::thread-based parallel loop.
// sccore_par.hpp defines non-inline functions, so it may be included in exactly one translation
// unit (parallel.cpp); every other file uses this declaration.
#ifndef CACOA_PARALLEL_H
#define CACOA_PARALLEL_H

#include <functional>

namespace cacoa {
// Runs task(i) for i in [start, end) on n_cores threads (serially when n_cores <= 1).
// Tasks must be independent; with verbose = true a progress bar is printed from the master thread.
void parallelFor(int start, int end, const std::function<void(int)> &task, int n_cores, bool verbose = false);
}

#endif
