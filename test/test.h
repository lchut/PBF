#include <omp.h>
#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <mutex>
#include <atomic>

inline void print(uint32_t* data, int n) {
    for (int i = 0; i < n; ++i) {
        printf("%d ", data[i]);
    }
    printf("\n");
}

void testPrefixSum();
void testSolver();