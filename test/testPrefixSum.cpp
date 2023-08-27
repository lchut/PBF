#include "test.h"

void computePrefixSum(uint32_t* data, uint32_t* prefixSum, int n) {
    auto start = std::chrono::system_clock::now();
    prefixSum[0] = 0;
    for (int i = 1; i < n; ++i) {
        prefixSum[i] = prefixSum[i - 1] + data[i - 1];
    }
    auto end = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den << std::endl;
}

void computePrefixSumAt(std::atomic<uint32_t>* data, std::atomic<uint32_t>* prefixSum, int n) {
    auto start = std::chrono::system_clock::now();
    prefixSum[0] = 0;
    for (int i = 1; i < n; ++i) {
        prefixSum[i] = prefixSum[i - 1] + data[i - 1];
    }
    auto end = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den << std::endl;
}


void computeParllelPrefixSum(uint32_t* data, uint32_t* prefixSum, int n, int pow2GreaterN) {
    auto start = std::chrono::system_clock::now();
    #pragma parallel for 
    for (int i = 0; i < n; ++i) { prefixSum[i] = data[i]; }
    
    // down-swap
    for (int d = 0, stride = 1; d < int(std::log2(pow2GreaterN)); ++d, stride *= 2) {
        int stride2 = stride * 2;
        #pragma parallel for
        for (int k = 0; k < pow2GreaterN; k += stride2) {
            prefixSum[k + stride2 - 1] += prefixSum[k + stride - 1];
        }
    }
    prefixSum[pow2GreaterN - 1] = 0;
    // up-swap
    for (int d = int(std::log2(pow2GreaterN)) - 1, stride = std::pow(2, int(std::log2(pow2GreaterN))); d >= 0; d--, stride = (stride >> 1)) {
        int halfStride = (stride >> 1);
        #pragma parallel for 
        for (int k = 0; k < pow2GreaterN; k += stride) {
            uint32_t tmp = prefixSum[k + halfStride - 1];
            prefixSum[k + halfStride - 1] = prefixSum[k + stride - 1];
            prefixSum[k + stride - 1] += tmp;
        }
    }
    auto end = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den << std::endl;
}

void check(uint32_t* prefixSum1, uint32_t* prefixSum2, int n) {
    for (int i = 0; i < n; i++) {
        if (prefixSum1[i] != prefixSum2[i]) {
            std::cout << "Wrong Answer" << std::endl;
            return;
        }
    }
    printf("Done!\n");
}

void testPrefixSum() {
    int N = 64 * 64 * 64;
    uint32_t* data = new uint32_t[N];
    std::atomic<uint32_t>* atData = new std::atomic<uint32_t>[N];
    uint32_t pow2GreaterN = 1;
    while (pow2GreaterN < N) { pow2GreaterN = (pow2GreaterN << 1); }
    uint32_t* prefixSum1 = new uint32_t[pow2GreaterN];
    uint32_t* prefixSum2 = new uint32_t[pow2GreaterN];
    std::atomic<uint32_t>* prefixSumat = new std::atomic<uint32_t>[pow2GreaterN];
    std::srand((unsigned(time(NULL))));
    for (int i = 0; i < N; ++i) {
        data[i] = rand() % 16;
        atData[i] = data[i];
    }
    computePrefixSum(data, prefixSum1, N);
    computePrefixSumAt(atData, prefixSumat, N);
    computeParllelPrefixSum(data, prefixSum2, N, pow2GreaterN);
    //print(data, N);
    //print(prefixSum1, N);
    //print(prefixSum2, N);
    check(prefixSum1, prefixSum2, N);
}