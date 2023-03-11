#ifndef COUNTVF_H
#define COUNTVF_H


#include<immintrin.h>
#include<pmmintrin.h>
//#include<zmmintrin.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <stdint.h>
#include<string>
#include <fstream>
#include <sstream>
#include<vector>
#include<time.h>
#include <chrono>
#include <algorithm>
#include<random>

#define _mm256_neg_epi32(x) _mm256_sub_epi32(MAX, x)

#define __int32 int
#define _mm256_and_epi32(x, y) _mm256_and_si256(x, y)
#define _mm256_storeu_epi32(a, b)  _mm256_storeu_si256((__m256i*)a, b)

const __m256i MAX = _mm256_set1_epi32(0xFFFFFFFF);
const __m256i ZERO = _mm256_setzero_si256();
const __m256i LEFT1 = _mm256_setr_epi32(8, 0, 1, 2, 3, 4, 5, 6);
const __m256i LEFT2 = _mm256_setr_epi32(8, 8, 0, 1, 2, 3, 4, 5);
const __m256i LEFT4 = _mm256_setr_epi32(8, 8, 8, 8, 0, 1, 2, 3);
const __m256i POS7 = _mm256_set1_epi32(7);
const __m256i NEG1 = _mm256_set1_epi32(-1);
const __m256i POS1 = _mm256_set1_epi32(1);
const __m256i MASK0 = _mm256_setr_epi32(
    0x0, 0xFFFFFFFF,
    0xFFFFFFFF, 0xFFFFFFFF,
    0xFFFFFFFF, 0xFFFFFFFF,
    0xFFFFFFFF, 0xFFFFFFFF
);

const __m256i MASK00 = _mm256_setr_epi32(
    0x0, 0x0,
    0xFFFFFFFF, 0xFFFFFFFF,
    0xFFFFFFFF, 0xFFFFFFFF,
    0xFFFFFFFF, 0xFFFFFFFF
);

const __m256i MASK0000 = _mm256_setr_epi32(
    0x0, 0x0,
    0x0, 0x0,
    0xFFFFFFFF, 0xFFFFFFFF,
    0xFFFFFFFF, 0xFFFFFFFF
);
 
int QUALIFY = 10000;

__m256i permute1(__m256i data);

__m256i permute2(__m256i data);

__m256i permute4(__m256i data);

__int32 sumAll(__m256i data);

uint64_t timeSinceEpochMillisec() {
    using namespace std::chrono;
    return duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();
}

int random_i32() {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist6(1, 6);
    return dist6(rng);
}

bool predicate_scalar(int val, char op, int c) {

    switch (op)
    {
    case '>':
        return val > c;
        break;
    case '<':
        return val < c;
        break;
    case '[':
        return val >= c;
        break;
    case ']':
        return val <= c;
        break;
    case '=':
        return val == c;
        break;
    default:
        return false;
        break;
    }
}

int log2(int x) {
    int ans = 0;
    while (x) {
        x >>= 1;
        ans++;
    }
    return ans;
}

__m256i predicates(__m256i vd, __m256i vcmp, char op) {
    __m256i flag;
    switch (op)
    {
    case '>':
        flag = _mm256_cmpgt_epi32(vd, vcmp);
        break;
    case '<':
        flag = _mm256_cmpgt_epi32(vcmp, vd);
        break;
    case '[':
        flag = _mm256_cmpgt_epi32(vcmp, vd);
        flag = _mm256_sub_epi32(MAX, flag);
        break;
    case ']':
        flag = _mm256_cmpgt_epi32(vd, vcmp);
        flag = _mm256_sub_epi32(MAX, flag);
        break;
    case '=':
        flag = _mm256_cmpeq_epi32(vd, vcmp);
        break;
    default:
        flag = MAX;
        break;
    }
    return flag;
}

uint32_t hsum_epi32_avx(__m128i x)
{
    __m128i hi64 = _mm_unpackhi_epi64(x, x);           // 3-operand non-destructive AVX lets us save a byte without needing a movdqa
    __m128i sum64 = _mm_add_epi32(hi64, x);
    __m128i hi32 = _mm_shuffle_epi32(sum64, _MM_SHUFFLE(2, 3, 0, 1));    // Swap the low two elements
    __m128i sum32 = _mm_add_epi32(sum64, hi32);
    return _mm_cvtsi128_si32(sum32);       // movd
}


uint32_t hsum_8x32(__m256i v)
{
    __m128i sum128 = _mm_add_epi32(
        _mm256_castsi256_si128(v),
        _mm256_extracti128_si256(v, 1));
    // silly GCC uses a longer AXV512VL instruction if AVX512 is enabled :/
    return hsum_epi32_avx(sum128);
}

__int32 sumAll(__m256i data) {
    // not a hotpot
    const __m128i hi = _mm256_extractf128_si256(data, 1);
    const __m128i lo = _mm256_castsi256_si128(data);
    const __m128i sum1 = _mm_add_epi32(lo, hi);
    return hsum_epi32_avx(sum1);
    /*int t[4];
    _mm_storeu_epi32(t, sum1);
    return t[0] + t[1] + t[2] + t[3];*/
}

__m256i permute4(__m256i data) {
    __m256i shift = _mm256_permutevar8x32_epi32(data, _mm256_setr_epi32(0, 0, 0, 0, 0, 1, 2, 3));
    return _mm256_and_epi32(shift, MASK0000);
}

__m256i permute2(__m256i data) {
    __m256i shift = _mm256_permutevar8x32_epi32(data, _mm256_setr_epi32(0, 0, 0, 1, 2, 3, 4, 5));
    return _mm256_and_epi32(shift, MASK00);
}

__m256i permute1(__m256i data) {
    __m256i shift = _mm256_permutevar8x32_epi32(data, _mm256_setr_epi32(0, 0, 1, 2, 3, 4, 5, 6));
    return _mm256_and_epi32(shift, MASK0);
}

void see_into(__m256i res);

__m256i rle_count(
    __m256i diff1ST, __m256i diff2ST,
    __m256i diff1ED, __m256i diff2ED,
    __m256i delta, __m256i rle,
    char op1, char op2
);

__m256i count_filter(
    int* deltas, int* rles,
    char op1, char op2,
    __m256i& vc1, __m256i& vc2
);


#endif