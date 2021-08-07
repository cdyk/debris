#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cassert>
#include <chrono>
#include <algorithm>
#include <cstring>
#include <cmath>

#include <vector>
#include <array>

void runNop(int32_t* data, size_t N)
{
}

void runInsertion(int32_t* data, size_t N)
{
  for (size_t i = 1; i < N; i++) {
    size_t j = i;
    const int32_t a = data[j];
    while (j && a < data[j - 1]) {
      data[j] = data[j - 1];
      j = j - 1;
    }
    data[j] = a;
  }
}

void runStdSort(int32_t* data, size_t N)
{
  std::sort(data, data + N);
}

void runQsort(int32_t* data, size_t N)
{
  std::qsort(data, N, sizeof(int32_t), [](const void* a, const void* b)
             {
               int aa = *(const int32_t*)a;
               int bb = *(const int32_t*)b;
               return aa < bb ? -1 : (aa == bb ? 0 : 1);
             });
}

__declspec(noinline) void run(int32_t* tmp, const int32_t* src, size_t N, void (*f)(int32_t*, size_t), size_t its)
{
  for (size_t i = 0; i < its; i++) {
    std::memcpy(tmp, src, sizeof(int32_t) * N);
    f(tmp, N);
  }
}


int main(int argc, char** argv)
{
  const size_t its = 100000;

  const std::array<int, 7> seeds = { 0, 1, 11, 42, 100, 12257, 1013546,  };

  std::vector<int32_t> A;
  std::vector<int32_t> B;
  std::vector<int32_t> C;

  fprintf(stdout, "Average over %zu iterations\n", its*seeds.size());
  fprintf(stdout, "| n   | overhead |     insertion     |     std::sort    |       qsort      |\n");
  fprintf(stdout, "|-----|----------|-------------------|------------------|------------------|\n");
  for (size_t N = 1; N <= 100; N++) {
    A.resize(N);
    B.resize(N);
    C.resize(N);


    long long overhead = 0;
    long long insertion = 0;
    long long stdsort = 0;
    long long qsort = 0;
    for (const auto seed : seeds) {

      if (seed == 0) {  // sorted data
        for (size_t i = 0; i < N; i++) {
          A[i] = (int32_t)i;
        }
      }
      else if (seed == 1) { // reverse sorted data
        std::srand(seed);
        for (size_t i = 0; i < N; i++) {
          A[i] = (int32_t)(N-i);
        }
      }
      else {  // random data
        std::srand(seed);
        for (size_t i = 0; i < N; i++) {
          A[i] = std::rand();
        }
      }


      // Verify sanity
      std::memcpy(C.data(), A.data(), sizeof(int32_t) * N);
      std::sort(C.begin(), C.end());
      run(B.data(), A.data(), N, runInsertion, 1);
      for (size_t i = 0; i < N; i++) assert(B[i] == C[i]);
      run(B.data(), A.data(), N, runStdSort, 1);
      for (size_t i = 0; i < N; i++) assert(B[i] == C[i]);

      run(B.data(), A.data(), N, runNop, 10);
      auto a = std::chrono::steady_clock::now();
      run(B.data(), A.data(), N, runNop, its);
      overhead += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - a).count();

      run(B.data(), A.data(), N, runInsertion, 10);
      a = std::chrono::steady_clock::now();
      run(B.data(), A.data(), N, runInsertion, its);
      insertion += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - a).count();

      run(B.data(), A.data(), N, runStdSort, 10);
      a = std::chrono::steady_clock::now();
      run(B.data(), A.data(), N, runStdSort, its);
      stdsort += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - a).count();

      run(B.data(), A.data(), N, runQsort, 10);
      a = std::chrono::steady_clock::now();
      run(B.data(), A.data(), N, runQsort, its);
      qsort += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - a).count();

    }

    double M = (double)(its * seeds.size());

    fprintf(stdout, "| %3zu |   %.1fns  | %6.1fns (C=%.2f) | %5.1fns (C=%.2f) | %6.1fns (C=%.2f) |\n",
            N, 
            overhead / M,
            insertion / M, (double)insertion / ((double)N * (double)N * M),
            stdsort / M, (double)stdsort / (std::log2((double)N) * (double)N * M),
            qsort / M, (double)qsort / (std::log2((double)N) * (double)N * M));
  }

  fprintf(stdout, "Done.\n");
  return 0;
}

