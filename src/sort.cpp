#include <cstdio>
#include <chrono>
#include <vector>
#include <algorithm>

template<typename T>
__declspec(noinline) void stlSort(T* data, size_t N, bool (*lessThan)(const T&, const T&))
{
  std::sort(data, data + N, lessThan);
}


namespace {

  bool lessThan(const unsigned& a, const unsigned& b)
  {
    return a < b;
  }

}


int main(int argc, char** argv)
{
  std::srand(42);
  std::vector<unsigned> values;
  std::vector<std::vector<unsigned>> tmp;

  unsigned its = 10;
  tmp.resize(its);

  size_t N = 1;
  for (unsigned i = 0; i < 6; i++) {

    auto a = values.size();
    values.resize(N);
    for (auto k = a; k < N; k++) {
      values[k] = std::rand();
    }

    for (unsigned it = 0; it < its; it++) {
      tmp[it].resize(N);
      std::memcpy(tmp[it].data(), values.data(), sizeof(unsigned)*N);
    }

    auto start = std::chrono::high_resolution_clock::now();
    for (unsigned it = 0; it < its; it++) {
      stlSort<unsigned>(tmp[it].data(), N, lessThan);
    }
    auto stop = std::chrono::high_resolution_clock::now();
    fprintf(stderr, "%5zu std::sort %10lldns\n", N , std::chrono::duration_cast<std::chrono::nanoseconds>((stop - start)).count()/its);
    N = 7 * N + N / 3;
  }



	
}