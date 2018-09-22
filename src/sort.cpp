#include <cstdio>
#include <chrono>
#include <vector>
#include <cassert>
#include <algorithm>

template<typename T>
__declspec(noinline) void insertionSort(T* data, size_t N, bool(*lessThan)(const T&, const T&))
{
  for (size_t j = 1; j < N; j++) {
    auto t(std::move(data[j]));
    size_t i = j;
    for (; 0 < i && lessThan(data[i - 1], t); i--) {
      data[i] = std::move(data[i - 1]);
    }
    data[i] = std::move(t);
  }
}


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

  std::vector<unsigned> values;
  std::vector<unsigned> gold;
  std::vector<std::vector<unsigned>> tmp;


  void runner(size_t N, const std::string& name, void(*impl)(unsigned* data, size_t N, bool(*)(const unsigned&, const unsigned&)))
  {
    unsigned its = 10;
    tmp.resize(10);
    for (unsigned it = 0; it < its; it++) {
      tmp[it].resize(N);
      std::memcpy(tmp[it].data(), values.data(), sizeof(unsigned)*N);
    }

    auto start = std::chrono::high_resolution_clock::now();
    for (unsigned it = 0; it < its; it++) {
      stlSort<unsigned>(tmp[it].data(), N, lessThan);
    }
    auto stop = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < N; i++) {
      assert(gold[i] == tmp.back()[i]);
    }
    fprintf(stderr, "%7zu %s %15lldns\n", N, name.c_str(), std::chrono::duration_cast<std::chrono::nanoseconds>((stop - start)).count() / its);
  }

}

int main(int argc, char** argv)
{
  std::srand(42);

  size_t N = 1;
  for (unsigned i = 0; i < 9; i++) {

    auto a = values.size();
    values.resize(N);
    for (auto k = a; k < N; k++) {
      values[k] = std::rand();
    }

    gold.resize(N);
    std::memcpy(gold.data(), values.data(), sizeof(unsigned)*N);
    std::sort(gold.begin(), gold.end());

    runner(N, "std::sort    ", stlSort);
    runner(N, "insertionSort", insertionSort);

    N = 7 * N + N / 3;
  }



	
}