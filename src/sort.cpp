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
    for (; 0 < i && !lessThan(data[i - 1], t); i--) {
      data[i] = std::move(data[i - 1]);
    }
    data[i] = std::move(t);
  }
}

template<typename T>
__declspec(noinline) void recursiveQuickSort(T* data, size_t N, bool(*lessThan)(const T&, const T&))
{
  if (N < 3) return;
  auto& a = data[0];
  auto& b = data[N / 2];
  auto& c = data[N - 1];



  auto & min_ab = lessThan(a, b) ? a : b;
  auto & max_ab = lessThan(a, b) ? b : a;
  auto & min_max_ab_c = lessThan(max_ab, c) ? max_ab : c;
  auto & pivot = lessThan(min_ab, min_max_ab_c) ? min_max_ab_c : min_ab;

  if (a != b && b != c && a != c) {
    assert(a < pivot || b < pivot || c < pivot);
    assert(pivot < a || pivot < b || pivot < c);
  }
  else {
    assert(a <= pivot || b <= pivot || c <= pivot);
    assert(pivot <= a || pivot <= b || pivot <= c);
  }
  std::swap(pivot, data[N - 1]);
  pivot = data[N - 1];

  size_t i0 = 0;
  size_t i1 = N - 1;
  while (i0 + 1 < i1) {
    for (; i0 < N && lessThan(data[i0], pivot); i0++);
    for (; i0 < i1 && lessThan(pivot, data[i1 - 1]); i1--);

    std::swap(data[i0], data[i1]);
    auto t(std::move(data[i0]));
    data[i0] = std::move(data[i1]);
    data[i1] = std::move(t);
  }
  recursiveQuickSort(data, i0, lessThan);
  recursiveQuickSort(data + i0, N - i0, lessThan);

  int g = 2;
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
    unsigned its = 1;// 10;
    tmp.resize(its);
    for (unsigned it = 0; it < its; it++) {
      tmp[it].resize(N);
      std::memcpy(tmp[it].data(), values.data(), sizeof(unsigned)*N);
    }

    auto start = std::chrono::high_resolution_clock::now();
    for (unsigned it = 0; it < its; it++) {
      impl(tmp[it].data(), N, lessThan);
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

    if (i != 0) {
      fprintf(stderr, "------------------------------------------------------------\n");
    }
    runner(N, "std::sort         ", stlSort<unsigned>);
    if (N < 5)runner(N, "insertionSort     ", insertionSort<unsigned>);
    runner(N, "recursiveQuickSort", recursiveQuickSort<unsigned>);

    N = 7 * N + N / 3;
  }



	
}