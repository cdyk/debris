#include <cstdio>
#include <chrono>
#include <vector>
#include <cassert>
#include <algorithm>

namespace {
  
  const size_t quickSortInsertionSortThreshold = 32;

  template<typename T>
  void insertionSort(T* data, size_t N, bool(*lessThan)(const T&, const T&))
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
  T& medianOfThree(T& a, T& b, T& c, bool(*lessThan)(const T&, const T&))
  {
    auto & min_ab = lessThan(a, b) ? a : b;
    auto & max_ab = lessThan(a, b) ? b : a;
    auto & min_max_ab_c = lessThan(max_ab, c) ? max_ab : c;
    auto & pivot = lessThan(min_ab, min_max_ab_c) ? min_max_ab_c : min_ab;

    if (a != b && b != c && a != c) {
      //assert(a < pivot || b < pivot || c < pivot);
      //assert(pivot < a || pivot < b || pivot < c);
    }
    else {
      //assert(a <= pivot || b <= pivot || c <= pivot);
      //assert(pivot <= a || pivot <= b || pivot <= c);
    }

    return pivot;
  }

  template<typename T>
  size_t partition(T* data, size_t N, bool(*lessThan)(const T&, const T&))
  {
    auto & pivot = data[N - 1];
    size_t i0 = 0;
    size_t i1 = N - 1;
    while (true) {

      // invariant: data[i<i0] < pivot
      // i0 < N should be unnecessary since pivot is at end and will terminate
      for (; lessThan(data[i0], pivot); i0++) {
        //assert(i0 < N && "Pivot is last element and can't be less than itself");
      }

      // invariant: data[i1 <= i] >= pivot
      for (; 0 < i1 && !lessThan(data[i1 - 1], pivot); i1--);

      //assert(i0 <= i1 && "i0 and i1 should not pass due to predicate");
      if (i0 == i1) break;
      else {
        std::swap(data[i0], data[i1 - 1]);
      }
    }
    return i0;
  }


  template<typename T>
  __declspec(noinline) void quickSort(T* data, size_t N, bool(*lessThan)(const T&, const T&))
  {
    auto & pivot = medianOfThree(data[0], data[N / 2], data[N - 1], lessThan);
    std::swap(pivot, data[N - 1]);

    auto m = partition(data, N, lessThan);
    std::swap(data[m], data[N - 1]);

    if (m) quickSort(data, m, lessThan);
    if (2 < N - m) quickSort(data + m + 1, N - m - 1, lessThan);
  }


  template<typename T>
  void quickSortInsertion(T* data, size_t N, bool(*lessThan)(const T&, const T&))
  {
    if (N < quickSortInsertionSortThreshold) {
      for (size_t j = 1; j < N; j++) {
        auto t(std::move(data[j]));
        size_t i = j;
        for (; 0 < i && !lessThan(data[i - 1], t); i--) {
          data[i] = std::move(data[i - 1]);
        }
        data[i] = std::move(t);
      }
    }
    else {
      auto & pivot = medianOfThree(data[0], data[N / 2], data[N - 1], lessThan);
      std::swap(pivot, data[N - 1]);

      auto m = partition(data, N, lessThan);
      std::swap(data[m], data[N - 1]);

      if (m) quickSortInsertion(data, m, lessThan);
      if (2 < N - m) quickSortInsertion(data + m + 1, N - m - 1, lessThan);
    }
  }


}


template<typename T>
void insertionSortRunner(T* data, size_t N, bool(*lessThan)(const T&, const T&))
{
  insertionSort(data, N, lessThan);
}

template<typename T>
__declspec(noinline) void quickSortRunner(T* data, size_t N, bool(*lessThan)(const T&, const T&))
{
  if (N <= 1) return;
  quickSort(data, N, lessThan);
}

template<typename T>
__declspec(noinline) void quickSortFatRunner(T* data, size_t N, bool(*lessThan)(const T&, const T&))
{
  if (N <= 1) return;
  quickSortFat(data, N, lessThan);
}

template<typename T>
__declspec(noinline) void quickSortInsertionRunner(T* data, size_t N, bool(*lessThan)(const T&, const T&))
{
  quickSortInsertion(data, N, lessThan);
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

    values.resize(N);
    gold.resize(N);
    for (unsigned l = 0; l < 2; l++) {
      if (l == 0) {
        for (auto k = 0; k < N; k++) {
          values[k] = std::rand();
        }
      }
      else if (l == 1) {
        for (auto k = 0; k < N; k++) {
          values[k] = values[k] & 1;
        }
      }

      std::memcpy(gold.data(), values.data(), sizeof(unsigned)*N);
      std::sort(gold.begin(), gold.end());

      fprintf(stderr, "l=%d ------------------------------------------------------------\n", l);
      runner(N, "std::sort         ", stlSort<unsigned>);
      if (i < 5)runner(N, "insertionSort     ", insertionSortRunner<unsigned>);
      runner(N, "quickSort         ", quickSortRunner<unsigned>);
      runner(N, "quickSortInsertion", quickSortInsertionRunner<unsigned>);
    }

    N = 7 * N + N / 3;
  }



	
}