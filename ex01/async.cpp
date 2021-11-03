#include <benchmark/benchmark.h>  // google benchmark
#include <algorithm>
#include <future>
#include <iostream>
#include <random>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

using ElementType = int;
using Ex1VectorType = std::vector<ElementType>;

//-------------------------------------------
// From ex1.1: For comparison with base line.
//-------------------------------------------
Ex1VectorType matVec(const Ex1VectorType& A, const Ex1VectorType& x) {
  if (A.size() % x.size() != 0)
    throw std::invalid_argument("wrong matrix size");

  std::size_t ni = x.size();
  std::size_t nj = A.size() / ni;
  Ex1VectorType result(nj, 0);  // initializes result to 0
  for (std::size_t j = 0; j < nj; ++j)
    for (std::size_t i = 0; i < ni; ++i)
      result[j] += A[ni * j + i] * x[i];

  // result is an rvalue as per C++14, §12.8, thus the move is implicit;
  // however, most compilers use NRVO here (i.e. the result object is
  // constructed in-place, i.e. neither copy nor move is necessary). C++17 even
  // mandates that the copy be elided.
  return result;
}

Ex1VectorType generateVectorOfRandomNumbers(std::size_t size, unsigned seed) {
  std::mt19937 gen(seed);
  std::uniform_int_distribution<std::int32_t> dis;

  Ex1VectorType result(size);
  for (auto& entry : result)
    entry = dis(gen);

  return result;
}

struct TestDataSet {
  Ex1VectorType A;
  Ex1VectorType x;

  TestDataSet(std::size_t nX, std::size_t nY, unsigned seed)
      : A(generateVectorOfRandomNumbers(nX * nY, seed)),
        x(generateVectorOfRandomNumbers(nX, seed)) {}
};
//--------------------------------------
// Ex1.3: Actual implementation of Ex1.3
//--------------------------------------

template <typename ResultIterT,
          typename MatIterT,
          typename VecIterT,
          typename BinaryOpT,
          typename AccumulatorT>
void matVec_generic_iterator(
    ResultIterT itResult,
    MatIterT firstA,
    MatIterT lastA,
    VecIterT firstX,
    VecIterT lastX,
    BinaryOpT binaryOp,        // T binaryOp(T lhs, T rhs)
    AccumulatorT accumulator,  // T accumulator(T sum, T part)
    typename std::iterator_traits<ResultIterT>::value_type seed = {}) {
  static_assert(
      std::is_same<
          typename std::iterator_traits<ResultIterT>::value_type,
          typename std::iterator_traits<MatIterT>::value_type>::value &&
          std::is_same<
              typename std::iterator_traits<ResultIterT>::value_type,
              typename std::iterator_traits<VecIterT>::value_type>::value,
      "all iterators should have the same value_type");

  // We can't easily check for dimension errors here because we don't know
  // whether itResult points to a large enough data structure, and because
  // MatIterT and VecIterT might be iterators without random access.
  // START WRITING YOUR CODE AFTER THIS LINE

  std::size_t ni = x.size();
  std::size_t nj = A.size() / ni;
  Ex1VectorType result(nj, 0);  // initializes result to 0
  for (std::size_t j = 0; j < nj; ++j)
    for (std::size_t i = 0; i < ni; ++i)
      result[j] += A[ni * j + i] * x[i];
            
  std::size_t nA = lastA - firstA;
  std::size_t nX = lastX - firstX;
            
  for (std::vector<int>::iterator itX = firstX; itX != lastX; ++itX) {
    for(std::vector<int>::iterator itA = firstA; itA != lastA; ++itA) {
      
    }
  }

// ----------------------------------
// ----------------------------------
// Happy Working wishes the ASC Team.
// ----------------------------------
// ----------------------------------


  // FINISH WRITING YOUR CODE BEFORE THIS LINE
}

template <typename ResultIterT,
          typename MatIterT,
          typename VecIterT,
          typename BinaryOpT,
          typename AccumulatorT>
void matVec_parallel(
    ResultIterT itResult,
    MatIterT firstA,
    MatIterT lastA,
    VecIterT firstX,
    VecIterT lastX,
    BinaryOpT binaryOp,        // T binaryOp(T lhs, T rhs)
    AccumulatorT accumulator,  // T accumulator(T sum, T part)
    unsigned parallelism = 1,
    typename std::iterator_traits<ResultIterT>::value_type seed = {}) {
  static_assert(
      std::is_same<
          typename std::iterator_traits<ResultIterT>::value_type,
          typename std::iterator_traits<MatIterT>::value_type>::value &&
          std::is_same<
              typename std::iterator_traits<ResultIterT>::value_type,
              typename std::iterator_traits<VecIterT>::value_type>::value,
      "all iterators should have the same value_type");

  // here we have to assume that MatIterT and VecIterT are random access
  // iterators anyway, so we can do at least some error checking
  std::size_t nA = lastA - firstA;
  std::size_t nX = lastX - firstX;
  if (nA % nX != 0)
    throw std::invalid_argument("wrong matrix size");
  std::size_t nY = nA / nX;

  // find out how many tasks we need
  unsigned numThreads = std::min(parallelism, static_cast<unsigned>(nY));
  std::size_t delta = (nY + numThreads - 1) / numThreads;

  for (std::size_t j = 0; j < nj; ++j)
    for (std::size_t i = 0; i < ni; ++i)
      result[j] += A[ni * j + i] * x[i];
            
  // set up tasks
  // START WRITING YOUR CODE AFTER THIS LINE
  std::future<void> fut;
  
  for(int i = 0; i < numThreads; ++i) {
    fut = std::async(std::launch::async, matVec_generic_iterator, itResult, firstA+i*, lastA, firstX, lastX, binaryOp, accumulator, seed);
  }
  fut.wait()


// ----------------------------------
// ----------------------------------
// Happy Working wishes the ASC Team.
// ----------------------------------
// ----------------------------------


  // FINISH WRITING YOUR CODE BEFORE THIS LINE
}

//--------------------------
// Google Benchmark Template
//--------------------------
static void Ex1Arguments(benchmark::internal::Benchmark* b) {
  unsigned int nthreads = std::thread::hardware_concurrency();
  const auto lowerLimit = 1;
  const auto upperLimit = 100000;
  // Generate x,y
  for (auto x = lowerLimit; x <= upperLimit; x *= 10) {
    for (auto y = lowerLimit; y <= upperLimit; y *= 10) {
      for (auto t = 1; t <= nthreads; t = t << 1) {
        b->Args({x, y, t});
      }
    }
  }
}

void setCustomCounter(benchmark::State& state) {
  state.counters["x"] = state.range(0);
  state.counters["y"] = state.range(1);
  state.counters["threads"] = state.range(2);
  auto BytesProcessed = int64_t(state.iterations()) * int64_t(state.range(0)) *
                        int64_t(state.range(1)) * sizeof(ElementType);
  state.counters["bytes"] = BytesProcessed;
  state.SetBytesProcessed(BytesProcessed);
  auto FLOP = int64_t(state.iterations()) * int64_t(state.range(0)) *
              int64_t(state.range(1)) * 2;
  state.counters["FLOP"] = FLOP;
  state.counters["FLOPperIter"] = FLOP / state.iterations();
  state.counters["sizeof(A+x+y)"] =
      sizeof(ElementType) * (int64_t(state.range(0)) * int64_t(state.range(1)) +
                             int64_t(state.range(0)) + int64_t(state.range(1)));
}

// Multithreaded
static void benchAsync(benchmark::State& state) {
  auto seed = 42;
  TestDataSet data{static_cast<std::size_t>(state.range(0)), static_cast<std::size_t>(state.range(1)), static_cast<unsigned>(seed)};
  auto y = matVec(data.A, data.x);

  // storage for result of parallel computation
  Ex1VectorType yp(data.A.size() / data.x.size());

  Ex1VectorType result(state.range(1));
  for (auto _ : state) {
    matVec_parallel(yp.begin(), data.A.begin(), data.A.end(), data.x.begin(),
                    data.x.end(), std::multiplies<ElementType>(),
                    std::plus<ElementType>(), state.range(2));
    benchmark::DoNotOptimize(yp.data());
    benchmark::ClobberMemory();
  }
  setCustomCounter(state);
  if (y != yp)
    throw std::runtime_error("wrong result.");
}
BENCHMARK(benchAsync)->Apply(Ex1Arguments)->UseRealTime();

BENCHMARK_MAIN();
