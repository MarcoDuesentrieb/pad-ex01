#include <benchmark/benchmark.h>  // google benchmark
#include <random>                 // random number generator
#include <vector>                 // vector
#include <iostream>
#include <numeric>
#include <stdexcept>

using ElementType = int;
using Ex1VectorType = std::vector<ElementType>;

//-----------------------------------------------
// Generate random number vector and TestDataSets
//-----------------------------------------------

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

//-------------------------------
// Actual assignment: ex01 part 1
//-------------------------------
Ex1VectorType matVec(const Ex1VectorType& A,
                                          const Ex1VectorType& x) {
  if (A.size() % x.size() != 0)
    throw std::invalid_argument("wrong matrix size");

// START WRITING YOUR CODE AFTER THIS LINE


    size_t m = A.size();
    size_t n = x.size();
    size_t o = m/n; // size of the product p = M * x
    
    // the solution is a m/n x 1 matrix
    Ex1VectorType p(o);
  
    for(int i=0; i<o; i++)
    {
        for(int j=0; j<n; j++)
            p.at(i) += A.at(i*n+j) * x.at(j);
    }
    return p;
    


// FINISH WRITING YOUR CODE BEFORE THIS LINE
}

void matVec(Ex1VectorType& y, const Ex1VectorType& A, const Ex1VectorType& x) {
  if (A.size() % x.size() != 0)
    throw std::invalid_argument("wrong matrix size");

// START WRITING YOUR CODE AFTER THIS LINE

    size_t m = A.size();
    size_t n = x.size();
    size_t o = m/n; // size of the product y = M * x

    for(int i=0; i<o; i++)
    {
        for(int j=0; j<n; j++)
            y.at(i) += A.at(i*n+j) * x.at(j);
    }

// FINISH WRITING YOUR CODE BEFORE THIS LINE
}

//--------------------------
// Google Benchmark Template
//--------------------------
static void Ex1Arguments(benchmark::internal::Benchmark* b) {
  const auto lowerLimit = 1;
  const auto upperLimit = 100000;
  // Generate x,y
  for (auto x = lowerLimit; x <= upperLimit; x *= 10) {
    for (auto y = lowerLimit; y <= upperLimit; y *= 10) {
      b->Args({x, y});
    }
  }
}

void setCustomCounter(benchmark::State& state) {
  state.counters["x"] = state.range(0);
  state.counters["y"] = state.range(1);
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

// move explicitly
static void benchMove(benchmark::State& state) {
  auto seed = 42;
  TestDataSet data{std::size_t(state.range(0)), std::size_t(state.range(1)), (unsigned int)seed};
  Ex1VectorType result;
  for (auto _ : state) {
    result = matVec(data.A, data.x);
    benchmark::DoNotOptimize(result.data());
    benchmark::ClobberMemory();
  }
  setCustomCounter(state);
}
BENCHMARK(benchMove)->Apply(Ex1Arguments);

// reference
static void benchReference(benchmark::State& state) {
  auto seed = 42;
  TestDataSet data{(std::size_t)state.range(0), (std::size_t)state.range(1), (unsigned int)seed};
  Ex1VectorType result(state.range(1));
  for (auto _ : state) {
    matVec(result, data.A, data.x);
    benchmark::DoNotOptimize(result.data());
    benchmark::ClobberMemory();
  }
  setCustomCounter(state);
}
BENCHMARK(benchReference)->Apply(Ex1Arguments);

// simple test if calc is correct
static void isCorr(benchmark::State& state) {
  auto seed = 42;
  TestDataSet data{std::size_t(state.range(0)), std::size_t(state.range(1)), (unsigned int)seed};
  std::size_t nj = data.x.size();
  std::size_t ni = data.A.size() / nj;
  Ex1VectorType check;
  Ex1VectorType result(ni,0);
  for (std::size_t j = 0; j < nj; ++j)
    for (std::size_t i = 0; i < ni; ++i)
      result[i] += data.A[nj * i + j] * data.x[j];
  check = matVec(data.A,data.x); 
  for(std::size_t i = 0; i < ni; ++i)
	  if(check[i] != result[i])
		  throw std::runtime_error("Wrong result!");
}

BENCHMARK(isCorr)->Apply(Ex1Arguments);

BENCHMARK_MAIN();
