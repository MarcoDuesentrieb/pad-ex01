#include <iostream>
#include <vector>

template <typename T, typename BinaryOpT, typename AccumulatorT>
auto matVec_generic(const std::vector<T>& A,
                    const std::vector<T>& x,
                    BinaryOpT binaryOp,        // T binaryOp(T lhs, T rhs)
                    AccumulatorT accumulator,  // T accumulator(T sum, T part)
                    T seed = {}) {
  if (A.size() % x.size() != 0)
    throw std::invalid_argument("wrong matrix size");

  // START WRITING YOUR CODE AFTER THIS LINE


// ----------------------------------
// ----------------------------------
// Happy Working wishes the ASC Team.
// ----------------------------------
// ----------------------------------


  // FINISH WRITING YOUR CODE BEFORE THIS LINE
}

int main() {
  std::vector<int> A(8 * 4, 0);
  std::vector<int> x = {
      0xDE, 0xAD, 0xBE, 0xEF, 0xBA, 0xAD, 0xF0, 0x0D,
  };
  static const int desiredCounts[] = {1, 3, 3, 7};

  // populate A
  for (std::size_t j = 0; j < sizeof(desiredCounts) / sizeof(int); ++j)
    for (int i = 0; i < desiredCounts[j]; ++i)
      A[j * 8 + i] = x[i] ^ 0x0F;

  auto binaryOp = [](int lhs, int rhs) { return lhs ^ rhs; };
  auto accumulator = [](int sum, int part) {
    return sum + (part == 0x0F ? 1 : 0);
  };
  // START WRITING YOUR CODE AFTER THIS LINE


// ----------------------------------
// ----------------------------------
// Happy Working wishes the ASC Team.
// ----------------------------------
// ----------------------------------


  // FINISH WRITING YOUR CODE BEFORE THIS LINE
  return 0;
}
