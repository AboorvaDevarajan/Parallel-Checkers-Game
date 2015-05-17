#include "Test.h"

Test::Test() {
  value = 0;
  for (int i = 0; i < 8; ++i) {
    for (int j = 0; j < 8; ++j) {
      checkerBoard[i][j] = 0;
    }
  }
}
Test::~Test() {}
