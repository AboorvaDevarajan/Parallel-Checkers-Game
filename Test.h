#ifndef CASE_H_
#define CASE_H_
#include <sstream>
#include <string>

using namespace std;

class Test {
public:
  int mFrom;
  int mTo;
  int bFrom;
  int bTo;
  int checkerBoard[8][8];
  int value;
  Test();
  virtual ~Test();
};

#endif
