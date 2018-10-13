#include <vector>
#include <cstdio>
#include <cassert>

// d,e,f are slack variables
//
// d =  5 -2a -3b  -c
// e = 11 -4a  +b +2c
// f =  8 -3a -4b -2c
// z =     5a +4b +3c
//

// find initial feasible solution: set a,b,c=0 -> d=5, e=11, f=8, z=0,

//              0   0   0
// d     =  5 -2a -3b  -c  =  5
//   e   = 11 -4a  +b +2c  = 11
//     f =  8 -3a -4b -2c  =  8
// z     =     5a +4b +3c  =  0

// try to increase the value of a:  d =  5 -2a -3*0  -0 > 0
//                                  a = 5/2   
//
// set this as new solution: a=5/2, b=0, c=0, d=0, e=1, f=1/2, z=25/2:
//
// move a to right-hand side:  d =  5 -2a -3b -c <=> a =  5/2 -d/2 -(3/2)b  -c/2
// and substitute expression for a into the other
//
// a     = (5/2) (-1/2)d  (-3/2)b  (-1/2)c
//   e   = 11 -4a  +b +2c  = 11
//     f =  8 -3a -4b -2c  =  8
// z     =     5a +4b +3c  =  0

// Tableau form
// a | b | c | d | e | f | |
// --+---+---+---+---+---++---
// 2 | 3 | 1 | 1 | 0 | 0 |=|  5
// 4 | 1 | 2 | 0 | 1 | 0 |=| 11
// 3 | 4 | 2 | 0 | 0 | 1 |=|  8
// --+---+---+---+---+---+-+---
// 5 | 4 | 3 | 0 | 0 | 0 |=|  0
//
// 1) Examine coefficients in last row, if all are neg or zero, tableau is optimal.
//    Find the largest. This is the pivot column. Above is the first (a)
// 2) For each row whose entry in the pivot column is positive, find the one with
//    smallest s/r, s is from the rightmost column. This is the pivot row. Above is the first.
// 3) Pivot number is in the intersection of pivot row and and pivot column, here is 2.
//    Divide every entry in pivot row by pivot number:
// a |  b  |  c  |  d  | e | f | |
// --+-----+-----+-----+---+---++---
// 1 | 3/2 | 1/2 | 1/2 | 0 | 0 |=|  5
// 4 |  1  |  2  |  0  | 1 | 0 |=| 11
// 3 |  4  |  2  |  0  | 0 | 1 |=|  8
// --+-----+-----+-----+---+---+-+---
// 5 |  4  |  3  |  0  | 0 | 0 |=|  0
// 4) For every other row, subtract a suitable number of pivot row s.t. entries in pivot row become zero.

struct Tableau
{
  Tableau() = delete;

  Tableau(std::initializer_list<std::initializer_list<double>> problem);

  double& operator()(int eq, int var);

  bool step();

  void dump();

  int decisionVariables;
  int equations;
  int rows;
  int cols;
  std::vector<int> variables;
  std::vector<double> data;
};

Tableau::Tableau(std::initializer_list<std::initializer_list<double>> problem)
{
  assert(problem.size() > 1);
  assert(problem.begin()->size()  > 0);
  equations = problem.size() - 1;
  decisionVariables = problem.begin()->size() - 1;

  rows = equations + 1;
  cols = decisionVariables + equations + 1;
  data.resize(rows*cols);
  variables.resize(cols-1);
  for (int i = 0; i < cols - 1; i++) {
    variables[i] = i;
  }

  int j = 0;
  for (auto & eq : problem) {
    int i = 0;
    for (auto & c : eq) {
      if (i < decisionVariables) {
        operator()(j, i) = c;
      }
      else {
        assert(i == decisionVariables);
        operator()(j, cols - 1) = c;
      }
      i++;
    }
    j++;
  }


  for (int j = 0; j < equations; j++) {
    operator()(j, decisionVariables+j) = 1.0;
  }
}

void Tableau::dump()
{
  for (int j = 0; j < rows; j++) {
    if (j + 1 == rows) {
      for (int i = 0; i < cols; i++) {
        if (i + 1 == cols) {
          printf("+-----");
        }
        else {
          printf("------");
        }
      }
      printf("\n");
    }
    for (int i = 0; i < cols; i++) {
      if (i + 1 == cols) {
        printf("|");
      }
      printf("%5.2f ", operator()(j, i));
    }
    printf("\n");
  }
}



double& Tableau::operator()(int eq, int var)
{
  return data[eq*cols + var];
}

bool Tableau::step()
{
  // Find pivot column
  int best_col_ix = -1;
  double best_col_val = std::numeric_limits<double>::epsilon();
  for (int k = 0; k < decisionVariables; k++) {
    int ix = variables[k];
    if (best_col_val < data[equations*cols + ix]) {
      best_col_ix = ix;
      best_col_val = data[equations*cols + ix];
    }
  }

  if (best_col_ix < 0) {
    printf("At optimal solution.");
    return false;
  }

  // Find pivot row
  int best_row_ix = -1;
  double smallest_s_over_r = std::numeric_limits<double>::max();
  for (int k = 0; k < equations; k++) {
    const auto r = data[equations*cols + best_col_ix];
    if (std::numeric_limits<double>::epsilon() < r) {
      const auto s_over_r = data[equations*cols + cols-1]/r;
      if (s_over_r < smallest_s_over_r) {
        best_row_ix = k;
        smallest_s_over_r = s_over_r;
      }
    }
  }

  if (best_row_ix < 0) {
    printf("Stalling.\n");
    return false;
  }

//  double pivot_val = data[equations]  


  printf("col ix=%d (%f), row ix=%d (%f).\n", best_col_ix, best_col_val, best_row_ix, smallest_s_over_r);

  return true;
}

#if 0
int main(int argc, char* argv[])
{
  Tableau t( {
    {2, 3, 1, 5},
    {4, 1, 2, 11},
    {3, 4, 2, 8},
    {5, 4, 3, 0}
  });

  t.dump();
  printf("---\n");
  t.step();
  printf("---\n");
  t.dump();

  return 0;
}
#endif
