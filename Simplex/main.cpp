#include <vector>
#include <cstdio>
#include <cassert>
#include <algorithm>

class LPProblem
{
public:
  enum struct TerminationType {
    OptimalSolution = 1,
    NoFeasibleSolution = 0,
    TooManyIterations = -1
  };

  LPProblem() = delete;

  LPProblem(std::initializer_list<std::initializer_list<double>> problem);

  TerminationType maximize(double* solution, const int maxIterations);

  void dump(const char* what);

private:

  void pivot(int j, int i);

  bool isSolutionFeasible();

  int basicVarIndex(int variable);

  int nonbasicVarIndex(int variable);

  void emitSolution(double* solution);

  int findPivotColumn();

  int findPivotRow(int pivotColumn);

  int findPivotColumnBland();

  int findPivotRowBland(int pivotColumn);

  TerminationType search(int maxIterations);

  int nonbasicVars;
  int basicVars;
  int stride;

  std::vector<int> basicVar;
  std::vector<int> nonbasicVar;
  std::vector<double> dictionary;
  std::vector<double> pivotColumnTmp;
  std::vector<double> objectiveTmp;       // nonbasicVars
};

LPProblem::LPProblem(std::initializer_list<std::initializer_list<double>> problem)
{
  assert(problem.size() > 1);
  assert(problem.begin()->size()  > 0);

  nonbasicVars = problem.begin()->size();
  basicVars = problem.size() - 1;
  stride = nonbasicVars + 2;

  dictionary.resize((basicVars+1)*stride);
  pivotColumnTmp.resize(basicVars + 1);
  objectiveTmp.resize(nonbasicVars);

  nonbasicVar.resize(nonbasicVars);
  for (int i = 0; i < nonbasicVars; i++) {
    nonbasicVar[i] = i+1;
  }

  basicVar.resize(basicVars);
  for (int i = 0; i < basicVars; i++) {
    basicVar[i] = nonbasicVars + i + 1;
  }

  int j = 0;
  for (auto & eq : problem) {
    assert(j <= basicVars);
    int i = 0;
    for (auto & c : eq) {
      assert(i <= nonbasicVars);
      if (j == 0) {
        dictionary[basicVars*stride + i + 1] = c;
      }
      else if(i<nonbasicVars) {
        dictionary[(j - 1)*stride + i + 1] = -c;
      }
      else {
        dictionary[(j - 1)*stride + 0] = c;
      }
      i++;
    }
    j++;
  }

}

void LPProblem::dump(const char* what)
{
  printf("%s:\n", what);

  std::vector<std::pair<int, int>> colPerm(nonbasicVars + 1);
  colPerm[0] = std::make_pair(0, 0);
  for (int i = 0; i < nonbasicVars; i++) {
    colPerm[i+1] = std::make_pair(nonbasicVar[i], i+1);
  }
  std::sort(colPerm.begin() + 1, colPerm.end(), [](auto & a, auto & b) -> bool {return a.first < b.first; });

  std::vector<std::pair<int, int>> rowPerm(basicVars);
  for (int i = 0; i < basicVars; i++) {
    rowPerm[i] = std::make_pair(basicVar[i], i);
  }
  std::sort(rowPerm.begin(), rowPerm.end(), [](auto & a, auto & b) -> bool {return a.first < b.first; });
  rowPerm.push_back(std::make_pair(0, basicVars));

  printf("           ");
  for (int i = 0; i < nonbasicVars; i++) {
    int l = colPerm[i + 1].second - 1;
    printf("  x%d  ", nonbasicVar[l]);
  }
  printf("\n");
  for (int j = 0; j <= basicVars; j++) {
    int k = rowPerm[j].second;
    if (j < basicVars) {
      printf("x%d = ", basicVar[k]);
    }
    else {
      printf(" z = ");
    }
    for (int i = 0; i < nonbasicVars + 1; i++) {
      int l = colPerm[i].second;
      printf("%5.2f ", dictionary[k*stride + l]);
    }
    printf("\n");
  }
}


bool LPProblem::isSolutionFeasible()
{
  for (int k = 0; k < basicVars; k++) {
    if (dictionary[k*stride] < 0.0) {
      return false;
    }
  }
  return true;
}


void LPProblem::pivot(int j, int i)
{
  // The pivot column will be a mix between the column
  // of the entering and leaving variable.

  // Store pivot column, set all entries except the one on pivot row to zero
  // such that these entries will accumulate the contents of the leaving
  // variable column. The pivot row entry will remain that of the entering
  // column such that its contents will be spread out into the leaving column's
  // entries
  for (int k = 0; k <= basicVars; k++) {
    pivotColumnTmp[k] = dictionary[k*stride + i];
    if (k != j) {
      dictionary[k*stride + i] = 0.0;
    }
  }

  // Conceptually make entries on the pivot column to zero by subtracting an
  // multiple of the pivot row. 
  for (int k = 0; k <= basicVars; k++) {
    if (k == j) continue;
    for (int l = 0; l < nonbasicVars + 1; l++) {
      if (l == i) continue;
      dictionary[k*stride + l] -= (pivotColumnTmp[k] * dictionary[j*stride + l]) / pivotColumnTmp[j];
    }
    // Existing value on the pivot column is in principle minus one (which it
    // isn't since we retain that value of the entering variable in the pivot column.
    dictionary[k*stride + i] += pivotColumnTmp[k] / pivotColumnTmp[j];
  }
  for (int l = 0; l < nonbasicVars + 1; l++) {
    if(l==i) continue;
    // Update pivot row. This is just scaled such that the pivot element becomes -1.
    dictionary[j*stride + l] = -dictionary[j*stride + l] / pivotColumnTmp[j];
  }
  // And finally, set the pivot element to the value of the leaving variable.
  dictionary[j*stride + i] = 1.0 / pivotColumnTmp[j];

  // And update book-keeping arrays.
  std::swap(basicVar[j], nonbasicVar[i-1]);
}

int LPProblem::basicVarIndex(int variable)
{
  for (int k = 0; k < basicVars; k++) {
    if (basicVar[k] == variable) {
      return k;
    }
  }
  return -1;
}


int LPProblem::nonbasicVarIndex(int variable)
{
  for (int l = 0; l < nonbasicVars; l++) {
    if (nonbasicVar[l] == variable) {
      return l;
    }
  }
  return -1;
}

int LPProblem::findPivotColumn()
{
  // Currently greedy. We might have change to Bland's rule.

  int column = -1;
  double best_val = std::numeric_limits<double>::epsilon();
  for (int k = 0; k < nonbasicVars; k++) {
    auto val = dictionary[basicVars*stride + k + 1];
    if (best_val < val) {
      column = k+1;
      best_val = val;
    }
  }

  // Returns -1 if solution is optimal.
  return column;
}


int LPProblem::findPivotRow(int pivotColumn)
{
  int pivotRow = -1;
  double worstRatio = std::numeric_limits<double>::max();

  for (int k = 0; k < basicVars; k++) {
    const auto r = dictionary[k*stride + pivotColumn];

    auto basicValue = dictionary[k*stride + 0];
    if (std::abs(basicValue) < std::numeric_limits<double>::epsilon()) {
      // We consider this a degeneracy.
      return -1;
    }

    if (r < -std::numeric_limits<double>::epsilon()) {
      const auto ratio = -dictionary[k*stride + 0] / r;
      if (ratio < worstRatio) {
        pivotRow = k;
        worstRatio = ratio;
      }
    }
  }

  return pivotRow;
}

int LPProblem::findPivotColumnBland()
{
  int column = -1;
  int lowest = std::numeric_limits<int>::max();
  for (int l = 0; l < nonbasicVars; l++) {
    auto i = nonbasicVar[l];
    auto v = dictionary[basicVars*stride + l + 1];
    if ((0.0 < v)&&(i < lowest)) {
      lowest = i;
      column = l + 1;
    }
  }
  return column;
}

int LPProblem::findPivotRowBland(int pivotColumn)
{
  int row = -1;
  int lowest = std::numeric_limits<int>::max();
  for (int k = 0; k < basicVars; k++) {
    auto j = basicVar[k];
    auto v = dictionary[k*stride + pivotColumn];
    if ((v < 0.0) && (j < lowest)) {
      lowest = j;
      row = k;
    }
  }
  return row;
}


void LPProblem::emitSolution(double* solution)
{
  for (int i = 0; i < basicVars; i++) {
    auto k = basicVar[i] - 1;
    if (k < nonbasicVars) {
      solution[k] = dictionary[i*stride + 0];
    }
  }
  for (int i = 0; i < nonbasicVars; i++) {
    auto k = nonbasicVar[i] - 1;
    if (k < nonbasicVars) {
      solution[k] = 0.0;
    }
  }
}

LPProblem::TerminationType LPProblem::search(int maxIterations)
{
  dump("Started search for optimum");
  for (int it = 0; it < maxIterations; it++) {
    int pivotColumn = findPivotColumn();
    if (pivotColumn < 0) {
      return TerminationType::OptimalSolution;
    }

    int pivotRow = findPivotRow(pivotColumn);
    if (pivotRow < 0) {
      // We're stalling, use Blands rule instead
      pivotColumn = findPivotColumnBland();
      assert(pivotColumn != -1);
      pivotRow = findPivotRowBland(pivotColumn);
      assert(pivotRow != -1);
    }

    pivot(pivotRow, pivotColumn);
    dump("Applied pivot");
  }
  return TerminationType::TooManyIterations;
}


LPProblem::TerminationType LPProblem::maximize(double* solution, const int maxIterations)
{
  if (isSolutionFeasible()) {
    dump("Initial solution is feasible");
  }
  else {
    dump("Initial solution is not feasible, must search for initial solution");
    
    // Store objective function
    for (int l = 0; l < nonbasicVars; l++) {
      objectiveTmp[l] = dictionary[basicVars*stride + l + 1];
    }

    // Extend dictionary with nonbasic variable x0 and add new objective function
    nonbasicVars++;
    nonbasicVar.push_back(0);
    for (int k = 0; k < basicVars; k++) {
      dictionary[k*stride + nonbasicVars] = 1.0;
    }
    dictionary[basicVars*stride + nonbasicVars] = -1.0;
    for (int l = 0; l < nonbasicVars; l++) {
      dictionary[basicVars*stride + l] = 0.0;
    }
 
    dump("Extended problem to find feasible solution");

    // Find the 'most infeasible' variable that should leave the basis for x0
    int leaving = -1;
    double mostInfeasibleValue = std::numeric_limits<double>::max();
    for (int k = 0; k < basicVars; k++) {
      auto value = dictionary[k*stride + 0];
      if (value < mostInfeasibleValue) {
        leaving = k;
        mostInfeasibleValue = value;
      }
    }
    // And swap x0 into the basis
    pivot(leaving, nonbasicVars);
    dump("Applied initial pivot, solution should be feasible");

    auto termination = search(maxIterations);
    if (termination != TerminationType::OptimalSolution) {
      return termination;
    }

    if (!isSolutionFeasible()) {
      return TerminationType::NoFeasibleSolution;
    }

    int x0_index = nonbasicVarIndex(0);
    if (x0_index < 0) {
      // x0 is in basis while we reached an optimum,
      // we conclude that problem is infeasible.
      return TerminationType::NoFeasibleSolution;
    }
      
    // Remove x0 from dictionary
    nonbasicVars--;
    for (int l = x0_index; l < nonbasicVars; l++) {
      nonbasicVar[l] = nonbasicVar[l + 1];
    }
    nonbasicVar.pop_back();
    for (int k = 0; k < basicVars; k++) {
      for (int l = x0_index; l < nonbasicVars; l++) {
        dictionary[k*stride + l + 1] = dictionary[k*stride + l + 2];
      }
    }

    // form original objective function in new dictionary
    dictionary[basicVars*stride + 0] = 0.0;
    for (int l = 0; l < nonbasicVars; l++) {
      auto i = nonbasicVar[l] - 1;
      dictionary[basicVars*stride + l + 1] = i < nonbasicVars ? objectiveTmp[i] : 0.0;
    }

    for (int k = 0; k < basicVars; k++) {
      auto j = basicVar[k] - 1;
      if (j < nonbasicVars) {
        auto s = objectiveTmp[j];
        for (int l = 0; l <= nonbasicVars; l++) {
          dictionary[basicVars*stride + l] += s*dictionary[k*stride + l];
        }
      }
    }
  }

  auto termination = search(maxIterations);
  if (termination == TerminationType::OptimalSolution) {
    emitSolution(solution);
    printf("At optimal solution:  %5.2f", solution[0]);
    for (int i = 1; i < nonbasicVars; i++) {
      printf(", %5.2f", solution[i]);
    }
    printf("\n");
  }
  return termination;
}



int main(int argc, char* argv[])
{
  std::vector<double> solution(100);
  if (true) {
    // Initial solution feasible
    // Solution is x1=32/29, x2=8/29, x3=30/29.
    LPProblem t({
      { 5,5,3 },
      { 1,3,1,3 },
      { -1,0,3,2 },
      { 2,-1,2,4 },
      { 2,3,-1,2 }
    });
    auto rv = t.maximize(solution.data(), 160);
    assert(rv == LPProblem::TerminationType::OptimalSolution);
    assert(std::abs(solution[0] - (32.0 / 29.0)) <= std::numeric_limits<double>::epsilon());
    assert(std::abs(solution[1] - (8.0 / 29.0)) <= std::numeric_limits<double>::epsilon());
    assert(std::abs(solution[2] - (30.0 / 29.0)) <= std::numeric_limits<double>::epsilon());
  }
  if (true) {
    // Initial solution not feasible.
    LPProblem t({
      { 1, -1,  1 },
      { 2, -1,  2,  4 },
      { 2, -3,  1, -5 },
      { -1,  1, -2, -1 }
    });
    auto rv = t.maximize(solution.data(), 160);
    assert(rv == LPProblem::TerminationType::OptimalSolution);
    assert(std::abs(solution[0] - 0.0) <= std::numeric_limits<double>::epsilon());
    assert(std::abs(solution[1] - (14.0 / 5.0)) <= 2.0*std::numeric_limits<double>::epsilon());
    assert(std::abs(solution[2] - (17.0 / 5.0)) <= std::numeric_limits<double>::epsilon());
  }
  if (true) {
    // Problem which cycles
    LPProblem t({
      { 10, -57, -9, -24 },
      { 0.5, -5.5, -2.5, 9, 0 },
      { 0.5, -1.5, -0.5, 1, 0 },
      { 1.0, 0.0, 0.0, 0, 1 }
    });

    auto rv = t.maximize(solution.data(), 160);
    assert(rv == LPProblem::TerminationType::OptimalSolution);
    assert(std::abs(solution[0] - 1.0) <= std::numeric_limits<double>::epsilon());
    assert(std::abs(solution[1] - 0.0) <= std::numeric_limits<double>::epsilon());
    assert(std::abs(solution[2] - 1.0) <= std::numeric_limits<double>::epsilon());
    assert(std::abs(solution[3] - 0.0) <= std::numeric_limits<double>::epsilon());
  }

  return 0;
}