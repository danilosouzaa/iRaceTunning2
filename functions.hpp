#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <cfloat>
#include <iostream>
//#include <time.h>
#include <limits.h>
#include <list>
#include <sys/time.h>
#include <omp.h>
#include <vector>
#include <map>
#include <set>
#include <string.h>

using namespace std;

typedef struct{
    int numberVariables;
    int numberConstraints;
    int cont;
    double *Coefficients;
    int *Elements;
    int *ElementsConstraints;
    double *rightSide;
    double *xAsterisc;
    int *typeVariables;
    double *ub;
    double *lb;
}constraintsReal;


int contCasasDecimais(double num);

int cutMaxDivisorCommonRec(int m, int n);

int cutMaxDivisorCommonVector(std::vector<double> coefs, int nElem);

int heuristicConstraintsCoefFrac(std::vector<double> &temp, int szConstraints);

double valueViolation(constraintsReal *cCover, constraintsReal *constraintsSmall, int idCover, int ogConstraint, int precision);

int verifyRepeatedIncidency(int **matrizIncidencia, constraintsReal *originalConstraints, int posCover);

constraintsReal *createCutsCover(constraintsReal *cutsCover, constraintsReal *constraintsOriginal, constraintsReal *constraintsSmall, int constraint, int nCuts, int precision);

void quicksortCof(double *values, int *idc, int began, int end);

constraintsReal *AllocStrConstraintsReal(int cont, int nConstraints, int nVariables);

void showStructFull(constraintsReal *constraintsFull, std::vector<std::string> nameConstraints, std::vector<std::string> nameVariables);

constraintsReal *prepareConstraintsOfLifting(constraintsReal *constraintsInitial, int precision);

void freeStrConstraintsReal(constraintsReal *cut);

constraintsReal *removeNegativeCoefficientsAndSort(constraintsReal *constraintsOriginal, int *convertVector);

void SortByCoefficients(constraintsReal *h_cut);

void quicksortDouble(double *values, int began, int end);

int *LCIBallas(int *coverSolution, constraintsReal *constraintsCover, int constraint);

int *LCIAdam(int *coverSolution, constraintsReal *constraintsCover, int constraint);

int verifyViolationGreedy(int *solutionFinal, constraintsReal *constraitsUsed, int constraint, int precision);

double avaliaSolution(constraintsReal *knapsackConstraints, int *solution, int constraint, int precision, int typeLift);

int *greedyInitialSolutionSA(constraintsReal *knapsackConstraints, int constraint, int precision, int typeLift, int typeSolutionInitial);

void copyAndVerifyPoolSolution(int *solution, int sz, int *poolSolution, int *numberSolutionAtual, double violation, double *violationFinal);

double fRand(double fMin, double fMax);

void shuffleVectorInt(int *vec, int sz);

int *neighborhoodOne(constraintsReal *knapsackConstraints, int precision, int constraint, int *solutionCurrent);

int *neighbohoodTwo(constraintsReal *knapsackConstraints, int precision, int constraint, int *solutionCurrent);

int *neighbohoodThree(constraintsReal *knapsackConstraints, int precision, int constraint, int *solutionCurrent);

int *neighbohoodFour(constraintsReal *knapsackConstraints, int precision, int constraint, int *solutionCurrent);

int *neighbohoodFive(constraintsReal *knapsackConstraints, int precision, int constraint, int *solutionCurrent);

constraintsReal *copyStrConstraintsReal(constraintsReal *originalRows);

constraintsReal *runCCwithSA(constraintsReal *binaryRows, int precision, int szPoolCutsMax, float TempInitial, float FactorResf, int typeLift, int minimal, double v1, double v2, double v3, double v4, double *violationFinal, int typeSolutionInitial, double timeLeft, int fm, double timeConstraints);

constraintsReal *returnVariablesOriginals(constraintsReal *constraintsOriginal, int *convertVector, int precision, int nVariablesInitial);

constraintsReal *convertBinaryOfOriginalConstraints(constraintsReal *constraintsOriginal, constraintsReal *constraintsBinary, int nInitialBinary);

std::vector<std::string> renamedNameConstraints(std::vector<std::string> nameConstraints, int typeContraints, int szNewConstraints, int szCuts, int lastCut);

//constraintsReal *LiftWithSA(constraintsReal *constraintsInitial /*,std::vector<std::string> nameConstraints, std::vector<std::string> nameVariables*/);
constraintsReal *LiftWithSA(constraintsReal *constraintsInitial, 
                                double tempInitial, double fatorResf, double pv1,double pv2, double pv3, double pv4, double pv5,
                                double timeConstraints, int fm, double timeLeft, int typeLift, int typeSolutionInitial,
                                int precision, int minimal, int szPoolCutsMaxCC);

constraintsReal *useMatrixConflict(constraintsReal *constraintsInitial, std::vector<std::vector<int>> matrixConflict, std::map<int, int> &varXclique);

constraintsReal *returnOriginalConflict(constraintsReal *cOriginal, constraintsReal *cNew, std::vector<std::vector<int>> matrixConflict, std::map<int, int> varXclique, int nCuts);

constraintsReal *readPseudoInstance(const char* nameFile, std::map< int, std::set<int> > &conflict);

bool conflicting(std::map<int, std::set<int>> conflictMapSet, int el1, int el2);

int *vectorNonRepeteadNonDominated(constraintsReal *constraintsOriginal, int nConstraintsInitial);

int verifyRepeatCuts(constraintsReal *constraintsOriginal, int cutOriginal, int cutCreate);

int verifyDominanceCG(constraintsReal *constraintsInitial, int c1, int c2);

double countViolation(constraintsReal *constraintsInitial, int *constraintsValited, int nConstraints);
