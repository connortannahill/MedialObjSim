#ifndef PRESSURE_SOLVER_H
#define PRESSURE_SOLVER_H

#include "SparseItObj.h"
#include "Pool2D.h"

using namespace SparseItObj;

/**
 * Class which will implement the second order method for the pressure equation using Justin's code.
 * 
 * Very much a work in progress. Should interact with Pool well.
*/
class PressureSolver {
private:
    // TODO: decide which variables will go in here.
    int nx;
    int ny;
    double hx;
    double hy;
    int n;
    MatrixIter *matrix;
    double *x;
    double *y;
    double *pFlat;
    double *tol;

    double **reconMatrix;

    // Helpers
    void setUpPressureMatrix(Pool2D &pool);
    void setUpPressureRHS(double dt, double **FU, double **FV, double **p, Pool2D &pool);

    // Discretization procedures depending on boundary information.
    void discSW(int xInd, int yInd, int row, MatrixIter *matrix);
    void discS(int xInd, int yInd, int row, MatrixIter *matrix);
    void discSE(int xInd, int yInd, int row, MatrixIter *matrix);

    void discW(int xInd, int yInd, int row, MatrixIter *matrix);
    void discC(int xInd, int yInd, int row, MatrixIter *matrix);
    void discE(int xInd, int yInd, int row, MatrixIter *matrix);

    void discNW(int xInd, int yInd, int row, MatrixIter *matrix);
    void discN(int xInd, int yInd, int row, MatrixIter *matrix);
    void discNE(int xInd, int yInd, int row, MatrixIter *matrix);
public:
    PressureSolver() = default;
    PressureSolver(int nx, double *x, int ny, double *y);
    ~PressureSolver();
    // void solveSOR(int niter, double omega, double dt, int nx, double *x, int ny, double *y,
    //                  double **FU, double **FV, double **p);
    void solvePCG(double dt, Pool2D &pool, double **FU, double **FV, double **p, bool iluRecompute, ParamIter &params);

};

#endif