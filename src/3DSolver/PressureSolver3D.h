#ifndef PRESSURE_SOLVER_3D_H
#define PRESSURE_SOLVER_3D_H

#include "../LASolver/SparseItObj.h"
#include "Pool3D.h"

using namespace SparseItObj;

/**
 * Class which will implement the second order method for the pressure equation using Justin's code.
 * 
 * Very much a work in progress. Should interact with Pool well.
*/
class PressureSolver3D {
private:
    // TODO: decide which variables will go in here.
    int nx, ny, nz;
    double hx, hy, hz;
    int n;
    MatrixIter *matrix;
    double *x;
    double *y;
    double *z;
    double *pFlat;
    double *tol;

    // Helpers
    void setUpPressureMatrix(Pool3D *pool);
    void setUpMatrixStruct(int nx, int ny, int nz);
    void setUpPressureRHS(double dt, double ***FU, double ***FV, double ***FW, double ***p, Pool3D *pool);

    void discC(int row);
    void discNeumannBC(int row);

    // Discretization procedures depending on boundary information.
    void poolDisc(int row, int xInd, int yInd, int zInd, Pool3D *pool, MatrixIter *matrix);
public:
    PressureSolver3D() = default;
    PressureSolver3D(int nx, double *x, int ny, double *y, int nz, double *z);
    ~PressureSolver3D();
    void solvePCG(double dt, Pool3D *pool, double ***FU, double ***FV,
                  double ***FW, double ***p, bool matrixRecompute,
                  ParamIter &params);

};

#endif