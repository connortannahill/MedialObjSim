#ifndef MOMENTUM_SOLVER_H
#define MOMENTUM_SOLVER_H

#include <iostream>
#include <functional>
#include <unordered_map>
#include "Pool2D.h"
#include "SolidObject.h"
#include "PressureSolver.h"
#include <string>
#include "../LASolver/SparseItObj.h"
#include "SimParams.h"

using namespace std;

/**
 * Abstract base class for 2D fluid solvers.
 * TODO: lots of design improvements needed if this is going to be generally usable.
*/
class MomentumSolver2D {
protected:
    // The Pool: the environment wherein the fluid is being solved. Can be updated
    //           at each point in time.
    Pool2D *pool;
    int nStructs;

    // Information related to the meshes.
    int nx, ny;
    double dx, dy;
    double gx, gy;

    // Information related to the temporal domain.
    double t, dt, dtPrev;
    
    // Useful 1D meshes
    double *x;
    double *y;

    // Hashmap containing any model parameters.
    SimParams *params;

    // Boolean(s) indicating simulation state
    bool stepTaken;

    // The order of the numerical method. Influinces the number of ghost cells
    // required for the finite volume method.
    int methodOrd;

    // The quantities being solved for.
    double **u;
    double **v;

    // Interpolated velocities to the cell centers
    double **iu;
    double **iv;

    double **p;
    // TODO: decide the size of these and if they are needed
    double **px;
    double **py;

    // The values for the explicit part of the method.
    double **FU;
    double **FV;
    // TODO: decide the size of these and if they are needed
    double **FUx;
    double **FVy;

    // Solver for the Pressure equation and its parameters
    PressureSolver *pressureSolver;
    ParamIter *pressureParams;


    // The barycentric weights used for flexible interpolation.
    double *baryWeights;

    // Use boundary condition callback function
    void (*applyFluidBCs)(int,int,double**,double**);

    // Functions related to implementation of the explicit-implicit fluid model.
    void updateP();
    void updateU();
    void interpolateVelocities();
    virtual void updateF(Pool2D *pool) = 0;
    virtual void applyInterfaceBCs() = 0;
    virtual double getDt() = 0;
    void test_setInternal();
public:
    MomentumSolver2D(Boundary &boundary,
                     vector<SolidObject> &solidObjects,
                     SimParams &params,
                     void (*initialConditions)(int,int,int,double*,double*,double**,double**),
                     void (*boundaryConditions)(int,int,double**,double**));
    MomentumSolver2D() = default;
    ~MomentumSolver2D();
    double step(double tEnd, double safetyFactor);
    void writeToFile(const char *fname);
    virtual double** getU() = 0;
    virtual double** getV() = 0;
    void writePoolToFile(const char *poolFName, const char *poolVelFName);
    // virtual double** getP() = 0;
};

#endif