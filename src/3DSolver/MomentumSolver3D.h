#ifndef MOMENTUM_SOLVER_3D_H
#define MOMENTUM_SOLVER_3D_H

#include <iostream>
#include <functional>
#include <unordered_map>
#include "Pool3D.h"
#include "SolidObject3D.h"
#include "PressureSolver3D.h"
#include <string>
#include "../LASolver/SparseItObj.h"
#include "SimParams3D.h"

using namespace std;

/**
 * Abstract base class for 2D fluid solvers.
 * TODO: lots of design improvements needed if this is going to be generally usable.
*/
class MomentumSolver3D {
protected:
    // The Pool: the environment wherein the fluid is being solved. Can be updated
    //           at each point in time.
    Pool3D *pool;
    int nStructs;

    // Information related to the meshes.
    int nx, ny, nz;
    double dx, dy, dz;
    double gx, gy, gz;

    // Information related to the temporal domain.
    double t, dt, dtPrev;
    
    // Useful 1D meshes
    double *x;
    double *y;
    double *z;

    // Any model parameters.
    SimParams3D *params;

    // Boolean(s) indicating simulation state
    bool stepTaken;

    // The order of the numerical method. Influinces the number of ghost cells
    // required for the finite volume method.
    int methodOrd;

    // The quantities being solved for.
    double ***u;
    double ***v;
    double ***w;

    // Interpolated velocities to the cell centers
    double ***iu;
    double ***iv;
    double ***iw;

    double ***p;
    // // TODO: decide the size of these and if they are needed
    // double **px;
    // double **py;

    // The values for the explicit part of the method.
    double ***FU;
    double ***FV;
    double ***FW;
    // TODO: decide the size of these and if they are needed
    // double **FUx;
    // double **FVy;

    // Solver for the Pressure equation and its parameters
    PressureSolver3D *pressureSolver;
    ParamIter *pressureParams;


    // The barycentric weights used for flexible interpolation.
    double *baryWeights;

    // Use boundary condition callback function
    void (*applyFluidBCs)(int,int,int,double***,double***,double***);

    // Functions related to implementation of the explicit-implicit fluid model.
    void updateP();
    void updateU();
    void interpolateVelocities();
    void test_setInternal();
    virtual void updateF(Pool3D *pool) = 0;
    virtual void applyInterfaceBCs() = 0;
    virtual double getDt() = 0;
public:
    MomentumSolver3D(Boundary3D &boundary,
                     vector<SolidObject3D> &solidObjects,
                     SimParams3D &params,
                     std::function<void (int,int,int,int,double*,double*,double*,double***,double***,double***)> initialConditions,
                     void (*boundaryConditions)(int,int,int,double***,double***,double***));
    MomentumSolver3D() = default;
    ~MomentumSolver3D();
    double step(double tEnd, double safetyFactor);
    void writeToFile(const char *fname);
    void writePoolToFile(const char *poolFName, const char *poolVelFName);
    virtual double*** getU() = 0;
    virtual double*** getV() = 0;
    virtual double*** getW() = 0;
};

#endif