#ifndef NS_SOLVER_H
#define NS_SOLVER_H

#include "Boundary.h"
#include "Pool2D.h"
#include "PressureSolver.h"
#include "SimParams.h"
#include "SolidObject.h"
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

class NSSolver {
private:
    // The Pool: the environment wherein the fluid is being solved. Can be updated
    //           at each point in time.
    Pool2D *pool;
    int nStructs;

    // Information related to the meshes.
    int nx, ny;
    double dx, dy;

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
    // double **px;
    // double **py;

    // The values for the explicit part of the method.
    double **FU;
    double **FV;
    // TODO: decide the size of these and if they are needed
    // double **FUx;
    // double **FVy;

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
    // virtual void updateF(Pool2D *pool) = 0;
    // virtual void applyInterfaceBCs() = 0;
    // virtual double getDt() = 0;
    void test_setInternal();

    void updateF(Pool2D *pool);
    void applyInterfaceBCs();
    double getDt();
    double eno_usqx(int i, int j);
    double eno_vsqy(int i, int j);
    double eno_uvx(int i, int j);
    double eno_uvy(int i, int j);
    void interpolateCellCenterVelocities(int i, int j, double outU[2]);
public:
    NSSolver(Boundary &boundary,
                vector<SolidObject> &solidObjects,
                SimParams &params,
                std::function<void (int,int,int,double*,double*,double**,double**)> initialConditions,
                // void (*initialConditions)(int,int,int,double*,double*,double**,double**),
                void (*boundaryConditions)(int,int,double**,double**));
    int eno_get_upwind(bool sten[7], double ySten[7]);
    void (*initialConditions)(int,int,int,double*,double*,double**,double**);
    void (*boundaryConditions)(int,int,double**,double**);
    double getKinematicBoundaryLC(int i, int j, double velObj,
         double velNeigh, double cPnt[2], double stagPnt[2], int nDir[2]);
    NSSolver() = default;
    ~NSSolver();
    double step(double tEnd, double safetyFactor);
    void writeToFile(const char *fname);
    double** getU();
    double** getV();
    void writePoolToFile(const char *poolFName, const char *poolVelFName);

    void outputStructure(int structNum, const char *fname);
    void outputStructureNodes(int structNum, const char *fname);
    void outputStructureVels(int structNum, const char *fname);
    void outputTracers(const char *fname);
    void outputAllStructures(const char *fname);
    void outputAllStructureNodes(const char *fname);
    void outputAllStructureVels(const char *fname);
    void outputMedialAxis(const char *fname);

    double** testGetU();
    double** testGetV();
};

#endif
