#ifndef NS_SOLVER_H
#define NS_SOLVER_H

#include "MomentumSolver2D.h"
#include "Boundary.h"
#include "SolidObject.h"
#include "Pool2D.h"
#include <unordered_map>
#include "SimParams.h"
#include <string>
#include <vector>

using namespace std;

class NSSolver : MomentumSolver2D {
private:
    void updateF(Pool2D *pool);
    void applyFluidBCs();
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
                void (*initialConditions)(int,int,int,double*,double*,double**,double**));
    int eno_get_upwind(bool sten[7], double ySten[7]);
    double getKinematicBoundaryLC(int i, int j, double velObj,
         double velNeigh, double cPnt[2], double stagPnt[2], int nDir[2]);
    NSSolver() = default;
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