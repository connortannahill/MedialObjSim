#ifndef NS_SOLVER_3D_H
#define NS_SOLVER_3D_H

#include "MomentumSolver3D.h"
#include "Boundary3D.h"
#include "SolidObject3D.h"
#include "Pool3D.h"
#include <unordered_map>
#include <string>
#include "SimParams3D.h"

using namespace std;

class NSSolver3D : MomentumSolver3D {
private:
    void updateF(Pool3D *pool);
    void resetBoundaryConditions();
    void applyInterfaceBCs();
    double getDt();

    // Helpers
    double getKinematicBoundaryLC(int i, int j, int k, double velObj, double velNeigh,
                                    double cPnt[3], double stagPnt[3], int nDir[3]);
public:
    NSSolver3D(Boundary3D &boundary,
            vector<SolidObject3D> &solidObjects,
            SimParams3D &params,
            void (*initialConditions)(int,int,int,int,double*,double*,double*,double***,double***,double***),
            void (*boundaryConditions)(int,int,int,double***));
    NSSolver3D() = default;
    double step(double tEnd, double safetyFactor);
    void writeToFile(const char *fname);
    void writePoolToFile(const char *poolName, const char *poolVelName);
    double*** getU();
    double*** getV();
    double*** getW();

    void outputStructure(int structNum, const char *fname);
    void outputStructureNodes(int structNum, const char *fname);
    void outputStructureVels(int structNum, const char *fname);
    void outputAllStructures(const char *fname);
    void outputAllStructureNodes(const char *fname);
    void outputAllStructureVels(const char *fname);


    double*** testGetU();
    double*** testGetV();
    double*** testGetW();
};

#endif