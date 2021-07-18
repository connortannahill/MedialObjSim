#ifndef POOL_2D_H
#define POOL_2D_H

#include "Boundary.h"
#include "SolidObject.h"
#include "../Utils/SimUtilities.h"
#include <utility>
#include <map>
#include "../Utils/GridVal2D.h"
#include "../Utils/GridValComparator.h"
#include <queue>
#include "MassSpring2D.h"
#include <vector>
#include "SimParams.h"
#include <nanoflann.hpp>
#include "./KDPointCloud.h"

using namespace mass_spring;
using namespace nanoflann;

/**
 * Enumeration which labels each point according to type in the FSI
 * model.
 * Note interfaces are defined as follows:
 * 
 *   NW      N      NE
 *   |---------------|
 *   |               |
 *   |               |
 * W |       C       | E
 *   |               |
 *   |---------------|
 *   SW      S       SE
*/
namespace objects {
    enum FSIObject {
        NORTH=2,
        EAST=3,
        SOUTH=5,
        WEST=7,
        FLUID_SW = EAST*NORTH,
        FLUID_S  = WEST*NORTH*EAST,
        FLUID_SE = WEST*NORTH,
        FLUID_W  = EAST*SOUTH*NORTH,
        FLUID_E  = WEST*NORTH*SOUTH, 
        FLUID_NW = SOUTH*EAST, 
        FLUID_N  = EAST*SOUTH*WEST, 
        FLUID_NE = SOUTH*WEST, 
        FLUID_C = -1,
        STRUCTURE = -2,
        UNLABELLED_INTERFACE = -3,
        UNLABELLED = -4
    };
}

extern int DOMAIN_UNDESCOVERED;
extern int DOMAIN_FLUID;
extern int DOMAIN_INTERSECTION;

/**
 * This class is the model for the fluid and all of the solids with it.
 * Holds onto the meshes, and keeps track of what is in each grid via an
 * enumeration.
 * 
 * TODO: going to have to refactor this whole thing at some point for better safety
 *       and less memory copying.
*/
class Pool2D {
protected:
    /* KDTree variables and collision detection */
    typedef KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<double, KDPointCloud> ,
                                    KDPointCloud, 2> KDTree;
    KDPointCloud *kdPointCloud;
    KDTree *kdTree;

    void detectCollisions();

    /* Essential parameters */
    int nx;
    int ny;
    int methodOrd;

    /* Variables for geometric algorithms */
    double *x;
    double *y;
    objects::FSIObject **pool;
    int **domainTracker;
    int **fastMarchingState;
    std::priority_queue<GridVal2D, std::vector<GridVal2D>, GridValComparator> heap;
    std::vector<std::tuple<double, double>> *medialAxisPnts;

    // The mass-spring systems that are used for representation of the deformable solids.
    std::vector<MassSpring2D> *solids;

    // Fluid viscosity
    double mu;

    int nSteps;

    // Number of structs
    int nStructs;

    // Externally generated velocity field for the pool
    double **poolU;
    double **poolV;
    
    // The level set equation.
    double **phi;
    double **phiReInit;
    double **phiRk1;
    double **phiRk2;

    // Indicator of whether we are using the mass-spring system
    bool isDeformable;

    // Tracer particles for each of the solid objects.
    simstructs::tracer2D *tracers;

    // Function to create the 2D pool from the given objects
    void create2DPool(Boundary &boundary,
                      std::vector<SolidObject> &structures,
                      SimParams &params);
    void embedShape(SolidObject &structure, int structNum);
    void labelInterface(int j, int i);
    void updateTracer(int structNum, double dt, int mode);
    void updatePoolVelocities(double dt, double **u, double **v, double **p, int ng);
    double buildSqeezeField();

    /* Functions for discretizations of HJ equations */
    // Optimal TVD third order TVD Runge-Kutta time-stepper for general HJ equation using ENO discretization
    // TODO: refactor code using a functor to bind the method to a particular instance instead of passing the pool
    friend void tvdRK3HJ(double dt, double **arr, Pool2D *pool, int bcType,
                         double (Pool2D::*rhs)(int,int,double*,double*,double*,double*,int),
                         void (Pool2D::*bcFun)(double**));
    friend void rhsHJENO(double hx, double hy, double **in, double **out, int bcType, Pool2D *pool,
              double (Pool2D::*rhs)(int,int,double*,double*,double*,double*,int),
              void (Pool2D::*bcFun)(double**));

    /* RHS functions using ENO discretizations. These all take in a collocation of m(x|y)Vals - the relative distances between
    * mesh points on a uniformly spaced grid in each space dimension, and (x|y)Vals, the values used to construct the interpolant in the
    * ENO schemes. The size of each array is arr[2*methodOrd+1], where methodOrd is the order of accuracy of the method. */
    double levelSetRHS_ENO3(int i,int j, double *mxVals, double *xVals,
                            double *myVals, double *yVals, int maxOrd);
    double signedDistanceReinitializationRHS_ENO3(int i,int j, double *mxVals,
                                                  double *xVals, double *myVals,
                                                  double *yVals);
    double signedDistanceReinitializationRHS_CorrectUpwind(int i,int j, double *mxVals,
                                                  double *xVals, double *myVals,
                                                  double *yVals);
    double velocityExtrapolationRHS_ENO3(int i,int j, double *mxVals, double *xVals,
                                         double *myVals, double *yVals);
    void fastMarchPool(bool nExtrap, int mode);
    void bfsSearchForEdges(int structNum);
    void fastMarchSetNVal(int i, int j, bool nExtrap, int mode);

    /* Boundary conditions */
    void applyLevelSetPeriodicBCs(double **phi);
    void applyReinitializationBCs(double **phi);
    void applyHomogenousNeumann(double **arr);

    /* Helper functions */
    double distanceFromInterface(int i, int j);
    void getWidestStencil(int minIndex, int maxIndex, int curIndex, bool *stencil);
    void levelSetUnitNormal(int i, int j, double n[2]);
    double smearSign(double h, double phi);
    double sign(double phi);
    double phiSmearSign(double h, double phi, double phi_x, double phi_y);
    bool isStructConnected(int i, int j);
    bool isFluidConnected(int i, int j);
    void assignDomainMemberships(int i, int j, double val, int mode);
    bool indInRange(int i, int j);
    bool shouldRefitSDF(double tol);


public:
    bool enoInRangeX(int val);
    bool enoInRangeY(int val);
    int getPossibleStencil(int i, int j, int axis, int methodOrd, bool *stencil);//, bool ghostSkip);
    // Collision params, TODO use proper OOP here.
    double repulseDist;
    double collisionStiffness;
    double collisionDist;
    double hx, hy;
    int repulseMode;
    // Get the closest interface point (approximately) to an input point
    Pool2D(Boundary &boundary, 
            std::vector<SolidObject> &structures,
            SimParams &params);
    Pool2D() = default;
    ~Pool2D();

    double gx, gy;

    /* Control method */
    void updatePool(double dt, double **u, double **v, double **p, int ng, bool reinitialize);

    /* Output methods */
    void printPool();
    void printTracers();
    void outputPool(const char *fname);
    void outputPoolVelocity(const char *fname);
    void outputStructure(int structNum, const char *fname);
    void outputAllStructures(const char *fname);
    void outputStructureNodes(int structNum, const char *fname);
    void outputAllStructureNodes(const char *fname);
    void outputStructureVels(int structNum, const char *fname);
    void outputAllStructureVels(const char *fname);
    void outputTracers(const char *fname);
    void outputMedialAxisApprox(const char *fname);

    /* Helper methods */
    void levelSetGradient(int i, int j, double g[2]);
    objects::FSIObject objAtIndex(int xInd, int yInd);
    bool hasStructInDir(objects::FSIObject obj, objects::FSIObject dir);
    bool isInterface(objects::FSIObject obj);
    bool isNormalInterface(objects::FSIObject obj);
    bool isUpdateableU(int i, int j);
    bool isUpdateableV(int i, int j);
    bool isUsableU(int i, int j);
    bool isUsableV(int i, int j);
    int buildDxStencil(int i, int j, bool sten[7]);
    int buildDyStencil(int i, int j, bool sten[7]);
    void enumeratePool();
    void resetLevelSet();
    void getClosestInterfacePoint(double P[2], double XC[2]);
    bool oneGridFromInterface(int i, int j);
    bool oneGridFromInterfaceStructure(int i, int j);
    void getNormalDir(objects::FSIObject obj, int nDir[2]);
    double interpolatePhi(double x, double y);
    void interpolatePhiGrad(double x, double y, double phiGrad[2]);
    double getLevelSetVal(int i, int j);
    void setUpDomainArray();
    void bfsFromTracer(int structNum);
    void printDomainMatrix();
    int domainMembership(int i, int j);

    /* Solid coupling methods */
    void calculateLocalFNet(int i, int j, objects::FSIObject obj, int ng,
                            double **u, double **v, double **p);
    void computeBoundaryStress(int i, int j, objects::FSIObject obj, int ng,
                                    double **u, double **v, double **p);
    void multiStructureVelocitySet(double **u, double **v, int ng);

    /* Public vars TODO encapsulate */
    double dtFix;

    /* Getters and setters */
    int getNx();
    int getNy();
    double getXMeshVal(int i);
    double getYMeshVal(int j);
    double getPhiVal(int i, int j);
    double getObjU(int i, int j);
    double getObjV(int i, int j);
    double *getXMesh();
    double *getYMesh();
    void refitToSolids(int ng);

};


#endif