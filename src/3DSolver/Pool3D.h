#ifndef POOL_3D_H
#define POOL_3D_H

#include "Boundary3D.h"
#include "SolidObject3D.h"
#include "../Utils/SimUtilities.h"
#include <utility>
#include <map>
#include "../Utils/GridVal3D.h"
#include "../Utils/GridValComparator.h"
#include <queue>
#include "MassSpring3D.h"
#include <vector>
#include "SimParams3D.h"
#include <nanoflann.hpp>
#include "KDPointCloud3D.h"

using namespace mass_spring;
using namespace nanoflann;
using namespace std;

/**
 * Enumeration which labels each point according to type in the FSI
 * model.
 * Note interfaces are defined as follows on each level. With 
 * 
 *   (U|D| ))NW      (U|D| )) N      (U|D| ))NE
 *               |---------------|
 *               |               |
 *               |               |
 *   (U|D| ))W   |   (U|D| ))C   |   (U|D| ))E
 *               |               |
 *               |---------------|
 *   (U|D| ))SW      (U|D| )) S      (U|D| )) SE
*/
namespace objects {
    enum FSIObject {
        NORTH=2,
        EAST=3,
        SOUTH=5,
        WEST=7,
        UP=11,
        DOWN=13,

        // All 27(!!!) second cases
        // all cases * UP*DOWN*NORTH*EAST*SOUTH*WEST

        // With up direction
        FLUID_USW = DOWN * EAST*NORTH,
        FLUID_US  = DOWN * WEST*NORTH*EAST,
        FLUID_USE = DOWN * WEST*NORTH,
        FLUID_UW  = DOWN * EAST*SOUTH*NORTH,
        FLUID_UE  = DOWN * WEST*NORTH*SOUTH, 
        FLUID_UNW = DOWN * SOUTH*EAST, 
        FLUID_UN  = DOWN * EAST*SOUTH*WEST, 
        FLUID_UNE = DOWN * SOUTH*WEST, 

        // The 3D cases
        FLUID_SW = UP*DOWN * EAST*NORTH,
        FLUID_S  = UP*DOWN * WEST*NORTH*EAST,
        FLUID_SE = UP*DOWN * WEST*NORTH,
        FLUID_W  = UP*DOWN * EAST*SOUTH*NORTH,
        FLUID_E  = UP*DOWN * WEST*NORTH*SOUTH, 
        FLUID_NW = UP*DOWN * SOUTH*EAST, 
        FLUID_N  = UP*DOWN * EAST*SOUTH*WEST, 
        FLUID_NE = UP*DOWN * SOUTH*WEST, 

        // With the down direction
        FLUID_DSW = UP * EAST*NORTH,
        FLUID_DS  = UP * WEST*NORTH*EAST,
        FLUID_DSE = UP * WEST*NORTH,
        FLUID_DW  = UP * EAST*SOUTH*NORTH,
        FLUID_DE  = UP * WEST*NORTH*SOUTH, 
        FLUID_DNW = UP * SOUTH*EAST, 
        FLUID_DN  = UP * EAST*SOUTH*WEST, 
        FLUID_DNE = UP * SOUTH*WEST, 

        FLUID_U = DOWN*EAST*NORTH*WEST*SOUTH,
        FLUID_D = UP*EAST*NORTH*WEST*SOUTH,

        // The easy cases
        FLUID_C = -1,
        UNLABELLED_INTERFACE = -3,
        STRUCTURE = -2
    };
}

extern int DOMAIN_UNDESCOVERED_3D;
extern int DOMAIN_FLUID_3D;

/**
 * This class is the model for the fluid and all of the solids with it.
 * Holds onto the meshes, and keeps track of what is in each grid via an
 * enumeration.
 * 
 * TODO: going to have to refactor this whole thing at some point for better safety
 *       and less memory copying.
*/
class Pool3D {
protected:
    /* KDTree variables */
    typedef KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<double, KDPointCloud3D> ,
                                    KDPointCloud3D, 3> KDTree3D;
    KDPointCloud3D *kdPointCloud;
    KDTree3D *kdTree;
    void detectCollisions();

    vector<tuple<double, double, double>> *medialAxisPnts;

    // Number of mesh points in each direction
    int nx;
    int ny;
    int nz;

    double *x;
    double *y;
    double *z;

    bool isDeformable;

    objects::FSIObject ***pool;

    // Flag to indicate whether the enumeration has changed
    bool enumChanged;

    int methodOrd;

    // Embedding tolerance
    double eps;

    // The mass-spring systems that are used for representation of the deformable solids.
    vector<MassSpring3D> *solids;

    // Number of structs
    int nStructs;

    // Fluid viscosity
    double mu;

    // Externally generated velocity field for the pool
    double ***poolU;
    double ***poolV;
    double ***poolW;

    // Tracers
    simstructs::tracer3D *tracers;

    // Domain tracking tensor
    int ***domainTracker;

    // Fast marching tracker
    int ***fastMarchingState;

    // Heap for fast marching method
    priority_queue<GridVal3D, vector<GridVal3D>, GridValComparator> heap;
    
    // The level set equation.
    double ***phi;
    double ***phiReInit;
    double ***phiRk1;
    double ***phiRk2;

    // Function to create the 2D pool from the given objects
    void create3DPool(Boundary3D &boundary,
                        vector<SolidObject3D> &structures,
                        SimParams3D &params);
    void embedShape(SolidObject3D &structure, int structNum);
    void labelInterface(int i, int j, int k);
    void updatePoolVelocities(double dt, double ***u, double ***v, double ***w, double ***p, int ng);

    void fastMarch(bool nExtrap, int mode);
    void fastMarchSetNVal(int i, int j, int k, bool nExtrap, int mode);
    void assignDomainMemberships(int i, int j, int k, int mode);
    bool indInRange(int i, int j, int k);
    double buildSqeezeField();
    bool shouldRefitSDF(double tol);
    void refitToSolids(int ng);

    /* Functions for discretizations of HJ equations */

    // Optimal TVD third order TVD Runge-Kutta time-stepper for general HJ equation using ENO discretization
    // TODO: refactor code using a functor to bind the method to a particular instance instead of passing the pool
    friend void tvdRK3HJ(double dt, double ***arr, Pool3D *pool, int bcType,
                         double (Pool3D::*rhs)(int,int,int,double*,double*,double*,double*,double*,double*),
                         void (Pool3D::*bcFun)(double***));
    friend void rhsHJENO(double hx, double hy, double hz, double ***in,
              double ***out, int bcType, Pool3D *pool,
              double (Pool3D::*rhs)(int,int,int,double*,double*,double*,double*,double*,double*),
              void (Pool3D::*bcFun)(double***));

    // RHS functions using ENO discretizations. These all take in a collocation of m(x|y)Vals - the relative distances between
    // mesh points on a uniformly spaced grid in each space dimension, and (x|y)Vals, the values used to construct the interpolant in the
    // ENO schemes. The size of each array is arr[2*methodOrd+1], where methodOrd is the order of accuracy of the method.
    double levelSetRHS_ENO3(int i,int j, int k, double *mxVals, double *xVals,
                            double *myVals, double *yVals, double *mzVals, double *zVals);
    double velocityExtrapolationRHS_ENO3(int i, int j, int k, double *mxVals, double *xVals,
                                        double *myVals, double *yVals, double *mzVals, double *zVals);
    double signedDistanceReinitializationRHS_CorrectUpwind(int i, int j, int k, double *mxVals,
                                                        double *xVals, double *myVals, 
                                                        double *yVals, double *mzVals, double *zVals);

    // Boundary conditition functions
    void applyLevelSetPeriodicBCs(double ***phi);
    void applyReinitializationBCs(double ***phi);
    void applyReinitializationSlice(int k, double ***phi);

    // Helper functions
    void levelSetUnitNormal(int i, int j, int k, double n[3]);
    double smearSign(double h, double phi);
    double sign(double phi);
    double phiSmearSign(double h, double phi, double phi_x, double phi_y, double phi_z);
    bool isStructConnected(int i, int j, int k);
    bool isFluidConnected(int i, int j, int k);
    double distanceFromInterface(int i, int j, int k);
    void bfsFromTracer(int structNum);
    void updateTracer(int structNum, double dt, int mode);
    void calculateLocalFNet(int i, int j, int k, objects::FSIObject obj, int ng,
                            double ***u, double ***v, double ***w, double ***p);
    void computeBoundaryStress(int i, int j, int k, objects::FSIObject obj, int ng,
                                    double ***u, double ***v, double ***w, double ***p);
    void multiStructureVelocitySet(double ***u, double ***v, double ***w, int ng);

    // Misc. helpers
    void applyPeriodicLevelSetSlice(int zi);

public:
    double dtFix;

    double repulseDist;
    double collisionStiffness;
    double collisionDist;
    int repulseMode;
    
    bool oneGridFromInterface(int i, int j, int k);
    bool oneGridFromInterfaceStructure(int i, int j, int k);
    int domainMembership(int i, int j, int k);
    bool enumHasChanged();
    Pool3D(Boundary3D &boundary,
            vector<SolidObject3D> &structures, SimParams3D &params);
    Pool3D() = default;
    ~Pool3D();
    void levelSetGradient(int i, int j, int k, double g[3]);
    void updatePool(double dt, double ***u, double ***v, double ***w,
                    double ***p, int ng, bool reinitialize);
    void setUpDomainArray();
    void printPool();
    void outputPool(const char *fname);
    void outputPoolVelocity(const char *fname);
    void outputSurfaceCentroids(int structNum, const char *fname);
    void outputStructure(int structNum, const char *fname);
    void outputAllStructures(const char *baseName);
    void outputStructureNodes(int structNum, const char *fname);
    void outputAllStructureNodes(const char *baseName);
    void outputStructureVels(int structNum, const char *fname);
    void outputAllStructureVels(const char *baseName);
    objects::FSIObject objAtIndex(int xInd, int yInd, int zInd);
    bool hasStructInDir(objects::FSIObject obj, objects::FSIObject dir);
    bool isInterface(objects::FSIObject obj);
    bool isNormalInterface(objects::FSIObject obj);
    bool isUpdateableU(int i, int j, int k);
    bool isUpdateableV(int i, int j, int k);
    bool isUpdateableW(int i, int j, int k);
    void enumeratePool();
    void getNormalDir(objects::FSIObject obj, int nDir[3]);
    double interpolatePhi(double x, double y, double z);
    void interpolatePhiGrad(double x, double y, double z, double phiGrad[3]);
    double closestBoundaryDist(double inPnt[3]);


    // Getters and setters
    int getNx();
    int getNy();
    int getNz();
    double getXMeshVal(int i);
    double getYMeshVal(int j);
    double getZMeshVal(int k);
    double getPhiVal(int i, int j, int k);
    double getObjU(int i, int j, int k);
    double getObjV(int i, int j, int k);
    double getObjW(int i, int j, int k);

    double *getXMesh();
    double *getYMesh();
    double *getZMesh();
};

#endif