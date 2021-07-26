#ifndef MASS_SPRING_3D_H
#define MASS_SPRING_3D_H

#include "SolidObject3D.h"
#include "Diff3DComparator.h"
#include <vector>
#include <map>
#include <set>
#include <queue>
#include "../LASolver/SparseItObj.h"
#include <Eigen/Dense>
#include "../Utils/Assembly.h"
#include "../Utils/ADMMPG.h"

class Pool3D;

using namespace SparseItObj;
using namespace std;

namespace mass_spring {

/**
 * Struct defining a 3D mass-point in the mass-spring system
*/
struct massPoint3D {
    double mass; // the mass of this point partical
    double x; // Location of this point in space
    double y;
    double z;
    double u;
    double v;
    double w;
    double sigU;
    double sigV;
    double sigW;
    vector<int> edgeIds; // Offsets into the edge list of this points neighbours
    vector<int> faceIds;
    bool boundaryPnt=false; // Indicate if this point is on the boundary
    int structNum;
    int nodeId;
};

/**
 * Struct defining a 3D edge in the mass-spring system
*/
struct edge3D {
    double l0; // the resting length of the spring
    int pntIds[2];
    bool boundaryEdge=false; // Indicate if this edge connects two boundary points
    bool structureEdge=false;
};

/**
 * Face defining a triangle on the surface of the MSS.
*/
struct face3D {
    int pntIds[3];
};


class MassSpring3D : Assembly {
public:
    /* Data for MSS stucture */
    vector<massPoint3D> *pntList;
    vector<edge3D> *edgeList;
    vector<int> *boundaryEdgeIdList;
    vector<int> *boundaryNodeIdList;
    vector<face3D> *faceList;

    /* Data for collision detection / handling */
    vector<set<massPoint3D*>> *nodeCols;
    double collisionStiffness; 
    double collisionDist;
    double repulseDist;

    // Algorithmic variables
    Eigen::VectorXd *q;
    Eigen::VectorXd *qt;
    Eigen::VectorXd *qprev;
    Eigen::VectorXd *f;

    // Backups for the state variables in case we wish to revert the state
    Eigen::VectorXd *qBackup;
    Eigen::VectorXd *fBackup;
    Eigen::VectorXd *qtBackup;
    Eigen::VectorXd *qprevBackup;
    void applyBodyForces();

    double hx, hy, hz;
    double gx, gy, gz;
    double pntMass;

    double eta;
    double E;
    double w;
    double dt, dtFix;
    int iterCount;

    double mass;
    double density;

    // Variables related to linear implicit method
    MatrixIter *matrix;
    int nnz;
    ParamIter *params;
    double *pcgRHS;
    double *pcgTol;
    double admmTol = 1e-6;

    int updateMode;
    int elementMode;
    bool initMode;

    SolidObject3D::ObjectType objType;

    // A queue to keep track of 
    // priority_queue<simstructs::surfaceDiff, vector<simstructs::surfaceDiff>, Diff3DComparator> heap;

    // Bounding box parameters
    int bb_jl = 0, bb_jr = 0, bb_il = 0, bb_ir = 0, bb_kl = 0, bb_kr = 0;

    // Constructors
    MassSpring3D(Pool3D &pool, int structNum, SolidObject3D &obj, int updateMode, int elementMode);
    MassSpring3D(const MassSpring3D &oldMSS);

    void updateSolidVels(double dt, Pool3D &pool, double ****stress, double fNet[3],
                    int ng, bool initMode);
    void updateSolidLocs(Pool3D &pool, bool interp);
    void interpBoundary(Pool3D &pool, bool resetRestingLengths);
    void outputNodes(const char* fname);
    void outputEdges(const char* fname);
    void outputNodeVels(const char *fname);
    void outputSurfaceCentroids(const char *fname);
    void interpFaceVels(double inPnt[3], double out[3]);
    void setAdmmTol(double tol);
    double closestBoundaryDist(double inPnt[3]);
    void closestBoundaryPnt(double inPnt[3], double outPnt[3]);
    ~MassSpring3D();
protected:
    /* Helper functions */
    void applyBoundaryForces(Pool3D &pool, double ****stress, int ng, double fNet[3], double colNet[3]);
    void calcLocalElasticForce(edge3D edge, int pntId1, massPoint3D pnt1, int pntId2, massPoint3D pnt2);
    void calcElasticForce(double E, double l0, massPoint3D pnt1,
                            massPoint3D pnt2, double force[6]);
    void calcLocalElasticHessian(double dt, edge3D edge, int pntId1, massPoint3D pnt1, int pntId2, massPoint3D pnt2);
    void calcLocalKelvinForce(edge3D edge, int pntId1, massPoint3D pnt1, int pntId2, massPoint3D pnt2);
    void computeCollisionStress(int nodeId, double colStress[3], double dA);
    void interpolateBoundaryLocation(Pool3D &pool, int i, int j, int k, double X[3]);
    bool findNearestGridNormalInterface(Pool3D &pool, int pntId, int nearestInterfacePnt[3]);
    double distToEdge(double x[3], int pntId1, int pntId2, double &t);

    /* Projective Dynamics functions */
    ADMMPG *admmSolver;
    void prox(double dt, Eigen::VectorXd &x, Eigen::VectorXd &DXpU, Eigen::VectorXd &z) override;
    void updateAfterStep(double dt, Eigen::VectorXd &xPrev, Eigen::VectorXd &x) override;
    void copyX(Eigen::VectorXd &tar) override;
    void predictX(double dt, Eigen::VectorXd &xPrev, Eigen::VectorXd &x, Eigen::VectorXd &xBar) override;
    void buildMassMatrix(Eigen::VectorXd &m) override;
    void buildDMatrix() override;
    void buildWMatrix(double w) override;


    ParamIter* setUpParams();
    MatrixIter* createStiffnessMatrix();

    // Different time discretization methods
    void eulerSolve(double dt, int elementMode, bool initMode);
    void linearImplicitSolve(double dt, int elementMode, bool initMode);
    void verletSolve(double dt, int elementMode, bool initMode);
    bool redundantEdges();
    int getNumNeighs(mass_spring::massPoint3D pnt);

    double pointDiff(massPoint3D pnt1, massPoint3D pnt2);
    double projTriangleDist(double x[3], int pntId1, int pntId2, int pntId3, double baryCoords[3]);
    void projTriangle(double x[3], int pntId1, int pntId2, int pntId3, double baryCoords[3]);
    void createFaceList();
    tuple<int, int, int> createOrderedTuple(int a, int b, int c);
    void computeCentroid(int pntId1, int pntId2, int pntId3, double centroid[3]);
    bool pntAlongEdge(double pnt[3], double X0[3], double X1[3]);
    double distToClosestBoundaryPoint(double x[3]);
    double distToClosestBoundaryEdge(double x[3]);
};

}


#endif