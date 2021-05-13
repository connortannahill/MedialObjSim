#ifndef MASS_SPRING_2D_H
#define MASS_SPRING_2D_H

#include "SolidObject.h"
#include <vector>
#include <set>
#include "../LASolver/SparseItObj.h"
#include "Eigen/Dense"
#include "../Utils/Assembly.h"
#include "../Utils/ADMMPG.h"

using namespace SparseItObj;
using namespace std;

class Pool2D;

namespace mass_spring {

/**
 * Struct defining a 2D mass-point in the mass-spring system
*/
struct massPoint2D {
    double mass; // the mass of this point partical
    double x; // Location of this point in space
    double y;
    double u;
    double v;
    double sigU;
    double sigV;
    vector<int> edgeIds; // Offsets into the edge list of this points neighbours
    bool boundaryPnt=false; // Indicate if this point is on the boundary
    int structNum;
    int nodeId;
};

/**
 * Struct defining a 2D edge in the mass-spring system
*/
struct edge2D {
    double l0; // the resting length of the spring
    int pntIds[2];
    bool boundaryEdge=false; // Indicate if this edge connects two boundary points
};


class MassSpring2D : Assembly {
public:
    /* Data for the MSS structure */
    vector<massPoint2D> *pntList;
    vector<edge2D> *edgeList;
    vector<int> *boundaryEdgeIdList;
    vector<int> *boundaryNodeIdList;
    SolidObject::ObjectType objType;

    /* Data for collision detection / handling */
    vector<set<massPoint2D*>> *nodeCols;
    double collisionStiffness; 
    double collisionDist;
    double repulseDist;
    double admmTol = 1e-6;

    /* Algorithmic variables and control modes */
    Eigen::VectorXd *q;
    double dt, dtFix;
    Eigen::VectorXd *qt;
    Eigen::VectorXd *qprev;
    Eigen::VectorXd *f;

    // Backups for the state variables in case we wish to revert the state
    Eigen::VectorXd *qBackup;
    Eigen::VectorXd *qtBackup;
    Eigen::VectorXd *qprevBackup;

    int elementMode;
    bool initMode;

    /* Semi-implicit method parameters */
    MatrixIter *matrix;
    int nnz;
    ParamIter *params;
    double *pcgRHS;
    double *pcgTol;

    // bool *edgeChecked;
    int updateMode;
    double pntMass;

    double eta;
    double E;
    int iterCount;

    double mass;
    double density;
    double w;

    // Bounding box parameters
    int bb_jl = 0, bb_jr = 0, bb_il = 0, bb_ir = 0;

    // Constructors
    MassSpring2D(Pool2D &pool, int structNum, SolidObject &obj, int updateMode, int elementMode);
    MassSpring2D(const MassSpring2D &oldMSS);

    /* Control methods */
    void updateSolidVels(double dt, Pool2D &pool, double ***stress, double fNet[2],
                    int ng, bool initMode);
    void updateSolidLocs(Pool2D &pool, bool interp);
    void setAdmmTol(double tol);

    /* Helper methods */
    void outputNodes(const char* fname);
    void outputEdges(const char* fname);
    void outputNodeVels(const char *fname);
    void interpBoundary(Pool2D &pool, bool resetRestingLengths);
    void interpFaceVels(double inPnt[2], double out[2]);
    double closestBoundaryDist(double inPnt[2]);
    double closestBoundaryPnt(double inPnt[2], double outPnt[2]);
    bool intersectsBoundary(double pnt1[2], double pnt2[2]);

    /* Destructor */
    ~MassSpring2D();
protected:
    /* Energy functions */
    void applyBoundaryForces(Pool2D &pool, double ***stress, int ng, double fNet[2]);
    void calcLocalElasticForce(edge2D edge, int pntId1, massPoint2D pnt1, int pntId2, massPoint2D pnt2);
    void calcElasticForce(double E, double l0, massPoint2D pnt1,
                            massPoint2D pnt2, double force[4]);
    void calcLocalElasticHessian(double dt, edge2D edge, int pntId1, massPoint2D pnt1, int pntId2, massPoint2D pnt2);
    void calcLocalKelvinForce(edge2D edge, int pntId1, massPoint2D pnt1, int pntId2, massPoint2D pnt2);
    void computeCollisionStress(int nodeId, double colStress[2]);

    /* Helper functions */
    void interpolateBoundaryLocation(Pool2D &pool, int i, int j, double X[2]);
    bool findNearestGridNormalInterface(Pool2D &pool, int pntId, int nearestInterfacePnt[2]);
    double projectOntoEdge(int pntId1, int pntId2, double inPnt[2], double &t);
    int getNumNeighs(mass_spring::massPoint2D pnt);
    void structToLoc(mass_spring::massPoint2D pnt, double loc[2]);
    void computeInterpolationStress(Pool2D &pool, double x[2], double stiffness, double out[2]);
    void reset();

    /* Projective Dynamics functions */
    ADMMPG *admmSolver;
    void prox(double dt, Eigen::VectorXd &x, Eigen::VectorXd &DXpU, Eigen::VectorXd &z) override;
    void updateAfterStep(double dt, Eigen::VectorXd &xPrev, Eigen::VectorXd &x) override;
    void copyX(Eigen::VectorXd &tar) override;
    void predictX(double dt, Eigen::VectorXd &xPrev, Eigen::VectorXd &x, Eigen::VectorXd &xBar) override;
    void buildMassMatrix(Eigen::VectorXd &m) override;
    void buildDMatrix() override;
    void buildWMatrix(double w) override;

    /* Semi-implicit method data */
    ParamIter* setUpParams();
    MatrixIter* createStiffnessMatrix();

    /* Different integrations modes */
    void eulerSolve(double dt, int elementMode, bool initMode);
    void linearImplicitSolve(double dt, int elementMode, bool initMode);
    void verletSolve(double dt, int elementMode, bool initMode);
};

}


#endif