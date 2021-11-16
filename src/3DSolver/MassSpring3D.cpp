#include "MassSpring3D.h"
#include "Pool3D.h"
#include "SolidObject3D.h"
#include <vector>
#include <assert.h>
#include <fstream>
#include <math.h>
#include "../Utils/SimUtilities.h"
#include <cmath>
#include <queue>
#include <set>
#include "../LASolver/SparseItObj.h"

using namespace mass_spring;
using namespace std;

const double EPS = 1e-14;

/**
 * Just the copy ctor
*/
MassSpring3D::MassSpring3D(const MassSpring3D &cpy) {
    // cout << "In copy CTOR" << endl;
    this->pntList = new vector<massPoint3D>(*cpy.pntList);
    this->edgeList = new vector<edge3D>(*cpy.edgeList);
    this->boundaryEdgeIdList = new vector<int>(*cpy.boundaryEdgeIdList);
    this->boundaryNodeIdList = new vector<int>(*cpy.boundaryNodeIdList);
    this->faceList = new vector<face3D>(*cpy.faceList);
    this->w = cpy.w;
    
    this->u0 = cpy.u0;
    this->v0 = cpy.v0;
    this->w0 = cpy.w0;
    
    this->gx = cpy.gx;
    this->gy = cpy.gy;
    this->gz = cpy.gz;

    /* The collision detection points */
    nodeCols = new vector<set<massPoint3D*>>(*cpy.nodeCols);
    this->admmTol = cpy.admmTol;

    this->objType = cpy.objType;

    this->repulseDist = cpy.repulseDist;
    this->collisionStiffness = cpy.collisionStiffness;
    this->collisionDist = cpy.collisionDist;

    this->updateMode = cpy.updateMode;
    this->elementMode = cpy.elementMode;

    this->dt = cpy.dt;
    this->dtFix = cpy.dtFix;

    q = new Eigen::VectorXd(*(cpy.q));
    qt = new Eigen::VectorXd(*(cpy.qt));
    qprev = new Eigen::VectorXd(*(cpy.qprev));
    qBackup = new Eigen::VectorXd(*(cpy.qBackup));
    fBackup = new Eigen::VectorXd(*(cpy.fBackup));
    qtBackup = new Eigen::VectorXd(*(cpy.qtBackup));
    qprevBackup = new Eigen::VectorXd(*(cpy.qprevBackup));
    f = new Eigen::VectorXd(*(cpy.f));

    this->structNum = cpy.structNum;

    this->eta = cpy.eta;
    this->E = cpy.E;
    this->iterCount = cpy.iterCount;
    this->mass = cpy.mass;
    this->density = cpy.density;
    this->bb_jl = cpy.bb_jl;
    this->bb_jr = cpy.bb_jr;
    this->bb_il = cpy.bb_il;
    this->bb_ir = cpy.bb_ir;
    this->bb_kl = cpy.bb_kl;
    this->bb_kr = cpy.bb_kr;
    this->pntMass = cpy.pntMass;

    // Parameters for the linear implicit method
    if (this->updateMode == 1) {
        this->matrix = createStiffnessMatrix();
        this->params = setUpParams();
        this->nnz = cpy.nnz;
        this->pcgTol = new double[3*pntList->size()];
        this->pcgRHS = new double[3*pntList->size()];
        for (int i = 0; i < 3*pntList->size(); i++) {
            this->pcgTol[i] = cpy.pcgTol[i];
            this->pcgRHS[i] = cpy.pcgRHS[i];
        }
    } else if (this->updateMode == 2) {
        this->D = new Eigen::SparseMatrix<double>(*cpy.D);
        this->W = new Eigen::SparseMatrix<double>(*cpy.W);
        this->M = new Eigen::SparseMatrix<double>(*cpy.M);

        setNPnts(pntList->size());
        setD(3);
        this->admmSolver = new ADMMPG((cpy.admmSolver)->dt, *this);
    }
    // cout << "FINSIHED In copy CTOR" << endl;
}

MassSpring3D::MassSpring3D(Pool3D &pool, int structNum, SolidObject3D &obj,
        int updateMode, int elementMode) : Assembly() {
    
    this->u0 = obj.getU0();
    this->v0 = obj.getV0();
    this->w0 = obj.getW0();

    // cout << "in CTOR" << endl;
    // Create the vectors for the edge and point lists
    pntList = new vector<massPoint3D>();
    edgeList = new vector<edge3D>();
    faceList = new vector<face3D>();
    boundaryEdgeIdList = new vector<int>();
    boundaryNodeIdList = new vector<int>();

    this->gx = pool.gx;
    this->gy = pool.gy;
    this->gz = pool.gz;

    this->structNum = structNum;

    this->elementMode = elementMode;
    this->updateMode = updateMode;

    this->objType = obj.getObjType();

    // Iterate through the Pool, adding to the vector struct as we go.
    int domain;
    objects::FSIObject poolObj;

    iterCount = 0;

    this->dt = pool.dtFix;
    this->dtFix = pool.dtFix;

    this->repulseDist = pool.repulseDist;
    this->collisionStiffness = pool.collisionStiffness;
    this->collisionDist = pool.collisionDist;

    this->mass = obj.getMass();
    this->density = obj.getDensity();
    this->eta = obj.getEta();
    this->E = obj.getE();
    this->w = sqrt(this->E);

    int pntId = 0;
    int edgeId = 0;

    double x, y, z;

    this->hx = pool.getXMeshVal(1) - pool.getXMeshVal(0);
    this->hy = pool.getYMeshVal(1) - pool.getYMeshVal(0);
    this->hz = pool.getZMeshVal(1) - pool.getZMeshVal(0);

    // bounding box indices
    bb_jl = 0; bb_jr = 0; bb_il = 0; bb_ir = 0, bb_kl = 0, bb_kr = 0;

    // Map for the object point labels (keeps track of the id of the point within each grid cell)
    map<tuple<int, int, int>, int> objectPointLabels;

    // First loop. Build up point set and HashMap for point numbers
    // int numNotFluid = 0;
    for (int k = 0; k < pool.getNz(); k++) {
        for (int j = 0; j < pool.getNy(); j++) {
            for (int i = 0; i < pool.getNx(); i++) {
                domain = pool.domainMembership(i, j, k);
                poolObj = pool.objAtIndex(i, j, k);

                if (domain == structNum && poolObj != objects::FLUID_C) {
                    // x, y, z coordinate of the current point. Interpolate in the normal direction to
                    // find the interface if this is a boundary cell
                    if (pool.isInterface(poolObj)) {
                        double pnt[3];
                        this->interpolateBoundaryLocation(pool, i, j, k, pnt);

                        // Find the inward normal direction
                        int nDir[3];
                        pool.getNormalDir(pool.objAtIndex(i, j, k), nDir);

                        nDir[0] *= -1;
                        nDir[1] *= -1;
                        nDir[2] *= -1;

                        // Coordinates of the structure point in the neighbouring direction
                        double structPnt[3] = {simutils::midpoint(pool.getXMeshVal(i+nDir[0]), pool.getXMeshVal(i+nDir[0]+1)),
                                            simutils::midpoint(pool.getYMeshVal(j+nDir[1]), pool.getYMeshVal(j+nDir[1]+1)),
                                            simutils::midpoint(pool.getZMeshVal(k+nDir[2]), pool.getZMeshVal(k+nDir[2]+1))};
                        
                        // Compute the unit outward normal
                        double outN[3] = {-((double)nDir[0]), -((double)nDir[1]), -((double)nDir[2])};
                        simutils::normalize3D(outN);

                        double diff[3] = {pnt[0] - structPnt[0], pnt[1] - structPnt[1], pnt[2] - structPnt[2]};

                        if (simutils::eps_equal(simutils::eucNorm3D(diff), 0.0, 1e-16)) {
                            // Add a small purtabation to the outward normal if two points are too close.
                            // Should not contribute to the overall error.
                            const double ADJ = 1e-16;

                            x = pnt[0] + ADJ*outN[0];
                            y = pnt[1] + ADJ*outN[1];
                            z = pnt[2] + ADJ*outN[2];
                        } else {
                            x = pnt[0];
                            y = pnt[1];
                            z = pnt[2];
                        }
                    } else {
                        x = simutils::midpoint(pool.getXMeshVal(i), pool.getXMeshVal(i+1));
                        y = simutils::midpoint(pool.getYMeshVal(j), pool.getYMeshVal(j+1));
                        z = simutils::midpoint(pool.getZMeshVal(k), pool.getZMeshVal(k+1));
                    }

                    // Add the current point to the point list
                    pntList->push_back(massPoint3D());
                    pntList->at(pntId).x = x;
                    pntList->at(pntId).y = y;
                    pntList->at(pntId).z = z;
                    pntList->at(pntId).u = obj.getU0();
                    pntList->at(pntId).v = obj.getV0();
                    pntList->at(pntId).w = obj.getW0();
                    pntList->at(pntId).structNum = structNum;
                    pntList->at(pntId).nodeId = pntId;

                    // If this is a boundary point, it should be indicated
                    pntList->at(pntId).boundaryPnt = pool.isInterface(poolObj);

                    // If it is a boundary point, add to container of boundary lists for faster looping
                    if (pntList->at(pntId).boundaryPnt) {
                        boundaryNodeIdList->push_back(pntId);
                    }

                    // Add the current label to the map
                    objectPointLabels[make_tuple(i, j, k)] = pntId;

                    // Update the bounding box
                    bb_jl = simutils::imin(bb_jl, j);
                    bb_il = simutils::imin(bb_il, i);
                    bb_kl = simutils::imin(bb_kl, k);

                    bb_jr = simutils::imax(bb_jr, j);
                    bb_ir = simutils::imax(bb_ir, i);
                    bb_kr = simutils::imax(bb_kr, k);
                    
                    pntId++;
                }

            }
        }
    }


    // Compute and then assign the point masses. Additionally, create the q vector.
    this->pntMass = obj.getMass()/pntList->size();

    q     = new Eigen::VectorXd(Eigen::VectorXd::Constant(3*pntList->size(), 0.0));
    qprev = new Eigen::VectorXd(Eigen::VectorXd::Constant(3*pntList->size(), 0.0));
    qt    = new Eigen::VectorXd(Eigen::VectorXd::Constant(3*pntList->size(), 0.0));
    qBackup     = new Eigen::VectorXd(Eigen::VectorXd::Constant(3*pntList->size(), 0.0));
    fBackup     = new Eigen::VectorXd(Eigen::VectorXd::Constant(3*pntList->size(), 0.0));
    qprevBackup = new Eigen::VectorXd(Eigen::VectorXd::Constant(3*pntList->size(), 0.0));
    qtBackup    = new Eigen::VectorXd(Eigen::VectorXd::Constant(3*pntList->size(), 0.0));
    f     = new Eigen::VectorXd(Eigen::VectorXd::Constant(3*pntList->size(), 0.0));

    int qoff = 0;
    // cout << "b4 making the pnts size = " << pntList->size() << endl;
    for (int i = 0; i < pntList->size(); i++) {
        pntList->at(i).mass = pntMass;

        // Add coordinates to the q vector.
        (*q)[qoff]   = pntList->at(i).x;
        (*q)[qoff+1] = pntList->at(i).y;
        (*q)[qoff+2] = pntList->at(i).z;

        (*qBackup)[qoff]   = pntList->at(i).x;
        (*qBackup)[qoff+1] = pntList->at(i).y;
        (*qBackup)[qoff+2] = pntList->at(i).z;

        (*fBackup)[qoff] = 0.0;
        (*fBackup)[qoff+1] = 0.0;
        (*fBackup)[qoff+2] = 0.0;

        (*qt)[qoff]   = pntList->at(i).u;
        (*qt)[qoff+1] = pntList->at(i).v;
        (*qt)[qoff+2] = pntList->at(i).w;

        (*qtBackup)[qoff]   = pntList->at(i).u;
        (*qtBackup)[qoff+1] = pntList->at(i).v;
        (*qtBackup)[qoff+2] = pntList->at(i).w;

        qoff += 3;
    }

    // cout << "Finsihed all the q's" << endl;

    // Second loop, go through the bounding box and build the edge connections
    // Note: we take advantage of the fact that we start from the bottom left,
    //       and move toward the top right.
    edgeId = 0;
    int neighId;

    double pntLoc[3];
    double diff[3];

    double hx = pool.getXMeshVal(1) - pool.getXMeshVal(0);
    double hy = pool.getYMeshVal(1) - pool.getYMeshVal(0);
    double hz = pool.getZMeshVal(1) - pool.getZMeshVal(0);

    // In this lazy version of edge addition, we search all of the 1-connected locations
    for (int k = bb_kl; k <= bb_kr; k++) {
        for (int j = bb_jl; j <= bb_jr; j++) {
            for (int i = bb_il; i <= bb_ir; i++) {
                domain = pool.domainMembership(i, j, k);
                poolObj = pool.objAtIndex(i, j, k);

                // if (domain == structNum && poolObj == objects::STRUCTURE) {
                if (domain == structNum && poolObj != objects::FLUID_C) {
                    // Get the ID of the current point.
                    pntId = objectPointLabels[make_tuple(i, j, k)];

                    // Coordinates of the current point
                    pntLoc[0] = pntList->at(pntId).x;
                    pntLoc[1] = pntList->at(pntId).y;
                    pntLoc[2] = pntList->at(pntId).z;

                    // Create arrays for i, j values that must be checked
                    int niList[13] = {i+1, i+1, i,   i-1, i,   i,   i-1, i+1, i,   i-1, i-1, i+1, i+1};
                    int njList[13] = {j,   j+1, j+1, j+1, j,   j+1, j+1, j+1, j-1, j,   j-1, j,   j-1};
                    int nkList[13] = {k,   k,   k,   k,   k+1, k+1, k+1, k+1, k+1, k+1, k+1, k+1, k+1};
                    bool sEdge[13] = {false,   true,   false,   true,   false,   true,
                                      true,   true,   true,   true,   true,   true,   true};

                    // cout << "neighId = ";
                    for (int l = 0; l < 13; l++) {
                        int ni = niList[l];
                        int nj = njList[l];
                        int nk = nkList[l];

                        // if (pool.objAtIndex(ni, nj, nk) == objects::STRUCTURE) {
                        if (pool.objAtIndex(ni, nj, nk) != objects::FLUID_C) {
                            // Extract some information about the neighbour in this direction.
                            neighId = objectPointLabels[make_tuple(ni, nj, nk)];
                            // cout << neighId << ", ";
                            if (pntId == neighId) {
                                assert(false);
                            }

                            assert(pntId != neighId);

                            // Create the new edge
                            edgeList->push_back(edge3D());
                            edgeId = edgeList->size() - 1;

                            // Compute the relaxed length of this spring, taking the norm of the difference in their coordinates
                            diff[0] = pntLoc[0] - pntList->at(neighId).x;
                            diff[1] = pntLoc[1] - pntList->at(neighId).y;
                            diff[2] = pntLoc[2] - pntList->at(neighId).z;

                            double l0;
                            l0 = simutils::eucNorm3D(diff);

                            assert(l0 != 0.0);

                            edgeList->at(edgeId).l0 = l0;
                            edgeList->at(edgeId).pntIds[0] = pntId;
                            edgeList->at(edgeId).pntIds[1] = neighId;
                            edgeList->at(edgeId).boundaryEdge = pntList->at(pntId).boundaryPnt
                                                                && pntList->at(neighId).boundaryPnt;
                            edgeList->at(edgeId).structureEdge = sEdge[l];
                            
                            // If this is a boundary edge, set the ideal edge length to 0.8*hdiag with the idea
                            // that the internal bars must "shrink" to accomodate the boundary nodes.
                            double fScale = 0.8;

                            if (!edgeList->at(edgeId).boundaryEdge) {
                                // If this is not a boundary edge, set its target length to fScale*|| [hx \\ hy \\ hz] ||
                                edgeList->at(edgeId).l0 = fScale*sqrt(simutils::square(hx) + simutils::square(hy) + simutils::square(hz));
                            }

                            // If this is a boundary edge, add a reference to it to the boundaryEdge vector
                            if (edgeList->at(edgeId).boundaryEdge) {
                                boundaryEdgeIdList->push_back(edgeId);
                            }

                            // Add this edge to the pointer list in each struct
                            pntList->at(pntId).edgeIds.push_back(edgeId);
                            pntList->at(neighId).edgeIds.push_back(edgeId);
                        }
                    }
                }
            }
        }
    }

    // cout << "Finsihed edge addition" << endl;
    // cout << "before createFaceList " << pntList->size() << endl;
    // Create the list of triangular faces to be used in the surface of the structure.
    this->createFaceList();
    // cout << "Finsihed face creation" << endl;

    // Now, attempt to update the solid to a steady state
    double eps = 1e-6;
    const int MAX_ITERS = 200;
    double dt = 0.05*simutils::dmin(hx, simutils::dmin(hy, hz));
    int iters = 0;
    // cout << "init the thing" << endl;
    double fNet[3] = {0.0, 0.0, 0.0};

    do  {
        this->updateSolidVels(dt, pool, NULL, fNet, 1, true);
        this->updateSolidLocs(pool, false);

        iters ++;
    }
    while (qt->lpNorm<1>() > eps && iters < MAX_ITERS);

    // cout << "FINISHED init the thing" << endl;

    // After this is done, set the edge length of each spring to its current config to bring the
    // potential energy to 0.
    double barVec[3];
    for (auto edge = edgeList->begin(); edge != edgeList->end(); ++edge) {
        // Vector connecting this edge.
        barVec[0] = pntList->at(edge->pntIds[1]).x - pntList->at(edge->pntIds[0]).x;
        barVec[1] = pntList->at(edge->pntIds[1]).y - pntList->at(edge->pntIds[0]).y;
        barVec[2] = pntList->at(edge->pntIds[1]).z - pntList->at(edge->pntIds[0]).z;

        // Set this bars resting length to this value
        edge->l0 = simutils::eucNorm3D(barVec);
    }

    // Set the derivatives of the locations to what they should be for the uniform initial speed
    // of the solid (positions should be correct).
    // cout << "Setting the vals" << endl;
    for (int i = 0; i < pntList->size(); i++) {
        (*q)(3*i) = pntList->at(i).x;
        (*q)(3*i+1) = pntList->at(i).y;
        (*q)(3*i+2) = pntList->at(i).z;

        (*qBackup)(3*i) = pntList->at(i).x;
        (*qBackup)(3*i+1) = pntList->at(i).y;
        (*qBackup)(3*i+2) = pntList->at(i).z;

        (*qt)(3*i)   = obj.getU0();
        (*qt)(3*i+1) = obj.getV0();
        (*qt)(3*i+2) = obj.getW0();

        (*qtBackup)(3*i)   = obj.getU0();
        (*qtBackup)(3*i+1) = obj.getV0();
        (*qtBackup)(3*i+2) = obj.getW0();

        (*fBackup)(3*i)   = 0.0;
        (*fBackup)(3*i+1) = 0.0;
        (*fBackup)(3*i+2) = 0.0;
    }

    for (auto pnt = pntList->begin(); pnt != pntList->end(); ++pnt) {
        pnt->u = obj.getU0();
        pnt->v = obj.getV0();
        pnt->w = obj.getW0();
    }

    // cout << "FINISHED Setting the vals" << endl;
    f->setZero();

    /* Set vels to zero if this is a static object */
    if (objType == SolidObject3D::ObjectType::STATIC) {
        qt->setZero();
        // return;
    }

    // If using implicit method, setup the matrix.
    if (updateMode == 1) {
        this->matrix = this->createStiffnessMatrix();

        // Set up the parameters
        this->params = setUpParams();

        pcgTol = new double[3*pntList->size()];
        pcgRHS = new double[3*pntList->size()];
        simutils::set_constant(3*pntList->size(), 0.0, pcgTol);
        simutils::set_constant(3*pntList->size(), 0.0, pcgRHS);

        matrix->set_toler(pcgTol);
    } else if (this->updateMode == 2) {
        // If using ADMM, set Assembly parameters
        
        // Set number of points necissary
        setNPnts(pntList->size());
        setD(3);

        // Build the diagonal mass matrix
        Eigen::VectorXd m(Eigen::VectorXd::Constant(3*pntList->size(), pntMass));
        buildMassMatrix(m);

        // Build the diagonal W matrix
        buildWMatrix(0.5*this->w);

        // Build the D matrix
        buildDMatrix();

        // Create the solver
        if (pool.dtFix > 0) {
            admmSolver = new ADMMPG(pool.dtFix, *this);
        } else {
            assert(false);
        }
    }

    // if (structNum != 0) {
    //     for (int i = 0; i < qt->size(); i++) {
    //         cout << (*qt)(i) << ", ";

    //     }
    //     assert(false);
    // }

    // Data structures for collision handling
    nodeCols = new vector<set<massPoint3D*>>(pntList->size());
    // cout << "Leaving the CTOR" << endl;
}

/**
 * Build the D matrix of size R^{d*nedges x d*nedges}
*/
void MassSpring3D::buildDMatrix() {
    if (dAlloc) {
        cout << "D already set!!" << endl;
        assert(false);
    }

    // Temp matrix (we insert the transpose) and we keep it this way for storage reasons
    Eigen::SparseMatrix<double> Dt(d*pntList->size(), d*edgeList->size());

    // Reserve the memory
    Dt.reserve(Eigen::VectorXi::Constant(d*edgeList->size(), d*pntList->size()));

    // Add the spring constraints to the system
    int pntIds[2];
    int rowStart0, rowStart1, colStart;
    for (int col = 0; col < edgeList->size(); col++) {
        // Extract the ids
        pntIds[0] = edgeList->at(col).pntIds[0];
        pntIds[1] = edgeList->at(col).pntIds[1];

        colStart = col*d;

        rowStart0 = pntIds[0]*d;
        rowStart1 = pntIds[1]*d;

        // First block
        for (int n = 0; n < d; n++) {
            for (int m = 0; m < d; m++) {
                if (n == m) {
                    Dt.insert(n+rowStart0, m+colStart) = -1.0;
                    Dt.insert(n+rowStart1, m+colStart) = 1.0;
                } else {
                    Dt.insert(n+rowStart0, m+colStart) = 0.0;
                    Dt.insert(n+rowStart1, m+colStart) = 0.0;
                }
            }
        }
    }

    // Create the D matrix
    D = new Eigen::SparseMatrix<double>(Dt.transpose());

    dAlloc = true;
}

/**
 * Build the M matrix of size R^{d*npnts x d*npnts}
*/
void MassSpring3D::buildMassMatrix(Eigen::VectorXd &m) {
    if (mAlloc) {
        cout << "M already set!!" << endl;
        assert(false);
    }

    // Temp matrix (we insert the transpose) and we keep it this way for storage reasons
    M = new Eigen::SparseMatrix<double>(d*pntList->size(), d*pntList->size());

    M->reserve(Eigen::VectorXd::Constant(d*nPnts, 1));

    mAlloc = true;

    for (int i = 0; i < m.size(); i++) {
        M->insert(i, i) = m[i];
    }

    mAlloc = true;
}

/**
 * Setting the tolerance for the ADMM algorithm
*/
void MassSpring3D::setAdmmTol(double tol) {
    this->admmTol = tol;
}

/**
 * Build the W matrix of size R^{d*nedges x d*nedges}
*/
void MassSpring3D::buildWMatrix(double w) {
    if (wAlloc) {
        cout << "W already set!!" << endl;
        assert(false);
    }

    // Temp matrix (we insert the transpose) and we keep it this way for storage reasons
    W = new Eigen::SparseMatrix<double>(d*edgeList->size(), d*edgeList->size());

    W->reserve(Eigen::VectorXd::Constant(d*edgeList->size(), 1));

    for (int i = 0; i < d*edgeList->size(); i++) {
        W->insert(i, i) = w;
    }

    wAlloc = true;

}

void MassSpring3D::prox(double dt, Eigen::VectorXd &x, Eigen::VectorXd &DXpU, Eigen::VectorXd &z) {
    // First, make z the projection onto the constraints
    z = DXpU;

    // Note: z_i is interpreted as a spring length D_i x_i in the typical MSS way

    // For each edge, project the length constraint : TODO: make more general constraint projection... may want to bake this into the objects,
    // specifiying constraint objects which act on the input
    for (int i = 0; i < edgeList->size(); i++) {
        // Normalize each sub-vector to the resting length of the spring
        z.segment(i*d, d) *= edgeList->at(i).l0/(z.segment(i*d, d).norm());
    }

    z = (this->E*z + simutils::square(this->w)*DXpU)/(this->E + simutils::square(this->w));
}

void MassSpring3D::updateAfterStep(double dt, Eigen::VectorXd &xPrev, Eigen::VectorXd &x) {
    // Update the node positions
    *q = x;
    
    // Update the node velocities
    *qt = (x - xPrev)/dt;
}

void MassSpring3D::copyX(Eigen::VectorXd &tar) {
    tar = *q;
}

void MassSpring3D::predictX(double dt, Eigen::VectorXd &xPrev, Eigen::VectorXd &x, Eigen::VectorXd &xBar) {
    xBar = x + dt*(*qt) + (simutils::square(dt)/pntMass)*(*f);
}

ParamIter* MassSpring3D::setUpParams() {
    ParamIter *params = new ParamIter();

    // Set the default params
    const int ILU_LEVEL = 0;

    (params)->order = 0; // = 0 natural;
                                //  =1, RCM

    (params)->level = ILU_LEVEL; // level of ilu

    (params)->drop_ilu = 0; // = 0 level ilu
                                    // = 1 drop tol

    (params)->iscal = 0; // = 0 no scaling by inverse of diag
                                // = 1 scale by inverse of diag

    (params)->nitmax = 10000; // max number of iterations

    (params)->ipiv = 0; // 0 no pivoting
                                    // 1 pivoting
                                    // only if drop_ilu = 1

    (params)->resid_reduc = 1.e-6; // residual reduction toleranc
                                // iterations stop if any of the following true
                                //   l^2 residual reduced by ctol
                                //   number of itns = nitmax
                                //   update in each variable x[i]<toler[i]

    (params)->info = 0; // = 0 no iteration info written to "output"
                    // = 1 write iteration info

    (params)->drop_tol = 1.e-3; // drop tolerance used if drop_ilu = 1

    (params)->new_rhat = 0; // = 0 use r^0 as rhat
                        // = 1 use (LU)^{-1^0 for rhat

    // (params)->iaccel = -1; // = 0 cgstab
    (params)->iaccel = 0; // = 0 cgstab
                        // = 1 orthomin
                        // = -1 conj gradiend (SPD only!)

    (params)->north = 10; // number of orthogs for orthomin

    return params;

}

/**
 * Create the stiffness matrix for the mass-spring system when using semi-implicit methods.
*/
MatrixIter* MassSpring3D::createStiffnessMatrix() {
    const int d = 3;

    // Extract the number of points in the MSS
    int n = pntList->size();

    // Sum of the degrees of all points
    int sumDeg = 0;
    for (auto pnt = pntList->begin(); pnt != pntList->end(); ++pnt) {
        sumDeg += pnt->edgeIds.size();
    }

    // Compute the number of non-zero elements of the matrix
    this->nnz = simutils::square(d)*(sumDeg + n);

    // Now setup the information for the matrix structure
    MatrixStruc *matrixBuilder = NULL;

    try {
        matrixBuilder = new MatrixStruc(3*pntList->size() /*Number of unknowns*/,
                                        true /*Don't add non-zeros to diagonals so I don't get mindflooded*/);

        /* Create the data structure */
        
        // Loop through the point ids
        int neighId;
        int pntId0, pntId1, off, noff;
        for (int i = 0; i < pntList->size(); i++) {

            // Set the entry for this point (just hardcode for brain reasons)

            off = i*d;

            matrixBuilder->set_entry(off, off);
            matrixBuilder->set_entry(off, off+1);
            matrixBuilder->set_entry(off, off+2);

            matrixBuilder->set_entry(off+1, off);
            matrixBuilder->set_entry(off+1, off+1);
            matrixBuilder->set_entry(off+1, off+2);

            matrixBuilder->set_entry(off+2, off);
            matrixBuilder->set_entry(off+2, off+1);
            matrixBuilder->set_entry(off+2, off+2);

            // Loop through each of the neighbouring edges.
            for (auto edgeId = pntList->at(i).edgeIds.begin(); edgeId != pntList->at(i).edgeIds.end(); ++edgeId) {
                // Extract the ID of the neighbouring point
                pntId0 = edgeList->at(*edgeId).pntIds[0];
                pntId1 = edgeList->at(*edgeId).pntIds[1];
                neighId = (pntId0 == i) ? pntId1 : pntId0;

                // Set the entries
                noff = neighId*d;

                matrixBuilder->set_entry(off, noff);
                matrixBuilder->set_entry(off, noff+1);
                matrixBuilder->set_entry(off, noff+2);

                matrixBuilder->set_entry(off+1, noff);
                matrixBuilder->set_entry(off+1, noff+1);
                matrixBuilder->set_entry(off+1, noff+2);

                matrixBuilder->set_entry(off+2, noff);
                matrixBuilder->set_entry(off+2, noff+1);
                matrixBuilder->set_entry(off+2, noff+2);
            }
        }

        // Pack the data constructor: convert the LinkedList array to csr format
        matrixBuilder->pack();

    } catch(bad_alloc) {
        delete matrixBuilder; matrixBuilder = NULL;
        exit(1);
    } catch(General_Exception excep) {
        cout << "General exception" << endl;
        cout << excep.p << endl;
        delete matrixBuilder; matrixBuilder = NULL;
        exit(1);
    }

    // Create the actual matrix
    MatrixIter *matrix = NULL;

    try {
        matrix = new MatrixIter(*matrixBuilder);
        assert(matrix != NULL);
        delete matrixBuilder; 
        matrixBuilder = NULL;
    } catch(bad_alloc) {
        delete matrixBuilder; matrixBuilder = NULL;
        exit(1);
    } catch(General_Exception excep) {
        delete matrixBuilder; matrixBuilder= NULL;
        exit(1);
    } 

    return matrix;
}

/**
 * General elastic force calculation, as we use this in several places in the code
*/
void MassSpring3D::calcElasticForce(double E, double l0, massPoint3D pnt1,
                            massPoint3D pnt2, double force[6]) {
                                    
    double diff[3] = {pnt2.x - pnt1.x, pnt2.y - pnt1.y, pnt2.z - pnt1.z};
    double diffNorm = simutils::eucNorm3D(diff);

    force[0] = -E*((diffNorm - l0)/diffNorm)*(-diff[0]);
    force[1] = -E*((diffNorm - l0)/diffNorm)*(-diff[1]);
    force[2] = -E*((diffNorm - l0)/diffNorm)*(-diff[2]);

    force[3] = -E*((diffNorm - l0)/diffNorm)*(diff[0]);
    force[4] = -E*((diffNorm - l0)/diffNorm)*(diff[1]);
    force[5] = -E*((diffNorm - l0)/diffNorm)*(diff[2]);
}

/**
 * Computes the collision stress for a given boundary node, specified by its ID
*/
bool MassSpring3D::computeCollisionStress(int nodeId, double colStress[3], double dA) {
    // 0 out the stress vector
    colStress[0] = 0.0;
    colStress[1] = 0.0;
    colStress[2] = 0.0;

    // Compute the forces acting internally on this node by the other members of the MSS
    massPoint3D mPnt = pntList->at(nodeId);

    int neighId = -1;
    int *pntIds;
    double forces[6];
    double nodeForces[3] = {0.0, 0.0, 0.0};
    for (auto edgeId = mPnt.edgeIds.begin(); edgeId != mPnt.edgeIds.end(); ++edgeId) {
        // Find the id of the connected node
        pntIds = edgeList->at(*edgeId).pntIds;
        neighId = ( pntIds[0] == nodeId ) ? pntIds[1] : pntIds[0];

        // Calculate the elastic force
        calcElasticForce(this->E, edgeList->at(*edgeId).l0,
                    mPnt, pntList->at(neighId), forces);

        // Calculate the force being applied to the boundary node
        nodeForces[0] += forces[0];
        nodeForces[1] += forces[1];
        nodeForces[2] += forces[2];
    }

    // double f_i[3] = {(*fBackup)(3*nodeId), (*fBackup)(3*nodeId+1), (*fBackup)(3*nodeId+2)};
    double v_i[3] = {(*qt)(3*nodeId), (*qt)(3*nodeId+1), (*qt)(3*nodeId+2)};

    // Calculate the force to cancel the velocities
    double cancelStress[3] = {
        - (pntMass/dt)*v_i[0],
        - (pntMass/dt)*v_i[1],
        - (pntMass/dt)*v_i[2]
    };

    // For each of the collisions this node is involved in, calculate a resisting force
    massPoint3D colPnt;
    double pntDiff[3];
    double pntDist;
    bool colComp = false;
    int numNear = 0;
    for (auto near = nodeCols->at(nodeId).begin(); near != nodeCols->at(nodeId).end(); ++near) {
        colPnt = **near;

        // Compute the spring repulsion force for this node
        pntDiff[0] = mPnt.x - colPnt.x;
        pntDiff[1] = mPnt.y - colPnt.y;
        pntDiff[2] = mPnt.z - colPnt.z;
        pntDist = simutils::eucNorm3D(pntDiff);

        if (repulseDist > pntDist) {
            calcElasticForce(this->collisionStiffness, repulseDist, mPnt, colPnt, forces);
            numNear++;
            colComp |= true;
        } else {
            forces[0] = 0.0;
            forces[1] = 0.0;
            forces[2] = 0.0;
            colComp |= false;
        }

        colStress[0] += forces[0];
        colStress[1] += forces[1];
        colStress[2] += forces[2];
    }

    colStress[0] /= numNear;
    colStress[1] /= numNear;
    colStress[2] /= numNear;

    // Cancel stress ignores scaling
    cancelStress[0] /= dA;
    cancelStress[1] /= dA;
    cancelStress[2] /= dA;

    // Apply the velocity stop stress
    colStress[0] += cancelStress[0]; //- nodeForces[0];
    colStress[1] += cancelStress[1]; //- nodeForces[1];
    colStress[2] += cancelStress[2]; //- nodeForces[2];

    // cout << "ColStress: U = " <<  colStress[0] << " V = " << colStress[1] << " W = " << colStress[2] << endl;
    // cout << "colComp? " << colComp << endl;

    return colComp;
}

/** 
 * Helper function to apply external force to a boundary node.
 * 
*/
void MassSpring3D::applyBoundaryForces(Pool3D &pool, double ****stress, int ng, double fNet[3], double colNet[3]) {

    int ni, nj, nk; // Locations of the nodes

    massPoint3D mPnt;

    for (auto pnt = pntList->begin(); pnt != pntList->end(); ++pnt) {
        pnt->sigU = 0.0;
        pnt->sigV = 0.0;
        pnt->sigW = 0.0;
    }

    int iPnt[3];

    massPoint3D mPnt1, mPnt2, mPnt3;

    // Use the integral of the nodal forces on each triangular element using 2nd order
    // (midpoint method) using evaluations at the centroid.
    double a[3];
    double b[3];
    double aXB[3];
    double centroid[3];
    double baryCoords[3];
    double dA, alph, beta, gamm;

    int id1, id2, id3;
    double s1[3];
    double s2[3];
    double s3[3];
    for (auto face = faceList->begin(); face != faceList->end(); ++face) {
        // Extract the points from this triangle
        id1 = face->pntIds[0];
        id2 = face->pntIds[1];
        id3 = face->pntIds[2];

        mPnt1 = pntList->at(id1);
        mPnt2 = pntList->at(id2);
        mPnt3 = pntList->at(id3);

        /* Compute the area of this triangle */

        // a = pnt2 - pnt1
        a[0] = mPnt2.x - mPnt1.x;
        a[1] = mPnt2.y - mPnt1.y;
        a[2] = mPnt2.z - mPnt1.z;

        // b = pnt3 - pnt1
        b[0] = mPnt3.x - mPnt1.x;
        b[1] = mPnt3.y - mPnt1.y;
        b[2] = mPnt3.z - mPnt1.z;

        // Take the cross product of these vectors
        simutils::cross_product_3D(a, b, aXB);

        // The area of the triangle is 0.5*||a x b||_2
        dA = simutils::eucNorm3D(aXB)/2.0;

        // Find stress acting on first node
        bool found1 = findNearestGridNormalInterface(pool, id1, iPnt);
        ni = iPnt[0];
        nj = iPnt[1];
        nk = iPnt[2];

        if (found1) {
            for (int i = 0; i < 3; i++) {
                s1[i] = 0.0;
            }

            // Apply the stress to this point and its connected neighbours
            if (nodeCols->at(id1).size() > 0) {
                // There is a collision on this node, compute the collision stress
                bool colApplied = computeCollisionStress(id1, s1, dA);

                if (!colApplied) {
                    s1[0] = stress[0][ng+nk][ng+nj][ng+ni];
                    s1[1] = stress[1][ng+nk][ng+nj][ng+ni];
                    s1[2] = stress[2][ng+nk][ng+nj][ng+ni];
                }

            } else {
                // assert(false);
                // Apply hydrodynamic stress if there is no collision
                s1[0] = stress[0][ng+nk][ng+nj][ng+ni];
                s1[1] = stress[1][ng+nk][ng+nj][ng+ni];
                s1[2] = stress[2][ng+nk][ng+nj][ng+ni];
            }
        } else {
            assert(false);
            s1[0] = 0.0;
            s1[1] = 0.0;
            s1[2] = 0.0;
        }

        bool found2 = findNearestGridNormalInterface(pool, id2, iPnt);
        ni = iPnt[0];
        nj = iPnt[1];
        nk = iPnt[2];

        if (found2) {
            for (int i = 0; i < 3; i++) {
                s2[i] = 0.0;
            }

            // Apply the stress to this point and its connected neighbours
            if (nodeCols->at(id2).size() > 0) {
                // There is a collision on this node, compute the collision stress
                bool colApplied = computeCollisionStress(id2, s2, dA);

                if (!colApplied) {
                    s2[0] = stress[0][ng+nk][ng+nj][ng+ni];
                    s2[1] = stress[1][ng+nk][ng+nj][ng+ni];
                    s2[2] = stress[2][ng+nk][ng+nj][ng+ni];
                }
            } else {
                // Apply hydrodynamic stress if there is no collision
                s2[0] = stress[0][ng+nk][ng+nj][ng+ni];
                s2[1] = stress[1][ng+nk][ng+nj][ng+ni];
                s2[2] = stress[2][ng+nk][ng+nj][ng+ni];
            }
        } else {
            assert(false);
            s2[0] = 0.0;
            s2[1] = 0.0;
            s2[2] = 0.0;
        }

        bool found3 = findNearestGridNormalInterface(pool, id3, iPnt);
        ni = iPnt[0];
        nj = iPnt[1];
        nk = iPnt[2];

        if (found3) {
            for (int i = 0; i < 3; i++) {
                s3[i] = 0.0;
            }

            // Apply the stress to this point and its connected neighbours
            if (nodeCols->at(id3).size() > 0) {
                // There is a collision on this node, compute the collision stress
                bool colApplied = computeCollisionStress(id3, s3, dA);

                if (!colApplied) {
                    s3[0] = stress[0][ng+nk][ng+nj][ng+ni];
                    s3[1] = stress[1][ng+nk][ng+nj][ng+ni];
                    s3[2] = stress[2][ng+nk][ng+nj][ng+ni];
                }
            } else {
                // Apply hydrodynamic stress if there is no collision
                s3[0] = stress[0][ng+nk][ng+nj][ng+ni];
                s3[1] = stress[1][ng+nk][ng+nj][ng+ni];
                s3[2] = stress[2][ng+nk][ng+nj][ng+ni];
            }
        } else {
            assert(false);
            s3[0] = 0.0;
            s3[1] = 0.0;
            s3[2] = 0.0;
        }

        // Compute the centroid of this triangle
        // computeCentroid(face->pntIds[0], face->pntIds[1], face->pntIds[2], centroid);

        // // Find the barycentric coordinates of the centroid.
        // projTriangleDist(centroid, face->pntIds[0], face->pntIds[1], face->pntIds[2], baryCoords);

        // // Extract the coordinates
        // alph = baryCoords[0];
        // beta = baryCoords[1];
        // gamm = baryCoords[2];

        // Update the stresses in the MSS
        // cout << "dA = " << dA << endl;
        // double s_x = (alph*s1[0] + beta*s2[0] + gamm*s3[0])*dA;
        // double s_y = (alph*s1[1] + beta*s2[1] + gamm*s3[1])*dA;
        // double s_z = (alph*s1[2] + beta*s2[2] + gamm*s3[2])*dA;

        double s_x = (s1[0] + s2[0] + s3[0])*(dA/3.0);
        double s_y = (s1[1] + s2[1] + s3[1])*(dA/3.0);
        double s_z = (s1[2] + s2[2] + s3[2])*(dA/3.0);

        (*f)[3*id1]   +=  s_x;
        (*f)[3*id1+1] +=  s_y;
        (*f)[3*id1+2] +=  s_z;

        (*f)[3*id2]   +=  s_x;
        (*f)[3*id2+1] +=  s_y;
        (*f)[3*id2+2] +=  s_z;

        (*f)[3*id3]   +=  s_x;
        (*f)[3*id3+1] +=  s_y;
        (*f)[3*id3+2] +=  s_z;
    }

    for (int id = 0; id < pntList->size(); id++) {
        pntList->at(id).sigU = (*f)[3*id];
        pntList->at(id).sigV = (*f)[3*id+1];
        pntList->at(id).sigW = (*f)[3*id+2];
    }

    // Compute the net force acting on the edges
    fNet[0] = 0.0;
    fNet[1] = 0.0;
    fNet[2] = 0.0;

    for (auto face = faceList->begin(); face != faceList->end(); ++face) {
        // Extract the points from this triangle
        mPnt1 = pntList->at(face->pntIds[0]);
        mPnt2 = pntList->at(face->pntIds[1]);
        mPnt3 = pntList->at(face->pntIds[2]);

        /* Compute the area of this triangle */

        // a = pnt2 - pnt1
        // cout << "Computing the cross product" << endl;
        a[0] = mPnt2.x - mPnt1.x;
        a[1] = mPnt2.y - mPnt1.y;
        a[2] = mPnt2.z - mPnt1.z;

        // b = pnt3 - pnt1
        b[0] = mPnt3.x - mPnt1.x;
        b[1] = mPnt3.y - mPnt1.y;
        b[2] = mPnt3.z - mPnt1.z;

        // Take the cross product of these vectors
        simutils::cross_product_3D(a, b, aXB);

        // The area of the triangle is 0.5*||a x b||_2
        dA = simutils::eucNorm3D(aXB)/2.0;

        // Compute the centroid of this triangle
        // computeCentroid(face->pntIds[0], face->pntIds[1], face->pntIds[2], centroid);

        // Find the barycentric coordinates for this triangle.
        // projTriangleDist(centroid, face->pntIds[0], face->pntIds[1], face->pntIds[2], baryCoords);

        // // Extract the coordinates
        // alph = baryCoords[0];
        // beta = baryCoords[1];
        // gamm = baryCoords[2];

        // Interpolate at the centroid
        fNet[0] += (mPnt1.sigU + mPnt2.sigU + mPnt3.sigU)*(dA/3.0);
        fNet[1] += (mPnt1.sigV + mPnt2.sigV + mPnt3.sigV)*(dA/3.0);
        fNet[2] += (mPnt1.sigW + mPnt2.sigW + mPnt3.sigW)*(dA/3.0);
    }

    // cout << "Net force detected in compute Forces:" << endl;
    // cout << "U = " << fNet[0] << " V = " << fNet[1] << " W " << fNet[2] << endl;
    // cout << "collisionStiffness = " << collisionStiffness << endl;
    // assert(false);
}

/**
 * Add to the local potential energy for simple spring system. (note, when summing everything
 * together is when we apply the (-) sign)!
 * 
 * Note: must have pntId1 < pntId2
*/
void MassSpring3D::calcLocalElasticForce(edge3D edge, int pntId1, massPoint3D pnt1, int pntId2, massPoint3D pnt2) {
    assert(pntId1 < pntId2);
    double force[6];

    calcElasticForce(E, edge.l0, pnt1, pnt2, force);

    // Compute the force for the first input point (note we must take the sum)
    (*f)(3*pntId1)   += force[0];
    (*f)(3*pntId1+1) += force[1];
    (*f)(3*pntId1+2) += force[2];

    (*f)(3*pntId2)   += force[3];
    (*f)(3*pntId2+1) += force[4];
    (*f)(3*pntId2+2) += force[5];
}

/**
 * Add to the local potential energy for spring-dashpot system. (note, when summing everything
 * together is when we apply the (-) sign)!
 * 
 * Note: must have pntId1 < pntId2
*/
void MassSpring3D::calcLocalKelvinForce(edge3D edge, int pntId1, massPoint3D pnt1, int pntId2, massPoint3D pnt2) {
    assert(pntId1 < pntId2);

    double diff[3]  = {pnt2.x - pnt1.x, pnt2.y - pnt1.y, pnt2.z - pnt1.z};
    double diffU[3] = {pnt2.u - pnt1.u, pnt2.v - pnt1.v, pnt2.w - pnt1.w};

    double diffNorm = simutils::eucNorm3D(diff);
    double diffNormCube = diffNorm*diffNorm*diffNorm;
    double innerProd = -pnt1.x*(pnt2.u - pnt1.u) - pnt1.y*(pnt2.v - pnt1.v) - pnt1.z*(pnt2.w - pnt1.w)
                    + pnt2.x*(pnt2.u - pnt1.u) + pnt2.y*(pnt2.v - pnt1.v) + pnt2.z*(pnt2.w - pnt1.w);
    
    if (simutils::eps_equal(diffNorm, 0.0, EPS)) {
        cout << "ERROR: calculating force for MSS" << endl;
        cout << "Points with ids " << pntId1 << " and " << pntId2 << " are too close!" << endl;
        cout << "diffNorm = " << diffNorm << endl;
        assert(false);
    }

    // Compute the force for the first input point (note we must take the sum)
    (*f)[3*pntId1] += -(E*((diffNorm - edge.l0)/diffNorm)*(-diff[0]) 
        + eta*( (innerProd/simutils::square(diffNorm))*(-diff[0])
                    + (diffNorm - edge.l0)/diffNorm*(-diffU[0])
                    - (diffNorm - edge.l0)*innerProd/diffNormCube*(-diff[0]) ));
    (*f)[3*pntId1+1] += -(E*((diffNorm - edge.l0)/diffNorm)*(-diff[1]) 
        + eta*( (innerProd/simutils::square(diffNorm))*(-diff[1])
                    + (diffNorm - edge.l0)/diffNorm*(-diffU[1])
                    - (diffNorm - edge.l0)*innerProd/diffNormCube*(-diff[1]) ));
    (*f)[3*pntId1+2] += -(E*((diffNorm - edge.l0)/diffNorm)*(-diff[2]) 
        + eta*( (innerProd/simutils::square(diffNorm))*(-diff[2])
                    + (diffNorm - edge.l0)/diffNorm*(-diffU[2])
                    - (diffNorm - edge.l0)*innerProd/diffNormCube*(-diff[2]) ));

    (*f)[3*pntId2] += -(E*((diffNorm - edge.l0)/diffNorm)*(diff[0]) 
        + eta*( (innerProd/simutils::square(diffNorm))*(diff[0])
                    + (diffNorm - edge.l0)/diffNorm*(diffU[0])
                    - (diffNorm - edge.l0)*innerProd/diffNormCube*(diff[0]) ));
    (*f)[3*pntId2+1] += -(E*((diffNorm - edge.l0)/diffNorm)*(diff[1]) 
        + eta*( (innerProd/simutils::square(diffNorm))*(diff[1])
                    + (diffNorm - edge.l0)/diffNorm*(diffU[1])
                    - (diffNorm - edge.l0)*innerProd/diffNormCube*(diff[1]) ));
    (*f)[3*pntId2+2] += -(E*((diffNorm - edge.l0)/diffNorm)*(diff[2]) 
        + eta*( (innerProd/simutils::square(diffNorm))*(diff[2])
                    + (diffNorm - edge.l0)/diffNorm*(diffU[2])
                    - (diffNorm - edge.l0)*innerProd/diffNormCube*(diff[2]) ));
}

/**
 * For a border cell (i, j), find the location of the interface using normal interpolation.
 * 
 * TODO: consider using the gradient of the level set function to find this.
 * 
 * The approach: find the structural cell that is in the normal direction of this boundary.
 *               Once this has been found, we find the vector pointing from the outside 
 *               point center to the inside point center. Solving for the root of the interpolating
 *               polynomial connecting these two points, in the direction of this vector gives an
 *               approximation to the interface location.
*/
void MassSpring3D::interpolateBoundaryLocation(Pool3D &pool, int i, int j, int k, double X[3]) {
    // Get the outward normal direction for this point
    int nDir[3];
    pool.getNormalDir(pool.objAtIndex(i, j, k), nDir);

    // Negate to find the inward normal
    nDir[0] *= -1;
    nDir[1] *= -1;
    nDir[2] *= -1;

    // Coordinate of the current boundary point
    double bndPnt[3] = {simutils::midpoint(pool.getXMeshVal(i), pool.getXMeshVal(i+1)),
                        simutils::midpoint(pool.getYMeshVal(j), pool.getYMeshVal(j+1)),
                        simutils::midpoint(pool.getZMeshVal(k), pool.getZMeshVal(k+1))};
    
    // Coordinates of the structure point in the neighbouring direction
    double structPnt[3] = {simutils::midpoint(pool.getXMeshVal(i+nDir[0]), pool.getXMeshVal(i+nDir[0]+1)),
                        simutils::midpoint(pool.getYMeshVal(j+nDir[1]), pool.getYMeshVal(j+nDir[1]+1)),
                        simutils::midpoint(pool.getZMeshVal(k+nDir[2]), pool.getZMeshVal(k+nDir[2]+1))};
    
    // Compute the unit inward normal
    double inN[3] = {(double)nDir[0], (double)nDir[1], (double)nDir[2]};
    simutils::normalize3D(inN);

    // Take the phi values at these points
    double phiOut = pool.getPhiVal(i, j, k);
    double phiIn = pool.getPhiVal(i+nDir[0], j+nDir[1], k+nDir[2]);


    // if (abs(phiOut - phiIn) < 1e-12) {
    //     cout << "diff = " << abs(phiOut - phiIn)  << endl;
    //     assert(abs(phiOut - phiIn) < 1e-12);

    // }
    // cout << abs(phiOut - phiIn) << endl;
    // cout << (abs(phiOut - phiIn) < 1e-12) << endl;
    // assert((abs(phiOut - phiIn) > 1e-12));

    // Compute the distance to the interface
    double distVec[3] = {bndPnt[0] - structPnt[0], bndPnt[1] - structPnt[1], bndPnt[2] - structPnt[2]};
    double d = abs(phiOut/(phiIn - phiOut))*simutils::eucNorm3D(distVec);
    // cout << "phiIn = " << phiIn << ", phiOut = " << phiOut << endl;
    // cout << "||distVec|| = " << simutils::eucNorm3D(distVec) << endl;
    // cout << "d = " << d << endl;

    // Finally, compute the interface point location
    X[0] = bndPnt[0] + d*inN[0];
    X[1] = bndPnt[1] + d*inN[1];
    X[2] = bndPnt[2] + d*inN[2];
}

/**
 * Update the locations and velocities of the spring nodes using Euler's method
*/
void MassSpring3D::eulerSolve(double dt, int elementMode, bool initMode) {
    int pntId1;
    int pntId2;
    massPoint3D pnt1;
    massPoint3D pnt2;

    double ETemp = E;
    double etaTemp = eta;
    
    if (initMode) {
        E = 5.0;
        eta = 0.0;
    }

    // Compute the force vector for the current configuration
    for (auto edge = edgeList->begin(); edge != edgeList->end(); ++edge) {
        // Extract the points connecting the current spring
        pntId1 = edge->pntIds[0];
        pntId2 = edge->pntIds[1];

        pnt1 = pntList->at(pntId1);
        pnt2 = pntList->at(pntId2);

        // Calculate the contribution of the current point to the total force vector.
        if (elementMode == 0) {
            calcLocalElasticForce(*edge, pntId1, pnt1, pntId2, pnt2);
        } else {
            calcLocalKelvinForce(*edge, pntId1, pnt1, pntId2, pnt2);
        }
    }

    // Now, the force vector has been obtained, we perform an Euler step
    for (int i = 0; i < pntList->size(); i++) {
        if (!(initMode && pntList->at(i).boundaryPnt)) {
            (*qt)(3*i) += (dt/pntMass)*(*f)(3*i);
            (*qt)(3*i+1) += (dt/pntMass)*(*f)(3*i+1);
            (*qt)(3*i+2) += (dt/pntMass)*(*f)(3*i+2);

            (*q)(3*i)  += dt*(*qt)(3*i);
            (*q)(3*i+1)  += dt*(*qt)(3*i+1);
            (*q)(3*i+2)  += dt*(*qt)(3*i+2);
        } else {
            (*qt)(3*i)   = 0.0;
            (*qt)(3*i+1) = 0.0;
            (*qt)(3*i+2) = 0.0;
        }
    }

    E = ETemp;
    eta = etaTemp;
}

/**
 * Add to the Hessian for simple spring system
 * 
 * Notes are on page 24 of fund var mesh adapt
 * 
 * // NOTE: assumes that the mass per node is constant (good assumption)
*/
void MassSpring3D::calcLocalElasticHessian(double dt, edge3D edge, int pntId1,
                                        massPoint3D pnt1, int pntId2,
                                        massPoint3D pnt2) {
    assert(pntId1 < pntId2);

    double diff[3] = {pnt2.x - pnt1.x, pnt2.y - pnt1.y, pnt2.z - pnt1.z};
    double diffNorm = simutils::eucNorm3D(diff);

    int rOff1 = 3*pntId1;
    int rOff2 = 3*pntId2;

    double coeff1 = E*edge.l0/(simutils::cube(diffNorm));
    double coeff2 = E*((diffNorm - edge.l0)/(diffNorm));
    int colIndex;

    // r = i(0), c = i(0) & r = i(0) c = i(1)
    for (int i = matrix->rowBegin(rOff1); i < matrix->rowEndPlusOne(rOff1); i++) {
        colIndex = matrix->getColIndex(i);

        if (colIndex == rOff1) {
            matrix->aValue(i) += coeff1*diff[0]*diff[0] + coeff2;
        } else if (colIndex == rOff1+1) {
            matrix->aValue(i) += coeff1*diff[0]*diff[1];
        } else if (colIndex == rOff1+2) {
            matrix->aValue(i) += coeff1*diff[0]*diff[2];
        } else if (colIndex == rOff2) {
            matrix->aValue(i) += coeff1*-diff[0]*diff[0] - coeff2;
        } else if (colIndex == rOff2+1) {
            matrix->aValue(i) += coeff1*-diff[0]*diff[1];
        } else if (colIndex == rOff2+2) {
            matrix->aValue(i) += coeff1*-diff[0]*diff[2];
        }
    }

    for (int i = matrix->rowBegin(rOff1+1); i < matrix->rowEndPlusOne(rOff1+1); i++) {
        colIndex = matrix->getColIndex(i);

        if (colIndex == rOff1) {
            matrix->aValue(i) += coeff1*diff[1]*diff[0];
        } else if (colIndex == rOff1+1) {
            matrix->aValue(i) += coeff1*diff[1]*diff[1] + coeff2;
        } else if (colIndex == rOff1+2) {
            matrix->aValue(i) += coeff1*diff[1]*diff[2];
        } else if (colIndex == rOff2) {
            matrix->aValue(i) += coeff1*-diff[1]*diff[0];
        } else if (colIndex == rOff2+1) {
            matrix->aValue(i) += coeff1*-diff[1]*diff[1] - coeff2;
        } else if (colIndex == rOff2+1) {
            matrix->aValue(i) += coeff1*-diff[1]*diff[2];
        }
    }

    for (int i = matrix->rowBegin(rOff1+2); i < matrix->rowEndPlusOne(rOff1+2); i++) {
        colIndex = matrix->getColIndex(i);

        if (colIndex == rOff1) {
            matrix->aValue(i) += coeff1*diff[2]*diff[0];
        } else if (colIndex == rOff1+1) {
            matrix->aValue(i) += coeff1*diff[2]*diff[1];
        } else if (colIndex == rOff1+2) {
            matrix->aValue(i) += coeff1*diff[2]*diff[2] + coeff2;
        } else if (colIndex == rOff2) {
            matrix->aValue(i) += coeff1*-diff[2]*diff[0];
        } else if (colIndex == rOff2+1) {
            matrix->aValue(i) += coeff1*-diff[2]*diff[1];
        } else if (colIndex == rOff2+1) {
            matrix->aValue(i) += coeff1*-diff[2]*diff[2] - coeff2;;
        }
    }

    // r = i(1), c = i(0) & r = i(1), c = i(1)
    for (int i = matrix->rowBegin(rOff2); i < matrix->rowEndPlusOne(rOff2); i++) {
        colIndex = matrix->getColIndex(i);

        if (colIndex == rOff1) {
            matrix->aValue(i) += coeff1*diff[0]*-diff[0] - coeff2;
        } else if (colIndex == rOff1+1) {
            matrix->aValue(i) += coeff1*diff[0]*-diff[1];
        } else if (colIndex == rOff1+2) {
            matrix->aValue(i) += coeff1*diff[0]*-diff[2];
        } else if (colIndex == rOff2) {
            matrix->aValue(i) += coeff1*diff[0]*diff[0] + coeff2;
        } else if (colIndex == rOff2+1) {
            matrix->aValue(i) += coeff1*diff[0]*diff[1];
        } else if (colIndex == rOff2+2) {
            matrix->aValue(i) += coeff1*diff[0]*diff[2];
        }
    }

    for (int i = matrix->rowBegin(rOff2+1); i < matrix->rowEndPlusOne(rOff2+1); i++) {
        colIndex = matrix->getColIndex(i);

        if (colIndex == rOff1) {
            matrix->aValue(i) += coeff1*diff[1]*-diff[0];
        } else if (colIndex == rOff1+1) {
            matrix->aValue(i) += coeff1*diff[1]*-diff[1] - coeff2;
        } else if (colIndex == rOff1+2) {
            matrix->aValue(i) += coeff1*diff[1]*-diff[2];
        } else if (colIndex == rOff2) {
            matrix->aValue(i) += coeff1*diff[1]*diff[0];
        } else if (colIndex == rOff2+1) {
            matrix->aValue(i) += coeff1*diff[1]*diff[1] + coeff2;
        } else if (colIndex == rOff2+2) {
            matrix->aValue(i) += coeff1*diff[1]*diff[2];
        }
    }

    for (int i = matrix->rowBegin(rOff2+2); i < matrix->rowEndPlusOne(rOff2+2); i++) {
        colIndex = matrix->getColIndex(i);

        if (colIndex == rOff1) {
            matrix->aValue(i) += coeff1*diff[2]*-diff[0];
        } else if (colIndex == rOff1+1) {
            matrix->aValue(i) += coeff1*diff[2]*-diff[1];
        } else if (colIndex == rOff1+2) {
            matrix->aValue(i) += coeff1*diff[2]*-diff[2] - coeff2;
        } else if (colIndex == rOff2) {
            matrix->aValue(i) += coeff1*diff[2]*diff[0];
        } else if (colIndex == rOff2+1) {
            matrix->aValue(i) += coeff1*diff[2]*diff[1];
        } else if (colIndex == rOff2+2) {
            matrix->aValue(i) += coeff1*diff[2]*diff[2] + coeff2;
        }
    }
}

/**
 * Update according to the linearly implicit time stepping method
 * 
 * Only really going to deal with this when it becomes nesissary. Would rather focus on simpler things 
 * for now.
*/
void MassSpring3D::linearImplicitSolve(double dt, int elementMode, bool initMode) {
    int pntId1;
    int pntId2;
    massPoint3D pnt1;
    massPoint3D pnt2;

    double ETemp = E;
    double etaTemp = eta;
    
    if (initMode) {
        E = 5.0;
        eta = 0.0;
    }

    *qprev = *qt;

    // 0 out the Hessian
    for (int r = 0; r < 3*pntList->size(); r++) {
        for (int i = matrix->rowBegin(r); i < matrix->rowEndPlusOne(r); i++) {
            matrix->aValue(i) = 0.0;
        }
    }

    // Initialize the linear system RHS
    simutils::set_constant(3*pntList->size(), 0.0, pcgRHS);

    // Compute the force vector and Hessian matrix for the current configuration
    for (auto edge = edgeList->begin(); edge != edgeList->end(); ++edge) {
        // Extract the points connecting the current spring
        pntId1 = edge->pntIds[0];
        pntId2 = edge->pntIds[1];

        pnt1 = pntList->at(pntId1);
        pnt2 = pntList->at(pntId2);

        // Calculate the contribution of the current point to the total force vector.
        if (elementMode == 0) {
            calcLocalElasticForce(*edge, pntId1, pnt1, pntId2, pnt2);
            calcLocalElasticHessian(dt, *edge, pntId1, pnt1, pntId2, pnt2);
        } else {
            cout << "Have not finished the Hessian for the Kelvin force" << endl;
            assert(false);
            calcLocalKelvinForce(*edge, pntId1, pnt1, pntId2, pnt2);
        }
    }

    /* Setup the remainder of the matrix by negating K everywhere, and adding I along the main diagonal */
    double coeff = simutils::square(dt)/pntList->at(0).mass;
    int colIndex;

    for (int r = 0; r < 3*pntList->size(); r++) {
        for (int i = matrix->rowBegin(r); i < matrix->rowEndPlusOne(r); i++) {
            colIndex = matrix->getColIndex(i);

            if (colIndex == r) {
                matrix->aValue(i) = 1 + coeff*matrix->aValue(i);
            } else {
                matrix->aValue(i) = coeff*matrix->aValue(i);
            }
        }
    }

    /* Setup the RHS vector (use Euler step as initial guess) */
    for (int i = 0; i < 3*pntList->size(); i++) {
        matrix->bValue(i) = (*qt)[i] + (dt/(pntList->at(0).mass))*(*f)[i] + 1e-16;
        pcgRHS[i] = matrix->bValue(i);
    }

    /* Solve the linear system */
    try {
        int niter=0; // The number of iterations done in the solve

        // Perform the ILU factorization
        if (params->drop_ilu == 0) {
            matrix->sfac(*params);
        }

        simutils::set_constant(3*pntList->size(), 0.0, pcgTol);

        matrix->set_toler(this->pcgTol); // Set the update tolerance

        // Solve the linear system. Solution is stored in provided vector
        matrix->solve(*this->params, pcgRHS, niter);

        cout << "solved MSS in " << niter << " iterations" << endl;

        if (niter < 0) {
            cout << "assert line 1518" << endl;
        }
        assert(niter >= 0);
    }
    catch(bad_alloc) {
        cout << "Threw exception 1 in solve!" << endl;
        exit(1);
    }
    catch(General_Exception excep) {
        cout << "Threw exception 2 in solve!" << endl;
        exit(1);
    }

    // Copy the output of the linear solve to the derivative vector
    // simutils::copyVals(3*pntList->size(), pcgRHS, qt->data());
    for (int i = 0; i < 3*pntList->size(); i++) {
        (*qt)(i) = pcgRHS[i];
    }

    // Do the Euler step to get the new node positions
    *q += dt*(*qt);

    E = ETemp;
    eta = etaTemp;
}

/**
 * Update according to the verlet method
*/
void MassSpring3D::verletSolve(double dt, int elementMode, bool initMode) {
    int pntId1;
    int pntId2;
    massPoint3D pnt1;
    massPoint3D pnt2;

    double ETemp = E;
    double etaTemp = eta;

    if (initMode) {
        E = 2.0;
        eta = 0.0;
    }

    // Compute the force vector for the current configuration
    for (auto edge = edgeList->begin(); edge != edgeList->end(); ++edge) {
        // Extract the points connecting the current spring
        pntId1 = edge->pntIds[0];
        pntId2 = edge->pntIds[1];

        pnt1 = pntList->at(pntId1);
        pnt2 = pntList->at(pntId2);

        // Calculate the contribution of the current point to the total force vector.
        if (elementMode == 0) {
            calcLocalElasticForce(*edge, pntId1, pnt1, pntId2, pnt2);
        } else {
            calcLocalKelvinForce(*edge, pntId1, pnt1, pntId2, pnt2);
        }
    }

    // Now, the force vector has been obtained, we perform an Euler step
    // Note: for the verlet solver, qt is not the velocity, it is the 
    double temp1, temp2, temp3;
    for (int i = 0; i < pntList->size(); i++) {
        if (iterCount == 0) {
            // Backup current position
            temp1 = (*q)[3*i];
            temp2 = (*q)[3*i+1];
            temp3 = (*q)[3*i+2];

            // Update the current position
            if (!(initMode && pntList->at(i).boundaryPnt))  {
                (*q)[3*i]   += dt*(*qt)[3*i]   + 0.5*simutils::square(dt)*((*f)[3*i]/pntList->at(i).mass);
                (*q)[3*i+1] += dt*(*qt)[3*i+1] + 0.5*simutils::square(dt)*((*f)[3*i+1]/pntList->at(i).mass);
                (*q)[3*i+2] += dt*(*qt)[3*i+2] + 0.5*simutils::square(dt)*((*f)[3*i+2]/pntList->at(i).mass);
            }

            // Update the solution history.
            (*qprev)[3*i]   = temp1;
            (*qprev)[3*i+1] = temp2;
            (*qprev)[3*i+2] = temp3;

            // Update the velocities
            (*qt)[3*i]   = ((*q)[3*i]   - (*qprev)[3*i])/dt;
            (*qt)[3*i+1] = ((*q)[3*i+1] - (*qprev)[3*i+1])/dt;
            (*qt)[3*i+2] = ((*q)[3*i+2] - (*qprev)[3*i+2])/dt;

        } else {
            temp1 = (*q)[3*i];
            temp2 = (*q)[3*i+1];
            temp3 = (*q)[3*i+2];

            // Update the current positions
            if (!(initMode && pntList->at(i).boundaryPnt)) {
                (*q)[3*i]   = 2.0*(*q)[3*i]   - (*qprev)[3*i]   + (simutils::square(dt)/pntList->at(i).mass)*(*f)[3*i];
                (*q)[3*i+1] = 2.0*(*q)[3*i+1] - (*qprev)[3*i+1] + (simutils::square(dt)/pntList->at(i).mass)*(*f)[3*i+1];
                (*q)[3*i+2] = 2.0*(*q)[3*i+2] - (*qprev)[3*i+2] + (simutils::square(dt)/pntList->at(i).mass)*(*f)[3*i+2];
            }

            // Update the velocities in each component
            (*qt)[3*i]   = ((*q)[3*i]   - (*qprev)[3*i]  )/(2.0*dt);
            (*qt)[3*i+1] = ((*q)[3*i+1] - (*qprev)[3*i+1])/(2.0*dt); 
            (*qt)[3*i+2] = ((*q)[3*i+2] - (*qprev)[3*i+2])/(2.0*dt); 

            // Update the solution history
            (*qprev)[3*i]   = temp1;
            (*qprev)[3*i+1] = temp2;
            (*qprev)[3*i+2] = temp3;
        }
    }

    E = ETemp;
    eta = etaTemp;
}

/**
 * Apply the boundary forces
*/
void MassSpring3D::applyBodyForces() {
    for (int id = 0; id < pntList->size(); id++) {
        (*f)[3*id] += gx;
        (*f)[3*id+1] += gy;
        (*f)[3*id+2] += gz;
    }
}

/**
 * Update the state of the mass-spring system using the information provided by the fluid
 * 
 * NOTE: stress does not have any ghost cells in the first axis, stress = [PoolU, PoolV]
 *       inherited from the pool class.
*/
void MassSpring3D::updateSolidVels(double dt, Pool3D &pool,
        double ****stress, double fNet[3], int ng, bool initMode) {
    
    this->initMode = initMode;
    int temp;

    if (!initMode && objType == SolidObject3D::ObjectType::STATIC) {
        return;
    }

    if (this->initMode) {
        temp = this->updateMode;
        this->updateMode = 0;
    }

    if (updateMode == 2) {
        this->dt = this->dtFix;
    } else {
        this->dt = dt;
    }

    // Initialize the force vector to 0.
    f->setZero();

    // Set the previous solution to current
    *qprev = *q;

    // Loop through all of the points. For each boundary point, add to the force vector.
    double colNet[3] = {0.0, 0.0, 0.0};
    if (!initMode) {
        applyBoundaryForces(pool, stress, ng, fNet, colNet);
        applyBodyForces();
    }

    // Apply the collision net force to the whole object
    for (int i = 0; i < pntList->size(); i++) {
        (*f)[3*i] += colNet[0];
        (*f)[3*i+1] += colNet[1];
        (*f)[3*i+2] += colNet[2];
    }

    // Loop through all of the edges, using the potential energy to compute the displacement of the
    // nodes.
    if (updateMode == 0) {
        // eulerSolve(dt, elementMode, initMode);
        verletSolve(dt, elementMode, initMode);
    } else if (updateMode == 1) {
        linearImplicitSolve(dt, elementMode, initMode);
    } else if (updateMode == 2) {
        admmSolver->step(500, 1e-10);
    }

    // Update the (x, y, z) position of each point in the point list using the updated q.
    for (int i = 0; i < pntList->size(); i++) {
        pntList->at(i).u = (*qt)[3*i];
        pntList->at(i).v = (*qt)[3*i+1];
        pntList->at(i).w = (*qt)[3*i+2];
    }

    iterCount++;
    if (!initMode) {
        nSteps++;
    }

    if (this->initMode) {
        this->updateMode = temp;
    }
}

/**
 * Make sure that the boundary points interpolate the SDF of the input Pool
*/
void MassSpring3D::interpBoundary(Pool3D &pool, bool resetRestingLengths) {

    double x, y, z, phi;
    double phiGrad[3];
    for (auto pnt = boundaryNodeIdList->begin(); pnt != boundaryNodeIdList->end(); ++pnt) {
        x = pntList->at(*pnt).x;
        y = pntList->at(*pnt).y;
        z = pntList->at(*pnt).z;

        pool.interpolatePhiGrad(x, y, z, phiGrad);
        simutils::normalize3D(phiGrad);

        phi = pool.interpolatePhi(x, y, z);

        pntList->at(*pnt).x -= phi*phiGrad[0];
        pntList->at(*pnt).y -= phi*phiGrad[1];
        pntList->at(*pnt).z -= phi*phiGrad[2];

        // Update the node locations.
        (*q)(3*(*pnt))       = pntList->at(*pnt).x;
        (*q)(3*(*pnt)+1)     = pntList->at(*pnt).y;
        (*q)(3*(*pnt)+2)     = pntList->at(*pnt).z;

        (*qprev)(3*(*pnt))   = pntList->at(*pnt).x;
        (*qprev)(3*(*pnt)+1) = pntList->at(*pnt).y;
        (*qprev)(3*(*pnt)+2) = pntList->at(*pnt).z;
    }

    // Set the resting lengths on the edges to what they are with the updated boundary positions
    if (resetRestingLengths) {
        double pnt1[3];
        double pnt2[3];
        for (auto edge = edgeList->begin(); edge != edgeList->end(); ++edge) {
            pnt1[0] = pntList->at(edge->pntIds[0]).x;
            pnt1[1] = pntList->at(edge->pntIds[0]).y;
            pnt1[2] = pntList->at(edge->pntIds[0]).z;

            pnt2[0] = pntList->at(edge->pntIds[1]).x;
            pnt2[1] = pntList->at(edge->pntIds[1]).y;
            pnt2[2] = pntList->at(edge->pntIds[1]).z;

            edge->l0 = simutils::eucDiff3D(pnt1, pnt2);
        }
    }
}


void MassSpring3D::updateSolidLocs(Pool3D &pool, bool interp) {
    for (int i = 0; i < pntList->size(); i++) {
        pntList->at(i).x = (*q)[3*i];
        pntList->at(i).y = (*q)[3*i+1];
        pntList->at(i).z = (*q)[3*i+2];
    }

    // Let the boundary points interpolate the pool.
    if (interp) {
        cout << "Do not have interpolation in 3D case yet" << endl;
        assert(false);
    }
}

/**
 * Output the nodes of the MSS with no explicit knowledge of the connectivity between them.
*/
void MassSpring3D::outputNodes(const char* fname) {
    this->outputSurfaceCentroids(fname);
}

/**
 * Output the node locations and velocities.
*/
void MassSpring3D::outputNodeVels(const char* fname) {

    ofstream outFile;
    outFile.open(fname);

    for (auto pnt = pntList->begin(); pnt != pntList->end(); ++pnt) {
        int id = pnt->nodeId;
        if (!pnt->boundaryPnt)
            continue;

        outFile << pnt->x << ", " << pnt->y << ", " << pnt->z << ", " << (*qt)(3*id)
            << ", " << (*qt)(3*id+1) << ", " << (*qt)(3*id+2) << endl;
    }

    outFile.close();
}

/**
 * Output the centroids of a triangle.
*/
void MassSpring3D::outputSurfaceCentroids(const char* fname) {

    ofstream outFile;
    outFile.open(fname);

    double centroid[3];

    for (auto face = faceList->begin(); face != faceList->end(); ++face) {
        this->computeCentroid(face->pntIds[0], face->pntIds[1], face->pntIds[2], centroid);

        outFile << centroid[0] << ", " << centroid[1] << ", " << centroid[2] << endl;
    }

    outFile.close();
}

/**
 * Output the pairs of nodes.
*/
void MassSpring3D::outputEdges(const char* fname) {
    ofstream outFile;
    outFile.open(fname);

    massPoint3D pnt1;
    massPoint3D pnt2;

    for (auto edge = edgeList->begin(); edge != edgeList->end(); ++edge) {
        if (!edge->boundaryEdge) {
            continue;
        }
        pnt1 = pntList->at(edge->pntIds[0]);
        pnt2 = pntList->at(edge->pntIds[1]);

        outFile << pnt1.x << ", " << pnt1.y << ", " << pnt1.z << endl;
        outFile << pnt2.x << ", " << pnt2.y << ", " << pnt2.z << endl;
    }

    outFile.close();
}

/**
 * Silly method to check for the number of redundant edges
*/
bool MassSpring3D::redundantEdges() {
    int numDup = 0;
    for (int i = 0; i < edgeList->size(); i++) {
        int edgeIPntId0 = edgeList->at(i).pntIds[0];
        int edgeIPntId1 = edgeList->at(i).pntIds[1];

        for (int j = 0; j < edgeList->size(); j++) {
            if (j == i)
                continue;

            int edgeJPntId0 = edgeList->at(j).pntIds[0];
            int edgeJPntId1 = edgeList->at(j).pntIds[1];

            if (((edgeIPntId0 == edgeJPntId0) && (edgeIPntId1 == edgeJPntId1))) {
                numDup += 1;
            }
        }
    }

    if (numDup <= 0) {
        cout << "There are " << numDup << " redundant edges in the MSS" << endl;
        assert(false);
    }

    return numDup > 0;
}

/**
 * Compute distance to input edge pntId1 -> pntId2
*/
double MassSpring3D::distToEdge(double x[3], int pntId1, int pntId2, double &t) {
    // Extract point coordinates
    double x1[3] = {pntList->at(pntId1).x, pntList->at(pntId1).y, pntList->at(pntId1).z};
    double x2[3] = {pntList->at(pntId2).x, pntList->at(pntId2).y, pntList->at(pntId2).z};

    // Edge vector
    double u[3] = {x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2]};

    // x - x1
    double w[3] = {x[0] - x1[0], x[1] - x1[1], x[2] - x1[2]};

    // Compute the vector projection from first edge point to line
    double coeff = simutils::ddot3d(u, w)/simutils::ddot3d(u, u);
    double proj[3] = {
        coeff*u[0],
        coeff*u[1],
        coeff*u[2]
    };

    // Compute the distance vector
    // proj - p
    double distVec[3];
    distVec[0] = proj[0] - w[0];
    distVec[1] = proj[1] - w[1];
    distVec[2] = proj[2] - w[2];

    // Compute the distance
    double d;
    if (simutils::sign(proj[0]) == simutils::sign(u[0])     &&
            simutils::sign(proj[1]) == simutils::sign(u[1]) &&
            simutils::sign(proj[2]) == simutils::sign(u[2])) {
        // Compute t
        t = simutils::eucNorm3D(proj) / simutils::eucNorm3D(u);
        if (t > 1.0) {
            t = 1.0;
            d = simutils::eucDiff3D(x, x2);
        } else {
            d = simutils::eucNorm3D(distVec);
        }
    } else {
        t = 0.0;
        d = simutils::eucDiff3D(x, x1);
    }

    return d;
}

/**
 * Find distance to closest point on MSS to input point
*/
double MassSpring3D::distToClosestBoundaryPoint(double x[3]) {
    double bPnt[3];
    double dMin = INFINITY;
    double d;
    for (auto pnt = boundaryNodeIdList->begin(); pnt != boundaryNodeIdList->end(); ++pnt) {
        bPnt[0] = pntList->at(*pnt).x;
        bPnt[1] = pntList->at(*pnt).y;
        bPnt[2] = pntList->at(*pnt).z;

        d = simutils::eucDiff3D(x, bPnt);

        if (d < dMin) {
            dMin = d;
        }
    }
    return dMin;
}

/**
 * Find distance to closest point on MSS to input point
*/
double MassSpring3D::distToClosestBoundaryEdge(double x[3]) {
    // double bPnt[3];
    double dMin = INFINITY;
    double d, t;
    for (auto edge = boundaryEdgeIdList->begin(); edge != boundaryEdgeIdList->end(); ++edge) {

        d = distToEdge(x, edgeList->at(*edge).pntIds[0], edgeList->at(*edge).pntIds[1], t);

        if (d < dMin) {
            dMin = d;
        }
    }
    return dMin;
}

/**
 * Compute the Barycentric coordinates of point x projected into the plane
 * defined by pnt1, pnt2, pnt3. Done as in Heindrich 2005.
 * 
*/
void MassSpring3D::projTriangle(double x[3], int pntId1, int pntId2,
                                int pntId3, double baryCoords[3]) {
    // Extract the vertices of the triangle, r1, r2, r3
    double r1[3] = {pntList->at(pntId1).x, pntList->at(pntId1).y, pntList->at(pntId1).z};

    double r2[3] = {pntList->at(pntId2).x, pntList->at(pntId2).y, pntList->at(pntId2).z};

    double r3[3] = {pntList->at(pntId3).x, pntList->at(pntId3).y, pntList->at(pntId3).z};

    // Convert into new variables to be inline with the formula I have written down.
    // u = r2 - r1
    double u[3] = {r2[0] - r1[0], r2[1] - r1[1], r2[2] - r1[2]};

    // v = r3 - r1
    double v[3] = {r3[0] - r1[0], r3[1] - r1[1], r3[2] - r1[2]};

    // w = x - r1
    double w[3] = { x[0] - r1[0],  x[1] - r1[1],  x[2] - r1[2]};

    // n = u x v
    double n[3];
    simutils::cross_product_3D(u, v, n);

    // u x w
    double uCrossW[3];
    simutils::cross_product_3D(u, w, uCrossW);

    // w x v
    double wCrossV[3];
    simutils::cross_product_3D(w, v, wCrossV);

    // Compute the barycentric coordinates for the projected point.
    double gamma = simutils::ddot3d(uCrossW, n) / simutils::ddot3d(n, n);
    double beta  = simutils::ddot3d(wCrossV, n) / simutils::ddot3d(n, n);
    double alpha = 1.0 - gamma - beta;

    baryCoords[0] = alpha;
    baryCoords[1] = beta;
    baryCoords[2] = gamma;
}


/**
 * Use Barycentric projection method Heindrich, 2005, Computing Barycentric
 * coordinates of a projected point. Output the Barycentric coordinates of
 * this point mapped onto the plane.
 * 
 * If the projected point is not within the triangle, we project it onto the nearest
 * edge.
 * 
 * Note: in the case where we instead use the closest edge, one of the coordinates
 * in the barycentric coordinates is 0.
 * 
 * baryCoords = [alph, beta, gam]
 * 
 * If ||baryCoords||_inf < 0 or ||baryCoords||_inf > 1, the point is not in
 * the triangle.
*/
double MassSpring3D::projTriangleDist(double x[3], int pntId1, int pntId2,
                                int pntId3, double baryCoords[3]) {
    
    projTriangle(x, pntId1, pntId2, pntId3, baryCoords);
    
    // Find the barycentric coordiantes of the projected point.
    // Extract the vertices of the triangle, r1, r2, r3
    double r1[3] = {pntList->at(pntId1).x, pntList->at(pntId1).y, pntList->at(pntId1).z};

    double r2[3] = {pntList->at(pntId2).x, pntList->at(pntId2).y, pntList->at(pntId2).z};

    double r3[3] = {pntList->at(pntId3).x, pntList->at(pntId3).y, pntList->at(pntId3).z};

    double alpha = baryCoords[0];
    double beta = baryCoords[1];
    double gamma = baryCoords[2];

    // if (false) {
    if (simutils::in_range(baryCoords[0], -EPS, 1.0+EPS)
            && simutils::in_range(baryCoords[1], -EPS, 1.0+EPS)
            && simutils::in_range(baryCoords[2], -EPS, 1.0+EPS)) {

        // Compute the coordinates of the projected point
        double projPnt[3] = {alpha*r1[0] + beta*r2[0] + gamma*r3[0],
                             alpha*r1[1] + beta*r2[1] + gamma*r3[1],
                             alpha*r1[2] + beta*r2[2] + gamma*r3[2]};
        return simutils::eucDiff3D(x, projPnt);
    } else {
        double t1, t2, t3;
        double d1 = distToEdge(x, pntId1, pntId2, t1);
        double d2 = distToEdge(x, pntId2, pntId3, t2);
        double d3 = distToEdge(x, pntId3, pntId1, t3);

        // Compute the barycentric coordinates based on the which distance is the least
        double dMin = simutils::dmin(d1, simutils::dmin(d2, d3)); 

        if (dMin == d1) {
            baryCoords[0] = 1.0 - t1;
            baryCoords[1] = t1;
            baryCoords[2] = 0.0;
        } else if (dMin == d2) {
            baryCoords[0] = 0.0;
            baryCoords[1] = 1.0 - t2;
            baryCoords[2] = t2;
        } else {
            baryCoords[0] = t3;
            baryCoords[1] = 0.0;
            baryCoords[2] = 1.0 - t3;
        }

        // Return least distance
        return dMin;
    }
}

/**
 * Compute the centroid of a triangle using three point ids.
*/
void MassSpring3D::computeCentroid(int pntId1, int pntId2, int pntId3, double centroid[3]) {
    centroid[0] = (pntList->at(pntId1).x + pntList->at(pntId2).x + pntList->at(pntId3).x)/3.0;
    centroid[1] = (pntList->at(pntId1).y + pntList->at(pntId2).y + pntList->at(pntId3).y)/3.0;
    centroid[2] = (pntList->at(pntId1).z + pntList->at(pntId2).z + pntList->at(pntId3).z)/3.0;
}

/**
 * Given an input point inPnt, project this to the nearest point on the MSS.
 * 
 * TODO: consider ways to speed up this algorithm.
*/
void MassSpring3D::interpFaceVels(double inPnt[3], double out[3]) {
    // Loop through all pairs of boundary edges along this point
    double baryCoords[3];
    int vertIds[3];

    double closestBaryCoords[3];
    int closestIds[3] = {-1, -1, -1};
    double closestDist = INFINITY;

    double d;

    for (auto face = faceList->begin(); face != faceList->end(); ++face) {
        vertIds[0] = face->pntIds[0];
        vertIds[1] = face->pntIds[1];
        vertIds[2] = face->pntIds[2];

        // Get the barycentric coordinates of the input point projected onto this face
        d = projTriangleDist(inPnt, vertIds[0], vertIds[1], vertIds[2], baryCoords);

        if (d < closestDist) {
            // Record the closest distance
            closestDist = d;

            // Record the point ids
            closestIds[0] = vertIds[0];
            closestIds[1] = vertIds[1];
            closestIds[2] = vertIds[2];

            // Record the barycentric coordaintes for performing interpolation.
            closestBaryCoords[0] = baryCoords[0];
            closestBaryCoords[1] = baryCoords[1];
            closestBaryCoords[2] = baryCoords[2];
        }
    }

    /* If we found a closest triangle, interpolate. Else interpolate on edges. */

    // Interpolate on this triangle.
    double alph = closestBaryCoords[0];
    double beta = closestBaryCoords[1];
    double gam  = closestBaryCoords[2];

    out[0] = alph*(*qt)[3*closestIds[0]]   + beta*(*qt)[3*closestIds[1]]   + gam*(*qt)[3*closestIds[2]];
    out[1] = alph*(*qt)[3*closestIds[0]+1] + beta*(*qt)[3*closestIds[1]+1] + gam*(*qt)[3*closestIds[2]+1];
    out[2] = alph*(*qt)[3*closestIds[0]+2] + beta*(*qt)[3*closestIds[1]+2] + gam*(*qt)[3*closestIds[2]+2];
}

/**
 * Find the interface point by searching along the "grid normal" direction until an interface is 
 * found.
 * 
 * Returns whether a point was found, mutates nearestInterfacePont
*/
bool MassSpring3D::findNearestGridNormalInterface(Pool3D &pool, int pntId, int nearestInterfacePnt[3]) {
    int nx = pool.getNx();
    int ny = pool.getNy();
    int nz = pool.getNz();

    // Compute the location on the grid of the input point
    int i = simutils::findLimInfMeshPoint(pntList->at(pntId).x, pool.getXMesh(), nx+1);
    int j = simutils::findLimInfMeshPoint(pntList->at(pntId).y, pool.getYMesh(), ny+1);
    int k = simutils::findLimInfMeshPoint(pntList->at(pntId).z, pool.getZMesh(), nz+1);

    // Compute the normal vector
    double n[3];
    pool.levelSetGradient(i, j, k, n);
    simutils::normalize3D(n);

    int nb[3];

    // Use a trick to compute the block normal vector

    // Compute projection onto x, y plane
    double s[3] = {n[0], n[1], n[2] - n[2]};
    double r = simutils::eucNorm3D(s);

    // x
    if (simutils::in_range(n[0], r*cos(5.0*M_PI/8.0), r*cos(3.0*M_PI/8.0))) {
        nb[0] = 0;
    } else if (n[0] >= 0.0) {
        nb[0] = 1;
    } else {
        nb[0] = -1;
    }

    // y
    if (simutils::in_range(n[1], r*sin(9.0*M_PI/8.0), r*cos(7.0*M_PI/8.0))) {
        nb[1] = 0;
    } else if (n[0] >= 0.0) {
        nb[1] = 1;
    } else {
        nb[1] = -1;
    }

    // z
    double rho = sqrt(simutils::square(simutils::eucNorm3D(s)) + simutils::square(n[2]));
    if (simutils::in_range(n[2], rho*cos(5.0*M_PI/8.0), rho*cos(3.0*M_PI/8.0))) {
        nb[2] = 0;
    } else if (n[2] >= 0) {
        nb[2] = 1;
    } else {
        nb[2] = -1;
    }


    // Search along this direction until the interface point is found
    objects::FSIObject obj = pool.objAtIndex(i, j, k);
    bool intFound = false;

    int so = 0;
    int count = 0;
    int jo, io, ko;
    while (!intFound && simutils::int_in_range(i+so*nb[0], 0, nx-1)
                            && simutils::int_in_range(j+so*nb[1], 0, ny-1)
                            && simutils::int_in_range(k+so*nb[2], 0, nz-1)) {
        io = i+so*nb[0];
        jo = j+so*nb[1];
        ko = k+so*nb[2];

        obj = pool.objAtIndex(io, jo, ko);
        intFound = pool.isInterface(obj);

        /* Do 6-connected band-search (there will be some overlap but its ok) */

        // east
        if (!intFound && simutils::int_in_range(io-1, 0, nx-1)
                            && simutils::int_in_range(jo, 0, ny-1)
                            && simutils::int_in_range(ko, 0, nz-1)) {
            obj = pool.objAtIndex(io-1, jo, ko);
            intFound = pool.isInterface(obj);
        }

        // west
        if (!intFound && simutils::int_in_range(io+1, 0, nx-1)
                            && simutils::int_in_range(jo, 0, ny-1)
                            && simutils::int_in_range(ko, 0, nz-1)) {
            obj = pool.objAtIndex(io+1, jo, ko);
            intFound = pool.isInterface(obj);
        }

        // south
        if (!intFound && simutils::int_in_range(io, 0, nx-1)
                            && simutils::int_in_range(jo-1, 0, ny-1)
                            && simutils::int_in_range(ko, 0, nz-1)) {
            obj = pool.objAtIndex(io, jo-1, ko);
            intFound = pool.isInterface(obj);
        }

        // north
        if (!intFound && simutils::int_in_range(io, 0, nx-1)
                            && simutils::int_in_range(jo+1, 0, ny-1)
                            && simutils::int_in_range(ko, 0, nz-1)) {
            obj = pool.objAtIndex(io, jo+1, ko);
            intFound = pool.isInterface(obj);
        }

        // down
        if (!intFound && simutils::int_in_range(io, 0, nx-1)
                            && simutils::int_in_range(jo, 0, ny-1)
                            && simutils::int_in_range(ko-1, 0, nz-1)) {
            obj = pool.objAtIndex(io, jo, ko-1);
            intFound = pool.isInterface(obj);
        }

        // up
        if (!intFound && simutils::int_in_range(io, 0, nx-1)
                            && simutils::int_in_range(jo, 0, ny-1)
                            && simutils::int_in_range(ko+1, 0, nz-1)) {
            obj = pool.objAtIndex(io, jo, ko+1);
            intFound = pool.isInterface(obj);
        }

        // Top corners
        if (!intFound && simutils::int_in_range(io+1, 0, nx-1)
                            && simutils::int_in_range(jo+1, 0, ny-1)
                            && simutils::int_in_range(ko+1, 0, nz-1)) {
            obj = pool.objAtIndex(io+1, jo+1, ko+1);
            intFound = pool.isInterface(obj);
        }
        if (!intFound && simutils::int_in_range(io-1, 0, nx-1)
                            && simutils::int_in_range(jo+1, 0, ny-1)
                            && simutils::int_in_range(ko+1, 0, nz-1)) {
            obj = pool.objAtIndex(io-1, jo+1, ko+1);
            intFound = pool.isInterface(obj);
        }
        if (!intFound && simutils::int_in_range(io-1, 0, nx-1)
                            && simutils::int_in_range(jo-1, 0, ny-1)
                            && simutils::int_in_range(ko+1, 0, nz-1)) {
            obj = pool.objAtIndex(io-1, jo-1, ko+1);
            intFound = pool.isInterface(obj);
        }
        if (!intFound && simutils::int_in_range(io+1, 0, nx-1)
                            && simutils::int_in_range(jo-1, 0, ny-1)
                            && simutils::int_in_range(ko+1, 0, nz-1)) {
            obj = pool.objAtIndex(io+1, jo-1, ko+1);
            intFound = pool.isInterface(obj);
        }

        // Mid corners
        if (!intFound && simutils::int_in_range(io+1, 0, nx-1)
                            && simutils::int_in_range(jo+1, 0, ny-1)
                            && simutils::int_in_range(ko, 0, nz-1)) {
            obj = pool.objAtIndex(io+1, jo+1, ko);
            intFound = pool.isInterface(obj);
        }
        if (!intFound && simutils::int_in_range(io-1, 0, nx-1)
                            && simutils::int_in_range(jo+1, 0, ny-1)
                            && simutils::int_in_range(ko, 0, nz-1)) {
            obj = pool.objAtIndex(io-1, jo+1, ko);
            intFound = pool.isInterface(obj);
        }
        if (!intFound && simutils::int_in_range(io-1, 0, nx-1)
                            && simutils::int_in_range(jo-1, 0, ny-1)
                            && simutils::int_in_range(ko, 0, nz-1)) {
            obj = pool.objAtIndex(io-1, jo-1, ko);
            intFound = pool.isInterface(obj);
        }
        if (!intFound && simutils::int_in_range(io+1, 0, nx-1)
                            && simutils::int_in_range(jo-1, 0, ny-1)
                            && simutils::int_in_range(ko, 0, nz-1)) {
            obj = pool.objAtIndex(io+1, jo-1, ko);
            intFound = pool.isInterface(obj);
        }

        // Mid corners
        if (!intFound && simutils::int_in_range(io+1, 0, nx-1)
                            && simutils::int_in_range(jo+1, 0, ny-1)
                            && simutils::int_in_range(ko-1, 0, nz-1)) {
            obj = pool.objAtIndex(io+1, jo+1, ko-1);
            intFound = pool.isInterface(obj);
        }
        if (!intFound && simutils::int_in_range(io-1, 0, nx-1)
                            && simutils::int_in_range(jo+1, 0, ny-1)
                            && simutils::int_in_range(ko-1, 0, nz-1)) {
            obj = pool.objAtIndex(io-1, jo+1, ko-1);
            intFound = pool.isInterface(obj);
        }
        if (!intFound && simutils::int_in_range(io-1, 0, nx-1)
                            && simutils::int_in_range(jo-1, 0, ny-1)
                            && simutils::int_in_range(ko-1, 0, nz-1)) {
            obj = pool.objAtIndex(io-1, jo-1, ko-1);
            intFound = pool.isInterface(obj);
        }
        if (!intFound && simutils::int_in_range(io+1, 0, nx-1)
                            && simutils::int_in_range(jo-1, 0, ny-1)
                            && simutils::int_in_range(ko-1, 0, nz-1)) {
            obj = pool.objAtIndex(io+1, jo-1, ko-1);
            intFound = pool.isInterface(obj);
        }



        if (!intFound) {
            count++;
            so += ( (count % 2) ? 1 : -1 )*count;
        }
    }

    if (intFound) {
        nearestInterfacePnt[0] = i+so*nb[0];
        nearestInterfacePnt[1] = j+so*nb[1];
        nearestInterfacePnt[2] = k+so*nb[2];
    } else {
        nearestInterfacePnt[0] = 0;
        nearestInterfacePnt[1] = 0;
        nearestInterfacePnt[2] = 0;

        // This case does not make sense. Terminate the program execution
        cout << "Could not find stess to extrapolate on MSS" << endl;
        assert(false);
    }

    return intFound;
}

/**
 * Add a triple to the face map. The important part is that the triple (a, b, c)
 * has a < b < c
*/
tuple<int, int, int> MassSpring3D::createOrderedTuple(int a, int b, int c) {
    int temp[3] = {a, b, c};

    int t;

    //  <--> first swap
    // [a, b, c]
    if (temp[1] < temp[0]) {
        t = temp[1];
        temp[1] = temp[0];
        temp[0] = t;
    } 

    //     <--> 2nd swap
    // [a, b, c]
    if (temp[2] < temp[1]) {
        t = temp[2];
        temp[2] = temp[1];
        temp[1] = t;
    }

    //  <--> 3rd swap
    // [a, b, c]
    if (temp[1] < temp[0]) {
        t = temp[1];
        temp[1] = temp[0];
        temp[0] = t;
    } 

    // Make sure everything is in order.
    assert((temp[0] < temp[1]) && (temp[1] < temp[2]) && (temp[0] < temp[2]));

    return make_tuple(temp[0], temp[1], temp[2]);
}

/**
 * Find closest distance to the boundary faces of the MSS
*/
// double MassSpring3D::closestBoundaryDist(double inPnt[3], KDTree3D *tree, KDPointCloud3D *pntCloud) {
//     double baryCoords[3];
//     double closestDist = INFINITY;
//     double d;

//     const int MAXCANDS = 10;
//     vector<pair<size_t,double> > ret_matches;
//     std::vector<size_t> ret_index(MAXCANDS);
//     std::vector<double> out_dist_sqr(MAXCANDS);

//     // Candidate list
//     int numFound = tree->knnSearch(&inPnt[0],
//                     MAXCANDS, &ret_index[0], &out_dist_sqr[0]);
    
//     assert(numFound > 0);
    
//     for (int i = 0; i < numFound; i++) {
//         mass_spring::massPoint3D pnt = *pntCloud->points->at(ret_index.at(i));

//         // Must correspond to this one
//         if (pnt.structNum != this->structNum) {
//             continue;
//         }

//         // Iterate through faces
//         for (auto face = pnt.faceIds.begin(); face != pnt.faceIds.end(); ++face) {
//             face3D faceLoc = faceList->at(*face);
//             d = projTriangleDist(inPnt, faceLoc.pntIds[0], faceLoc.pntIds[1], faceLoc.pntIds[2], baryCoords);

//             if (d < closestDist) {
//                 // Record the closest distance
//                 closestDist = d;
//             }

//         }
//     }


//     // for (auto face = faceList->begin(); face != faceList->end(); ++face) {
//     //     // Get the distance of the input point to this triangle (baryCoords are irrelevant)
//     //     d = projTriangleDist(inPnt, face->pntIds[0], face->pntIds[1], face->pntIds[2], baryCoords);

//     //     if (d < closestDist) {
//     //         // Record the closest distance
//     //         closestDist = d;
//     //     }
//     // }
//     return closestDist;
// }

/**
 * Find closest distance to the boundary faces of the MSS
*/
// void MassSpring3D::closestBoundaryPnt(double inPnt[3], double outPnt[3], KDTree3D *tree, KDPointCloud3D *pntCloud) {
//     double baryCoords[3];
//     double closestDist = INFINITY;
//     double d;

//     int closestId = -1;

//     double baryCoordsClosest[3];
//     const int MAXCANDS = 10;
//     vector<pair<size_t,double> > ret_matches;
//     std::vector<size_t> ret_index(MAXCANDS);
//     std::vector<double> out_dist_sqr(MAXCANDS);

//     // Candidate list
//     int numFound = tree->knnSearch(&inPnt[0],
//                     MAXCANDS, &ret_index[0], &out_dist_sqr[0]);
    
//     assert(numFound > 0);
    
//     for (int i = 0; i < numFound; i++) {
//         mass_spring::massPoint3D pnt = *pntCloud->points->at(ret_index.at(i));

//         // Must correspond to this one
//         if (pnt.structNum != this->structNum) {
//             continue;
//         }

//         // Iterate through faces
//         for (auto face = pnt.faceIds.begin(); face != pnt.faceIds.end(); ++face) {
//             face3D faceLoc = faceList->at(*face);
//             d = projTriangleDist(inPnt, faceLoc.pntIds[0], faceLoc.pntIds[1], faceLoc.pntIds[2], baryCoords);

//             if (d < closestDist) {
//                 // Record the closest distance
//                 closestDist = d;
//                 closestId = *face;

//                 simutils::copyVals(3, baryCoords, baryCoordsClosest);
//             }

//         }
//     }

//     // for (auto face = faceList->begin(); face != faceList->end(); ++face) {
//     //     // Get the distance of the input point to this triangle (baryCoords are irrelevant)
//     //     d = projTriangleDist(inPnt, face->pntIds[0], face->pntIds[1], face->pntIds[2], baryCoords);

//     //     if (d < closestDist) {
//     //         // Record the closest distance
//     //         closestDist = d;

//     //         // Record the nearest IDs and coordinates
//     //         closestId = face - faceList->begin();

//     //         simutils::copyVals(3, baryCoords, baryCoordsClosest);
//     //     }
//     // }

//     // Make the output point the Barycentric reconstruction of the closest point recorded
//     for (int j = 0; j < 3; j++) {
//         outPnt[j] = 0.0;
//     }

//     double pnt[3];
//     for (int i = 0; i < 3; i++) {
//         pnt[0] = pntList->at(faceList->at(closestId).pntIds[i]).x;
//         pnt[1] = pntList->at(faceList->at(closestId).pntIds[i]).y;
//         pnt[2] = pntList->at(faceList->at(closestId).pntIds[i]).z;

//         for (int j = 0; j < 3; j++) {
//             outPnt[j] += baryCoordsClosest[i]*pnt[j];
//         }
//     }
// }

/**
 * Create the list of faces for the MSS.
 * 
 * First, we loop through the edges and remove intersecting springs.
 * 
 * Secondly, we use a brute-force algorithm to find all faces.
 * 
 * TODO: I have no clue why this works in the cases that it does. Using triangle intersection may help.
 *       My notes on this are in page 9 of Fund Var Mesh Adapt
*/
void MassSpring3D::createFaceList() {

    set<tuple<int, int, int>> faceSet;

    // Flag vector to indicate edges to be popped.
    int edgeLen = edgeList->size();
    bool *toBeRemoved = new bool[edgeLen];
    for (int i = 0; i < edgeLen; i++) {
        toBeRemoved[i] = false;
    }

    // cout << "size of pntList = " << pntList->size() << endl;

    /* Part 1: brute-force algorithm to create a set of unique faces */
    // This is a freaking mess lol, but it works
    for (int pntOneIdx = 0; pntOneIdx < pntList->size(); pntOneIdx++) {
        massPoint3D pntOne = pntList->at(pntOneIdx);
        for (auto edgeOneId = pntOne.edgeIds.begin(); edgeOneId != pntOne.edgeIds.end(); ++edgeOneId) {

            // If this is a removed edge, continue
            if (toBeRemoved[*edgeOneId] || !edgeList->at(*edgeOneId).boundaryEdge)
                continue;

            // ID of the second point
            int pntTwoIdx = edgeList->at(*edgeOneId).pntIds[0];
            pntTwoIdx = (pntTwoIdx == pntOneIdx) ? edgeList->at(*edgeOneId).pntIds[1] : pntTwoIdx;
            massPoint3D pntTwo = pntList->at(pntTwoIdx);
 
            // Iterate through all edges of the second point, adding all 3-cycles.
            for (auto edgeTwoId = pntTwo.edgeIds.begin(); edgeTwoId != pntTwo.edgeIds.end(); ++edgeTwoId) {

                // If this is a removed edge, or if it's the same edge as the first, continue
                if (toBeRemoved[*edgeTwoId] || *edgeTwoId == *edgeOneId || !edgeList->at(*edgeTwoId).boundaryEdge)
                    continue;

                // Check that this edge connects to the first point.
                int pntThreeIdx = edgeList->at(*edgeTwoId).pntIds[0];
                pntThreeIdx = (pntThreeIdx == pntTwoIdx) ? edgeList->at(*edgeTwoId).pntIds[1] : pntThreeIdx;
                massPoint3D pntThree = pntList->at(pntThreeIdx);

                for (auto edgeThreeId = pntThree.edgeIds.begin(); edgeThreeId != pntThree.edgeIds.end(); ++edgeThreeId) {
                    if (toBeRemoved[*edgeThreeId] || *edgeThreeId == *edgeOneId || *edgeThreeId == *edgeTwoId || !edgeList->at(*edgeThreeId).boundaryEdge) 
                        continue;

                    if (edgeList->at(*edgeThreeId).pntIds[0] == pntOneIdx || edgeList->at(*edgeThreeId).pntIds[1] == pntOneIdx) {
                        faceSet.insert(this->createOrderedTuple(pntOneIdx, pntTwoIdx, pntThreeIdx));
                    }
                }
            }
        }
    }

    // cout << "entering face assign loop" << endl;
    int faceId = 0;
    for (auto tup = faceSet.begin(); tup != faceSet.end(); ++tup) {
        faceList->push_back(face3D());
        faceList->at(faceId).pntIds[0] = get<0>(*tup);
        faceList->at(faceId).pntIds[1] = get<1>(*tup);
        faceList->at(faceId).pntIds[2] = get<2>(*tup);

        faceId++;
    }
    // cout << "FINSIHED face assign loop" << endl;

    // Use the Barycentric coordinates of the centroid to detect which faces are intersecting.
    // Remove those who have a barycentric coordinate too close to 0.
    bool finished = false;
    double centroidOne[3];
    double centroidTwo[3];
    double bCoor1[3];
    double bCoor2[3];
    vector<int> removedIdx;
    int curIdx = 0;
    while (!finished) {
        // cout << "comp roid" << endl;
        // cout << "size of faceList = " << faceList->size() << endl;
        face3D faceOne = faceList->at(curIdx);
        computeCentroid(faceOne.pntIds[0], faceOne.pntIds[1], faceOne.pntIds[2], centroidOne);
        // cout << "finsihed comp roid" << endl;

        for (int i = curIdx+1; i < faceList->size()-1; i++) {
            // Compute the centroid of the second face
            face3D faceTwo = faceList->at(i);
            computeCentroid(faceTwo.pntIds[0], faceTwo.pntIds[1], faceTwo.pntIds[2], centroidTwo);

            int num_matches = 0;
            for (int n = 0; n < 3; n++) {
                for (int m = 0; m < 3; m++) {
                    if (faceOne.pntIds[n] == faceTwo.pntIds[m]) {
                        num_matches++;
                        break;
                    }
                }
            }

            if (num_matches >= 2) {
                // Check if their centroids do not belong inside of each other
                projTriangle(centroidTwo, faceOne.pntIds[0],
                                faceOne.pntIds[1], faceOne.pntIds[2],
                                bCoor1);
                projTriangle(centroidOne, faceTwo.pntIds[0],
                                faceTwo.pntIds[1], faceTwo.pntIds[2],
                                bCoor2);

                if ((bCoor1[0] >= -EPS && bCoor1[1] >= -EPS && bCoor1[2] >= -EPS) ||
                        (bCoor2[0] >= -EPS && bCoor2[1] >= -EPS && bCoor2[2] >= -EPS)) {
                    removedIdx.push_back(i);
                }
            }
        }

        // Remove all of the elements from faceList which overlap
        // cout << "removing elements" << endl;
        int len = removedIdx.size();
        int nRem = 0;
        for (int i = 0; i < len; i++) {
            faceList->erase(faceList->begin()+removedIdx.at(i)-nRem);
            nRem++;
        }

        // Clear the vector
        removedIdx.clear();

        // Increment the index we are looking at
        curIdx++;

        // Update finished, based on whether we have looked at all of the points
        finished = (curIdx == faceList->size()-1);
    }
    
    // For each point in the face, indicate its membership
    // cout << "indicating membership" << endl;
    for (int i = 0; i < faceList->size(); i++) {
        pntList->at((faceList->at(i)).pntIds[0]).faceIds.push_back(i);
        pntList->at((faceList->at(i)).pntIds[1]).faceIds.push_back(i);
        pntList->at((faceList->at(i)).pntIds[2]).faceIds.push_back(i);
    }
    // cout << "FINISHED indicating membership" << endl;

    delete[] toBeRemoved;
}

/**
 * Get the number of boundary neighbours of a boundary mass point
*/
int MassSpring3D::getNumNeighs(mass_spring::massPoint3D pnt) {
    int count = 0;

    for (auto edge = pnt.edgeIds.begin(); edge != pnt.edgeIds.end(); ++edge) {
        if (edgeList->at(*edge).boundaryEdge) {
            count++;
        }
    }

    return count;
}

/**
 * Detect whether the input point is along the input edge
*/
bool MassSpring3D::pntAlongEdge(double pnt[3], double X0[3], double X1[3]) {

    // If the input point is equal to either of the end points, return true
    if (simutils::eps_equal(simutils::eucDiff3D(pnt, X0), 0.0, EPS)
        || simutils::eps_equal(simutils::eucDiff3D(pnt, X1), 0.0, EPS)) {

        return true;
    }

    // Vector to keep track of which coordinates are used to compute t. The ones
    // that are not chosen are those where X0[i] = X1[i], must check that
    // pnt[i] = X0[i]
    bool toBeChecked[3] = {
        simutils::eps_equal(X0[0], X1[0], EPS),
        simutils::eps_equal(X0[1], X1[1], EPS),
        simutils::eps_equal(X0[2], X1[2], EPS)
    };

    // Input point must have same value on constant edges.
    for (int i = 0; i < 3; i++) {
        if (toBeChecked[i]) {
            if (!simutils::eps_equal(X0[i], pnt[i], EPS)) {
                return false;
            }
        }
    }

    // Compute an estimate for t.
    double t;
    int t_ind = -1;
    for (int i = 0; i < 3; i++) {
        if (!toBeChecked[i]) {
            t = (pnt[i] - X0[i]) / (X1[i] - X0[i]);
            t_ind = i;
            break;
        }
    }

    // Check that it is consistant in the edges to be checked. If there is any
    // inconsistancy, return false.
    double t_new;
    for (int i = 0; i < 3; i++) {
        if (!toBeChecked[i] && i != t_ind) {
            t_new = (pnt[i] - X0[i]) / (X1[i] - X0[i]);
            // If there is an inconsistancy, return false
            if (!simutils::eps_equal(t, t_new, EPS)) {
                return false;
            }
        }
    }

    return true;
}

double MassSpring3D::pointDiff(massPoint3D pnt1, massPoint3D pnt2) {
    return sqrt(simutils::square(pnt1.x - pnt2.x) + simutils::square(pnt1.y - pnt2.y)
                    + simutils::square(pnt1.z - pnt2.z));
}

MassSpring3D::~MassSpring3D() {
    // cout << "In destructor MSS" << endl;
    // cout << "edgeList" << endl;
    delete edgeList;
    // cout << "pntList" << endl;
    delete pntList;
    // cout << "faceList" << endl;
    delete faceList;
    // cout << "bedgeIdList" << endl;
    delete boundaryEdgeIdList;
    // cout << "bNodeIdList" << endl;
    delete boundaryNodeIdList;
    // cout << "finsihed deleting boundary stuff" << endl;

    // cout << "deleting newtons homies" << endl;
    delete f;
    delete q;
    delete qt;
    delete qprev;
    delete qBackup;
    delete qtBackup;
    delete qprevBackup;
    delete fBackup;
    // cout << "finished deleting newtons homies" << endl;

    if (updateMode == 1) {
        // cout << "deleting CG" << endl;
        delete[] pcgRHS;
        delete[] pcgTol;
        delete matrix;
        delete params;
        // cout << "FINSIHED deleting CG" << endl;
    } else if (updateMode == 2) {
        delete admmSolver;
    }
    // cout << "FINISHED In destructor" << endl;
}