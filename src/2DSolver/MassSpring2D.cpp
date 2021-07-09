#include "MassSpring2D.h"
#include "Pool2D.h"
#include "SolidObject.h"
#include <vector>
#include <assert.h>
#include <fstream>
#include <math.h>
#include "../Utils/SimUtilities.h"
#include <cmath>
#include "../LASolver/SparseItObj.h"
#include <Eigen/Dense>

using namespace mass_spring;
using namespace std;

/**
 * Just the copy ctor
*/
MassSpring2D::MassSpring2D(const MassSpring2D &cpy) : Assembly() {
    // cout << "In MSS copy ctor" << endl;
    this->pntList = new vector<massPoint2D>(*cpy.pntList);
    this->edgeList = new vector<edge2D>(*cpy.edgeList);
    this->boundaryEdgeIdList = new vector<int>(*cpy.boundaryEdgeIdList);
    this->boundaryNodeIdList = new vector<int>(*cpy.boundaryNodeIdList);
    this->w = cpy.w;
    this->admmTol = cpy.admmTol;

    this->objType = cpy.objType;

    /* The collision detection points */
    nodeCols = new vector<set<massPoint2D*>>(*cpy.nodeCols);

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
    qtBackup = new Eigen::VectorXd(*(cpy.qtBackup));
    qprevBackup = new Eigen::VectorXd(*(cpy.qprevBackup));
    f = new Eigen::VectorXd(*(cpy.f));

    this->eta = cpy.eta;
    this->E = cpy.E;
    this->iterCount = cpy.iterCount;
    this->mass = cpy.mass;
    this->density = cpy.density;
    this->bb_jl = cpy.bb_jl;
    this->bb_jr = cpy.bb_jr;
    this->bb_il = cpy.bb_il;
    this->bb_ir = cpy.bb_ir;
    this->pntMass = cpy.pntMass;

    // Create the matrix (matrixStruc objects do not have a copy ctor)
    if (this->updateMode == 1) {
        this->matrix = createStiffnessMatrix();
        this->params = setUpParams();
        this->nnz = cpy.nnz;
        this->pcgTol = new double[2*pntList->size()];
        // this->edgeChecked = new bool[edgeList->size()];
        this->pcgRHS = new double[2*pntList->size()];
        for (int i = 0; i < 2*pntList->size(); i++) {
            this->pcgTol[i] = cpy.pcgTol[i];
            this->pcgRHS[i] = cpy.pcgRHS[i];
        }
    } else if (this->updateMode == 2) {
        this->D = new Eigen::SparseMatrix<double>(*cpy.D);
        this->W = new Eigen::SparseMatrix<double>(*cpy.W);
        this->M = new Eigen::SparseMatrix<double>(*cpy.M);

        setNPnts(pntList->size());
        setD(2);
        this->admmSolver = new ADMMPG((cpy.admmSolver)->dt, *this);
    }
}

/**
 * mode = 0: explicit time stepping
 * mode = 1: semi-implicit
 * mode = 2: optimization
*/
MassSpring2D::MassSpring2D(Pool2D &pool, int structNum,
        SolidObject &obj, int updateMode, int elementMode) : Assembly() {
    // Create the vectors for the edge and point lists
    pntList = new vector<massPoint2D>();
    edgeList = new vector<edge2D>();
    boundaryEdgeIdList = new vector<int>();
    boundaryNodeIdList = new vector<int>();

    this->updateMode = updateMode;
    this->elementMode = elementMode;

    // Iterate through the Pool, adding to the vector struct as we go.
    int domain;
    objects::FSIObject poolObj;
    objType = obj.getObjType();

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

    double x, y;

    // bounding box indices
    bb_jl = 0; bb_jr = 0; bb_il = 0; bb_ir = 0;

    // Map for the object
    map<pair<int, int>, int> objectPointLabels;

    // First loop. Build up point set and HashMap for point numbers
    for (int j = 0; j < pool.getNy(); j++) {
        for (int i = 0; i < pool.getNx(); i++) {
            domain = pool.domainMembership(i, j);

            poolObj = pool.objAtIndex(i, j);

            // if (domain == structNum && poolObj == objects::STRUCTURE) {
            if (domain == structNum && poolObj != objects::FLUID_C) {
                // x, y coordinate of the current point. Interpolate in the normal direction to
                // find the interface if this is a boundary cell
                if (pool.isInterface(poolObj)) {
                    double pnt[2];
                    this->interpolateBoundaryLocation(pool, i, j, pnt);

                    int nDir[2];
                    pool.getNormalDir(pool.objAtIndex(i, j), nDir);

                    // Negate to find the inward normal
                    nDir[0] *= -1;
                    nDir[1] *= -1;

                    // Coordinates of the structure point in the neighbouring direction
                    double structPnt[2] = {simutils::midpoint(pool.getXMeshVal(i+nDir[0]), pool.getXMeshVal(i+nDir[0]+1)),
                                        simutils::midpoint(pool.getYMeshVal(j+nDir[1]), pool.getYMeshVal(j+nDir[1]+1))};
                    
                    // // Compute the unit outward normal
                    double outN[2] = {-((double)nDir[0]), -((double)nDir[1])};
                    simutils::normalize2D(outN);

                    double diff[2] = {pnt[0] - structPnt[0], pnt[1] - structPnt[1]};

                    if (simutils::eps_equal(simutils::eucNorm2D(diff), 0.0, 1e-16)) {
                        // Add a small purtabation to the outward normal. Should not contribute signficantly to
                        // the overall error.
                        const double ADJ = 1e-16;

                        x = pnt[0] + ADJ*outN[0];
                        y = pnt[1] + ADJ*outN[1];
                    } else {
                        x = pnt[0];
                        y = pnt[1];
                    }
                } else {
                    x = simutils::midpoint(pool.getXMeshVal(i), pool.getXMeshVal(i+1));
                    y = simutils::midpoint(pool.getYMeshVal(j), pool.getYMeshVal(j+1));
                }

                // Add the current point to the point list
                pntList->push_back(massPoint2D());
                pntList->at(pntId).x = x;
                pntList->at(pntId).y = y;
                pntList->at(pntId).u = obj.getU0();
                pntList->at(pntId).v = obj.getV0();
                pntList->at(pntId).structNum = structNum;
                pntList->at(pntId).nodeId = pntId;

                // If this is a boundary point, it should be indicated
                // TODO: when we add proper boundary points, we should reconsider this.
                // pntList->at(pntId).boundaryPnt = pool.oneGridFromInterface(i, j);
                pntList->at(pntId).boundaryPnt = pool.isInterface(poolObj);

                // If it is a boundary point, add to container of boundary lists for faster looping
                if (pntList->at(pntId).boundaryPnt) {
                    boundaryNodeIdList->push_back(pntId);
                }

                // Add the current label to the map
                objectPointLabels[make_pair(i, j)] = pntId;

                // Update the bounding box
                bb_jl = simutils::imin(bb_jl, j);
                bb_il = simutils::imin(bb_il, i);

                bb_jr = simutils::imax(bb_jr, j);
                bb_ir = simutils::imax(bb_ir, i);
                
                pntId++;
            }
        }
    }

    // Compute and then assign the point masses. Additionally, create the q vector.
    this->pntMass = obj.getMass()/pntList->size();

    q     = new Eigen::VectorXd(Eigen::VectorXd::Constant(2*pntList->size(), 0.0));
    qprev = new Eigen::VectorXd(Eigen::VectorXd::Constant(2*pntList->size(), 0.0));
    qt    = new Eigen::VectorXd(Eigen::VectorXd::Constant(2*pntList->size(), 0.0));
    qBackup     = new Eigen::VectorXd(Eigen::VectorXd::Constant(2*pntList->size(), 0.0));
    qprevBackup = new Eigen::VectorXd(Eigen::VectorXd::Constant(2*pntList->size(), 0.0));
    qtBackup    = new Eigen::VectorXd(Eigen::VectorXd::Constant(2*pntList->size(), 0.0));
    f     = new Eigen::VectorXd(Eigen::VectorXd::Constant(2*pntList->size(), 0.0));

    int qoff = 0;
    for (int i = 0; i < pntList->size(); i++) {
        pntList->at(i).mass = pntMass;

        // Add coordinates to the q vector.
        (*q)[qoff] = pntList->at(i).x;
        (*q)[qoff+1] = pntList->at(i).y;
        (*qBackup)[qoff] = pntList->at(i).x;
        (*qBackup)[qoff+1] = pntList->at(i).y;

        (*qt)[qoff] = pntList->at(i).u;
        (*qt)[qoff+1] = pntList->at(i).v;
        (*qtBackup)[qoff] = pntList->at(i).u;
        (*qtBackup)[qoff+1] = pntList->at(i).v;

        qoff += 2;
    }

    // Second loop, go through the bounding box and build the edge connections
    // Note: we take advantage of the fact that we start from the bottom left,
    //       and move toward the top right.
    edgeId = 0;
    int neighId;
 
    double pntLoc[2];
    double diff[2];

    double hx = pool.getXMeshVal(1) -  pool.getXMeshVal(0);
    double hy = pool.getYMeshVal(1) -  pool.getYMeshVal(0);

    for (int j = bb_jl; j <= bb_jr; j++) {
        for (int i = bb_il; i <= bb_ir; i++) {
            domain = pool.domainMembership(i, j);
            poolObj = pool.objAtIndex(i, j);

            // if (domain == structNum && poolObj == objects::STRUCTURE) {
            if (domain == structNum && poolObj != objects::FLUID_C) {
                // Get the ID of the current point.
                pntId = objectPointLabels[make_pair(i, j)];

                // Coordinates of the current poi
                pntLoc[0] = pntList->at(pntId).x;
                pntLoc[1] = pntList->at(pntId).y;

                // Create arrays for i, j values that must be checked
                int niList[4] = {i+1, i-1, i, i+1};
                int njList[4] = {j, j+1, j+1, j+1};

                for (int k = 0; k < 4; k++) {
                    int ni = niList[k];
                    int nj = njList[k];

                    // if (pool.objAtIndex(ni, nj) == objects::STRUCTURE) {
                    if (pool.objAtIndex(ni, nj) != objects::FLUID_C) {
                        // Extract some information about the neighbour in this direction.
                        neighId = objectPointLabels[make_pair(ni, nj)];

                        // Create the new edge
                        edgeList->push_back(edge2D());
                        edgeId = edgeList->size() - 1;

                        // Compute the relaxed length of this spring, taking the norm of the difference in their coordinates
                        diff[0] = pntLoc[0] - pntList->at(neighId).x;
                        diff[1] = pntLoc[1] - pntList->at(neighId).y;

                        double l0;
                        l0 = simutils::eucNorm2D(diff);

                        assert(l0 != 0.0);
                        
                        edgeList->at(edgeId).l0 = l0;
                        edgeList->at(edgeId).pntIds[0] = pntId;
                        edgeList->at(edgeId).pntIds[1] = neighId;
                        edgeList->at(edgeId).boundaryEdge = pntList->at(pntId).boundaryPnt
                                                            && pntList->at(neighId).boundaryPnt;
                        
                        // If this is a boundary edge, set the ideal edge length to 1.2*hdiag
                        double fScale = 0.8;

                        if (!edgeList->at(edgeId).boundaryEdge) {
                            // If this is not a boundary edge, set its target length to fScale*|| [hx \\ hy] ||
                            edgeList->at(edgeId).l0 = fScale*sqrt(simutils::square(hx) + simutils::square(hy));
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

    // Now, attempt to update the solid to a steady state
    double eps = 1e-10;
    const int MAX_ITERS = 1000;
    double dt = 0.1*simutils::dmin(hx, hy);
    int iters = 0;

    double fNet[2] = {0.0, 0.0};
    do  {
        this->updateSolidVels(dt, pool, NULL, fNet, 1, true);
        this->updateSolidLocs(pool, false);
        iters ++;
    }
    while (qt->lpNorm<1>() > eps &&  iters < MAX_ITERS);

    this->iterCount = 0;

    // After this is done, set the edge length of each spring to its current config to bring the
    // potential energy to 0.
    double barVec[2];
    for (auto edge = edgeList->begin(); edge != edgeList->end(); ++edge) {
        // Vector connecting this edge.
        barVec[0] = pntList->at(edge->pntIds[1]).x - pntList->at(edge->pntIds[0]).x;
        barVec[1] = pntList->at(edge->pntIds[1]).y - pntList->at(edge->pntIds[0]).y;

        // Set this bars resting length to this value
        edge->l0 = simutils::eucNorm2D(barVec);
    }

    // Set the derivatives of the locations to what they should be for the uniform initial speed
    // of the solid (positions should be correct).
    for (int i = 0; i < pntList->size(); i++) {
        (*q)(2*i) = pntList->at(i).x;
        (*q)(2*i+1) = pntList->at(i).y;

        (*qt)[2*i]   = obj.getU0();
        (*qt)[2*i+1] = obj.getV0();

        (*qtBackup)[2*i]   = obj.getU0();
        (*qtBackup)[2*i+1] = obj.getV0();
    }

    // 0 out the forces after intiialization
    f->setZero();

    if (this->updateMode == 1) {
        // If using semi-implicit method, setup the matrix.
        this->matrix = this->createStiffnessMatrix();

        // Set up the parameters
        this->params = setUpParams();

        pcgTol = new double[2*pntList->size()];
        pcgRHS = new double[2*pntList->size()];

        simutils::set_constant(2*pntList->size(), 0.0, pcgTol);
        simutils::set_constant(2*pntList->size(), 0.0, pcgRHS);

        matrix->set_toler(pcgTol);
    } else if (this->updateMode == 2) {
        // If using ADMM mode, set up the assembly parameters

        // Set the number of points in the assembly
        setNPnts(pntList->size());
        setD(2);

        // Build the diagonal mass matrix
        Eigen::VectorXd m(Eigen::VectorXd::Constant(2*pntList->size(), pntMass));
        buildMassMatrix(m);

        // Build the W matrix
        buildWMatrix(this->w);

        // Build the D matrix (just the identity)
        buildDMatrix();

        // Create the solver
        if (pool.dtFix > 0) {
            admmSolver = new ADMMPG(pool.dtFix, *this);
        } else {
            assert(false);
        }
    }

    // If this is a static object, set the vels to 0.
    if (objType == SolidObject::ObjectType::STATIC) {
        cout << "STATIC OBJECT BEING USED" << endl;
        qt->setZero();
        qtBackup->setZero();
    }

    // Data structures for collision handling
    nodeCols = new vector<set<massPoint2D*>>(pntList->size());
}

/**
 * Build the D matrix of size R^{d*nedges x d*nedges}
*/
void MassSpring2D::buildDMatrix() {
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
void MassSpring2D::buildMassMatrix(Eigen::VectorXd &m) {
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
 * Build the W matrix of size R^{d*nedges x d*nedges}
*/
void MassSpring2D::buildWMatrix(double w) {
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

ParamIter* MassSpring2D::setUpParams() {
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
MatrixIter* MassSpring2D::createStiffnessMatrix() {
    const int d = 2;

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
        // matrixBuilder = new MatrixStruc(nnz /*Number of unknowns*/,
        //                                 true /*Don't add non-zeros to diagonals so I don't get mindflooded*/);
        matrixBuilder = new MatrixStruc(2*pntList->size() /*Number of unknowns*/,
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
            matrixBuilder->set_entry(off+1, off);
            matrixBuilder->set_entry(off+1, off+1);

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
                matrixBuilder->set_entry(off+1, noff);
                matrixBuilder->set_entry(off+1, noff+1);
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
 * Point projection algorithm.
*/
double MassSpring2D::projectOntoEdge(int pntId1, int pntId2, double x[2], double &t) {
    // Extract point coordinates
    double x1[2] = {pntList->at(pntId1).x, pntList->at(pntId1).y};
    double x2[2] = {pntList->at(pntId2).x, pntList->at(pntId2).y};

    // Edge vector
    double u[2] = {x2[0] - x1[0], x2[1] - x1[1]};

    // x - x1
    double w[2] = {x[0] - x1[0], x[1] - x1[1]};

    // Compute the vector projection from first edge point to line
    double coeff = simutils::ddot2d(u, w)/simutils::ddot2d(u, u);
    double proj[2] = {
        coeff*u[0],
        coeff*u[1]
    };

    // Compute the distance vector
    // proj - p
    double distVec[2];
    distVec[0] = proj[0] - w[0];
    distVec[1] = proj[1] - w[1];

    // Compute the distance
    double d;
    if (simutils::sign(proj[0]) == simutils::sign(u[0])     &&
            simutils::sign(proj[1]) == simutils::sign(u[1])) {
        // Compute t
        t = simutils::eucNorm2D(proj) / simutils::eucNorm2D(u);
        if (t > 1.0) {
            t = 1.0;
            d = simutils::eucDiff2D(x, x2);
        } else {
            d = simutils::eucNorm2D(distVec);
        }
    } else {
        t = 0.0;
        d = simutils::eucDiff2D(x, x1);
    }

    return d;
}

/** 
 * Find the nearest point on the boundary of the MSS
*/
double MassSpring2D::closestBoundaryDist(double inPnt[2]) {
    double d, t;
    double dmin = INFINITY;

    int pntIds[2];

    // Iterate through all of the boundary edges of the MSS
    for (int i = 0; i < boundaryEdgeIdList->size(); i++) {
        pntIds[0] = edgeList->at(boundaryEdgeIdList->at(i)).pntIds[0];
        pntIds[1] = edgeList->at(boundaryEdgeIdList->at(i)).pntIds[1];

        d = projectOntoEdge(pntIds[0], pntIds[1], inPnt, t);

        // If this is the new minimum distance, update the output velocity
        if (d <= dmin) {
            dmin = d;
        }
    }

    return dmin;
}

/** 
 * Find the nearest point on the boundary of the MSS
*/
double MassSpring2D::closestBoundaryPnt(double inPnt[2], double out[2]) {
    double d, t;
    double dmin = INFINITY;

    int pntIds[2];

    // Iterate through all of the boundary edges of the MSS
    for (int i = 0; i < boundaryEdgeIdList->size(); i++) {
        pntIds[0] = edgeList->at(boundaryEdgeIdList->at(i)).pntIds[0];
        pntIds[1] = edgeList->at(boundaryEdgeIdList->at(i)).pntIds[1];

        d = projectOntoEdge(pntIds[0], pntIds[1], inPnt, t);

        // If this is the new minimum distance, update the output velocity
        if (d <= dmin) {
            dmin = d;
            double pnt1[2] = {pntList->at(pntIds[0]).x, pntList->at(pntIds[0]).y};
            double pnt2[2] = {pntList->at(pntIds[1]).x, pntList->at(pntIds[1]).y};


            out[0] = ((1-t)*pnt1[0] + t*pnt2[0]);
            out[1] = ((1-t)*pnt1[1] + t*pnt2[1]);
        }
    }

    return dmin;
}

/**
 * Find the closest point on the boundary of the MSS to the input point.
*/
void MassSpring2D::interpFaceVels(double inPnt[2], double out[2]) {
    double d;
    double dmin = INFINITY;
    double t;

    double pnt1_u, pnt1_v, pnt2_u, pnt2_v;

    int pntIds[2];

    // Iterate through all of the boundary edges of the MSS
    for (int i = 0; i < boundaryEdgeIdList->size(); i++) {
        pntIds[0] = edgeList->at(boundaryEdgeIdList->at(i)).pntIds[0];
        pntIds[1] = edgeList->at(boundaryEdgeIdList->at(i)).pntIds[1];

        // Coordinates of the end points
        pnt1_u = pntList->at(pntIds[0]).u;
        pnt1_v = pntList->at(pntIds[0]).v;

        pnt2_u = pntList->at(pntIds[1]).u;
        pnt2_v = pntList->at(pntIds[1]).v;

        // Get the distance, and t value (by reference)
        d = projectOntoEdge(pntIds[0], pntIds[1], inPnt, t);

        // If this is the new minimum distance, update the output velocity
        if (d <= dmin) {
            out[0] = (1-t)*pnt1_u + t*pnt2_u;
            out[1] = (1-t)*pnt1_v + t*pnt2_v;
            dmin = d;
        }
    }
}

/**
 * Return whether the line segment defined by two imput points intersects the MSS boundary
*/
bool MassSpring2D::intersectsBoundary(double pnt1[2], double pnt2[2]) {
    double edgePnt1[2];
    double edgePnt2[2];
    int ids[2];

    for (auto edgeId = boundaryEdgeIdList->begin(); edgeId != boundaryEdgeIdList->end(); ++edgeId) {
        ids[0] = edgeList->at(*edgeId).pntIds[0];
        ids[1] = edgeList->at(*edgeId).pntIds[1];

        edgePnt1[0] = pntList->at(ids[0]).x;
        edgePnt1[1] = pntList->at(ids[0]).y;

        edgePnt2[0] = pntList->at(ids[1]).x;
        edgePnt2[1] = pntList->at(ids[1]).y;

        if (simutils::line_intersect(pnt1, pnt2, edgePnt1, edgePnt2)) {
            return true;
        }
    }

    return false;
}

/**
 * Find the interface point by searching along the "grid normal" direction until an interface is 
 * found.
 * 
 * Returns whether a point was found, mutates nearestInterfacePont
*/
bool MassSpring2D::findNearestGridNormalInterface(Pool2D &pool, int pntId, int nearestInterfacePnt[2]) {
    // cout << "Computing the gradient" << endl;
    int nx = pool.getNx();
    int ny = pool.getNy();

    // Compute the location on the grid of the input point
    int i = simutils::findLimInfMeshPoint(pntList->at(pntId).x, pool.getXMesh(), nx+1);
    int j = simutils::findLimInfMeshPoint(pntList->at(pntId).y, pool.getYMesh(), ny+1);

    // Compute the normal vector
    double n[2];
    pool.levelSetGradient(i, j, n);
    simutils::normalize2D(n);

    int nb[2];

    // Use a trick to compute the block normal vector
    if (simutils::in_range(n[0], cos(5.0*M_PI/8.0), cos(3.0*M_PI/8.0))) {
        nb[0] = 0;
    } else if (n[0] >= 0.0) {
        nb[0] = 1;
    } else {
        nb[0] = -1;
    }

    if (simutils::in_range(n[1], sin(9.0*M_PI/8.0), cos(7.0*M_PI/8.0))) {
        nb[1] = 0;
    } else if (n[0] >= 0.0) {
        nb[1] = 1;
    } else {
        nb[1] = -1;
    }

    // Search along this direction until the interface point is found
    objects::FSIObject obj = pool.objAtIndex(i, j);
    bool intFound = false;

    int so = 0;
    int count = 0;
    int jo, io;
    while (!intFound && simutils::int_in_range(i+so*nb[0], 0, nx-1) && simutils::int_in_range(j+so*nb[1], 0, ny-1)) {
        io = i+so*nb[0];
        jo = j+so*nb[1];

        obj = pool.objAtIndex(io, jo);
        intFound = pool.isInterface(obj);

        /* Do 4-connected band-search (there will be some overlap but its ok) */

        // left
        if (!intFound && simutils::int_in_range(io-1, 0, nx-1) && simutils::int_in_range(jo, 0, ny-1)) {
            obj = pool.objAtIndex(io-1, jo);
            intFound = pool.isInterface(obj);
        }

        // right
        if (!intFound && simutils::int_in_range(io+1, 0, nx-1) && simutils::int_in_range(jo, 0, ny-1)) {
            obj = pool.objAtIndex(io+1, jo);
            intFound = pool.isInterface(obj);
        }

        // down
        if (!intFound && simutils::int_in_range(io, 0, nx-1) && simutils::int_in_range(jo-1, 0, ny-1)) {
            obj = pool.objAtIndex(io, jo-1);
            intFound = pool.isInterface(obj);
        }

        // up
        if (!intFound && simutils::int_in_range(io, 0, nx-1) && simutils::int_in_range(jo+1, 0, ny-1)) {
            obj = pool.objAtIndex(io, jo+1);
            intFound = pool.isInterface(obj);
        }

        if (!intFound) {
            count++;
            so += ( (count % 2) ? 1 : -1 )*count;
        }
    }

    nearestInterfacePnt[0] = 0;
    nearestInterfacePnt[1] = 0;

    if (intFound) {
        nearestInterfacePnt[0] = i+so*nb[0];
        nearestInterfacePnt[1] = j+so*nb[1];
    } else {
        nearestInterfacePnt[0] = 0;
        nearestInterfacePnt[1] = 0;

        cout << "Could not find stess to extrapolate on MSS" << endl;
        assert(false);
    }

    return intFound;
}

/**
 * Get the number of boundary neighbours of a boundary mass point
*/
int MassSpring2D::getNumNeighs(mass_spring::massPoint2D pnt) {
    int count = 0;

    for (auto edge = pnt.edgeIds.begin(); edge != pnt.edgeIds.end(); ++edge) {
        if (edgeList->at(*edge).boundaryEdge) {
            count++;
        }
    }

    return count;
}

/**
 * Calculate the stress due to boundary points not satisfying the interpolation condition
*/
void MassSpring2D::computeInterpolationStress(Pool2D &pool, double x[2], double stiffness, double out[2]) {

    // Compute the location of the nearest interface point
    double p[2];
    pool.getClosestInterfacePoint(x, p);

    // double hx = pool.getXMeshVal(1) - pool.getXMeshVal(0);
    // double hy = pool.getYMeshVal(1) - pool.getYMeshVal(0);

    // double h = sqrt(hx*hx + hy*hy);

    // stiffness = (abs(simutils::eucDiff2D(p, x)) > h) ? exp(simutils::eucDiff2D(p, x)) : 0.0;

    // Compute the spring stress
    out[0] = -stiffness*(p[0] - x[0]);
    out[1] = -stiffness*(p[1] - x[1]);

}

/**
 * Convert mass point to structure location
*/
void MassSpring2D::structToLoc(mass_spring::massPoint2D pnt, double loc[2]) {
    loc[0] = pnt.x;
    loc[1] = pnt.y;
}

/**
 * Computes the collision stress for a given boundary node, specified by its ID
*/
void MassSpring2D::computeCollisionStress(int nodeId, double colStress[2]) {
    // 0 out the stress vector
    colStress[0] = 0.0;
    colStress[1] = 0.0;

    // Compute the forces acting internally on this node by the other members of the MSS
    massPoint2D mPnt = pntList->at(nodeId);

    int neighId = -1;
    int *pntIds;
    double forces[4];
    double nodeForces[2] = {0.0, 0.0};
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
    }

    // For each of the collisions this node is involved in, calculate a resisting force
    massPoint2D colPnt;
    double collisionNormal[2] = {0.0, 0.0};
    double pntDiff[2];
    double pntDist;
    int numNear = 0;
    for (auto near = nodeCols->at(nodeId).begin(); near != nodeCols->at(nodeId).end(); ++near) {
        colPnt = **near;

        // Compute the normal vector for this collision
        collisionNormal[0] = mPnt.x - colPnt.x;
        collisionNormal[1] = mPnt.y - colPnt.y;
        simutils::normalize2D(collisionNormal);

        // Compute the collision stress to stop the momentum in this direction
        double fac = simutils::ddot2d(nodeForces, collisionNormal);
        colStress[0] += - fac*collisionNormal[0];
        colStress[1] += - fac*collisionNormal[1];

        // Compute the spring repulsion force for this node
        pntDiff[0] = mPnt.x - colPnt.x;
        pntDiff[1] = mPnt.y - colPnt.y;
        pntDist = simutils::eucNorm2D(pntDiff);

        if (collisionDist > pntDist) {
            calcElasticForce(this->E, collisionDist, mPnt, colPnt, forces);
        } else {
            forces[0] = 0.0;
            forces[1] = 0.0;
        }

        colStress[0] += forces[0];
        colStress[1] += forces[1];

        numNear++;
    }

    colStress[0] /= numNear;
    colStress[1] /= numNear;
    
    cout << "colStress = (" << colStress[0] << ", " << colStress[1] << ")" << endl;
}

/**
 * Set the tolerance for the ADMM algorithm.
*/
void MassSpring2D::setAdmmTol(double tol) {
    this->admmTol = tol;
}

/** 
 * Helper function to apply external force to a boundary node.
*/
void MassSpring2D::applyBoundaryForces(Pool2D &pool, double ***stress, int ng, double fNet[2]) {

    int ni, nj; // Locations of the nodes

    for (auto pnt = pntList->begin(); pnt != pntList->end(); ++pnt) {
        pnt->sigU = 0.0;
        pnt->sigV = 0.0;
    }

    int iPnt[2];

    // Compute the net force acting on the edges
    massPoint2D mPnt1;
    massPoint2D mPnt2;
    int id1 = 0;
    int id2 = 0;
    double s1[2];
    double s2[2];
    double diffNorm;
    for (auto edgeId = boundaryEdgeIdList->begin(); edgeId != boundaryEdgeIdList->end(); ++edgeId) {
        id1 = edgeList->at(*edgeId).pntIds[0];
        id2 = edgeList->at(*edgeId).pntIds[1];
        mPnt1 = pntList->at(id1);
        mPnt2 = pntList->at(id2);

        // Distance between the two points.
        diffNorm = sqrt(simutils::square(mPnt2.x - mPnt1.x) + simutils::square(mPnt2.y - mPnt1.y));

        // Find stress acting on first node
        bool found1 = findNearestGridNormalInterface(pool, id1, iPnt);
        ni = iPnt[0];
        nj = iPnt[1];

        if (found1) {
            // Apply the stress to this point and its connected neighbours

            if (nodeCols->at(id1).size() > 0) {
                // There is a collision on this node, compute the collision stress
                computeCollisionStress(id1, s1);
            } else {
                // Apply the hydrodynamic stress if there is no collision
                s1[0] = stress[0][ng+nj][ng+ni];
                s1[1] = stress[1][ng+nj][ng+ni];
            }
        } else {
            assert(false);
            s1[0] = 0.0;
            s1[1] = 0.0;
        }

        bool found2 = findNearestGridNormalInterface(pool, id2, iPnt);
        ni = iPnt[0];
        nj = iPnt[1];

        if (found2) {
            // // Apply the stress to this point and its connected neighbours
            // s2[0] = stress[0][ng+nj][ng+ni];
            // s2[1] = stress[1][ng+nj][ng+ni];
            if (nodeCols->at(id2).size() > 0) {
                // There is a collision on this node, compute the collision stress
                computeCollisionStress(id2, s2);
            } else {
                // Apply the hydrodynamic stress
                // Hydrodynamic stress
                s2[0] = stress[0][ng+nj][ng+ni];
                s2[1] = stress[1][ng+nj][ng+ni];
            }
        } else {
            assert(false);
            s2[0] = 0.0;
            s2[1] = 0.0;
        }

        // Update the stresses in the MSS
        (*f)[2*id1] += 0.5*(s1[0] + s2[0])*diffNorm;
        (*f)[2*id1+1] += 0.5*(s1[1] + s2[1])*diffNorm;

        (*f)[2*id2] += 0.5*(s1[0] + s2[0])*diffNorm;
        (*f)[2*id2+1] += 0.5*(s1[1] + s2[1])*diffNorm;
    }

    for (int id = 0; id < pntList->size(); id++) {
        pntList->at(id).sigU = (*f)[2*id];
        pntList->at(id).sigV = (*f)[2*id+1];
    }

    // Compute the net force acting on the edges
    fNet[0] = 0.0;
    fNet[1] = 0.0;
    // double diffNorm;
    for (auto edgeId = boundaryEdgeIdList->begin(); edgeId != boundaryEdgeIdList->end(); ++edgeId) {
        mPnt1 = pntList->at(edgeList->at(*edgeId).pntIds[0]);
        mPnt2 = pntList->at(edgeList->at(*edgeId).pntIds[1]);

        // Distance between the two points.
        diffNorm = sqrt(simutils::square(mPnt2.x - mPnt1.x) + simutils::square(mPnt2.y - mPnt1.y));

        // Compute contribution to the net force acting on this object
        fNet[0] += 0.5*(mPnt1.sigU + mPnt2.sigU)*diffNorm;
        fNet[1] += 0.5*(mPnt1.sigV + mPnt2.sigV)*diffNorm;
    }
}

/**
 * Add to the Hessian for simple spring system
 * 
 * Notes are on page 24 of fund var mesh adapt
 * 
 * // NOTE: assumes that the mass per node is constant (good assumption)
*/
void MassSpring2D::calcLocalElasticHessian(double dt, edge2D edge, int pntId1,
                                           massPoint2D pnt1, int pntId2,
                                           massPoint2D pnt2) {
    assert(pntId1 < pntId2);

    double diff[2] = {pnt2.x - pnt1.x, pnt2.y - pnt1.y};
    double diffNorm = simutils::eucNorm2D(diff);

    int rOff1 = 2*pntId1;
    int rOff2 = 2*pntId2;

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
        } else if (colIndex == rOff2) {
            matrix->aValue(i) += coeff1*-diff[0]*diff[0] - coeff2;
        } else if (colIndex == rOff2+1) {
            matrix->aValue(i) += coeff1*-diff[0]*diff[1];
        }
    }

    for (int i = matrix->rowBegin(rOff1+1); i < matrix->rowEndPlusOne(rOff1+1); i++) {
        colIndex = matrix->getColIndex(i);

        if (colIndex == rOff1) {
            matrix->aValue(i) += coeff1*diff[1]*diff[0];
        } else if (colIndex == rOff1+1) {
            matrix->aValue(i) += coeff1*diff[1]*diff[1] + coeff2;
        } else if (colIndex == rOff2) {
            matrix->aValue(i) += coeff1*-diff[1]*diff[0];
        } else if (colIndex == rOff2+1) {
            matrix->aValue(i) += coeff1*-diff[1]*diff[1] - coeff2;
        }
    }

    // r = i(1), c = i(0) & r = i(1), c = i(1)
    for (int i = matrix->rowBegin(rOff2); i < matrix->rowEndPlusOne(rOff2); i++) {
        colIndex = matrix->getColIndex(i);

        if (colIndex == rOff1) {
            matrix->aValue(i) += coeff1*diff[0]*-diff[0] - coeff2;
        } else if (colIndex == rOff1+1) {
            matrix->aValue(i) += coeff1*diff[0]*-diff[1];
        } else if (colIndex == rOff2) {
            matrix->aValue(i) += coeff1*diff[0]*diff[0] + coeff2;
        } else if (colIndex == rOff2+1) {
            matrix->aValue(i) += coeff1*diff[0]*diff[1];
        }
    }

    for (int i = matrix->rowBegin(rOff2+1); i < matrix->rowEndPlusOne(rOff2+1); i++) {
        colIndex = matrix->getColIndex(i);

        if (colIndex == rOff1) {
            matrix->aValue(i) += coeff1*diff[1]*-diff[0];
        } else if (colIndex == rOff1+1) {
            matrix->aValue(i) += coeff1*diff[1]*-diff[1] - coeff2;
        } else if (colIndex == rOff2) {
            matrix->aValue(i) += coeff1*diff[1]*diff[0];
        } else if (colIndex == rOff2+1) {
            matrix->aValue(i) += coeff1*diff[1]*diff[1] + coeff2;
        }
    }
}

/**
 * General elastic force calculation, as we use this in several places in the code
*/
void MassSpring2D::calcElasticForce(double E, double l0, massPoint2D pnt1,
                                    massPoint2D pnt2, double force[4]) {
                                    
    double diff[2] = {pnt2.x - pnt1.x, pnt2.y - pnt1.y};
    double diffNorm = simutils::eucNorm2D(diff);

    force[0] = -E*((diffNorm - l0)/diffNorm)*(-diff[0]);
    force[1] = -E*((diffNorm - l0)/diffNorm)*(-diff[1]);

    force[2] = -E*((diffNorm - l0)/diffNorm)*(diff[0]);
    force[3] = -E*((diffNorm - l0)/diffNorm)*(diff[1]);
}

/**
 * Add to the local potential energy for simple spring system. (note, when summing everything
 * together is when we apply the (-) sign)!
 * 
 * Note: must have pntId1 < pntId2
*/
void MassSpring2D::calcLocalElasticForce(edge2D edge, int pntId1,
        massPoint2D pnt1, int pntId2, massPoint2D pnt2) {
    assert(pntId1 < pntId2);

    double force[4];

    calcElasticForce(E, edge.l0, pnt1, pnt2, force);

    (*f)(2*pntId1)   += force[0];
    (*f)(2*pntId1+1) += force[1];

    (*f)(2*pntId2)   += force[2];
    (*f)(2*pntId2+1) += force[3];

}

/**
 * Add to the local potential energy for spring-dashpot system. (note, when summing everything
 * together is when we apply the (-) sign)!
 * 
 * Note: must have pntId1 < pntId2
*/
void MassSpring2D::calcLocalKelvinForce(edge2D edge, int pntId1, massPoint2D pnt1, int pntId2, massPoint2D pnt2) {
    assert(pntId1 < pntId2);

    double diff[2] = {pnt2.x - pnt1.x, pnt2.y - pnt1.y};
    double diffU[2] = {pnt2.u - pnt1.u, pnt2.v - pnt1.v};

    double diffNorm = simutils::eucNorm2D(diff);
    double diffNormCube = diffNorm*diffNorm*diffNorm;
    double innerProd = -pnt1.x*(pnt2.u - pnt1.u) - pnt1.y*(pnt2.v - pnt1.v)
                    + pnt2.x*(pnt2.u - pnt1.u) + pnt2.y*(pnt2.v - pnt1.v);
    
    if (simutils::eps_equal(diffNorm, 0.0, 1e-16)) {
        cout << "ERROR: calculating force for MSS" << endl;
        cout << "Points with ids " << pntId1 << " and " << pntId2 << " are too close!" << endl;
        cout << "diffNorm = " << diffNorm << endl;
        assert(false);

    }

    // Compute the force for the first input point (note we must take the sum)
    (*f)[2*pntId1] += -(E*((diffNorm - edge.l0)/diffNorm)*(-diff[0]) 
        + eta*( (innerProd/simutils::square(diffNorm))*(-diff[0])
                    + (diffNorm - edge.l0)/diffNorm*(-diffU[0])
                    - (diffNorm - edge.l0)*innerProd/diffNormCube*(-diff[0]) ));
    (*f)[2*pntId1+1] += -(E*((diffNorm - edge.l0)/diffNorm)*(-diff[1]) 
        + eta*( (innerProd/simutils::square(diffNorm))*(-diff[1])
                    + (diffNorm - edge.l0)/diffNorm*(-diffU[1])
                    - (diffNorm - edge.l0)*innerProd/diffNormCube*(-diff[1]) ));

    (*f)[2*pntId2] += -(E*((diffNorm - edge.l0)/diffNorm)*(diff[0]) 
        + eta*( (innerProd/simutils::square(diffNorm))*(diff[0])
                    + (diffNorm - edge.l0)/diffNorm*(diffU[0])
                    - (diffNorm - edge.l0)*innerProd/diffNormCube*(diff[0]) ));
    (*f)[2*pntId2+1] += -(E*((diffNorm - edge.l0)/diffNorm)*(diff[1]) 
        + eta*( (innerProd/simutils::square(diffNorm))*(diff[1])
                    + (diffNorm - edge.l0)/diffNorm*(diffU[1])
                    - (diffNorm - edge.l0)*innerProd/diffNormCube*(diff[1]) ));
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
void MassSpring2D::interpolateBoundaryLocation(Pool2D &pool, int i, int j, double X[2]) {
    // Get the outward normal direction for this point
    int nDir[2];
    pool.getNormalDir(pool.objAtIndex(i, j), nDir);

    // Negate to find the inward normal
    nDir[0] *= -1;
    nDir[1] *= -1;

    // Coordinate of the current boundary point
    double bndPnt[2] = {simutils::midpoint(pool.getXMeshVal(i), pool.getXMeshVal(i+1)),
                        simutils::midpoint(pool.getYMeshVal(j), pool.getYMeshVal(j+1))};
    
    // Coordinates of the structure point in the neighbouring direction
    double structPnt[2] = {simutils::midpoint(pool.getXMeshVal(i+nDir[0]), pool.getXMeshVal(i+nDir[0]+1)),
                        simutils::midpoint(pool.getYMeshVal(j+nDir[1]), pool.getYMeshVal(j+nDir[1]+1))};
    
    // Compute the unit inward normal
    double inN[2] = {(double)nDir[0], (double)nDir[1]};
    simutils::normalize2D(inN);

    // Take the phi values at these points
    double phiOut = pool.getPhiVal(i, j);
    double phiIn = pool.getPhiVal(i+nDir[0], j+nDir[1]);

    // Compute the distance to the interface
    double distVec[2] = {bndPnt[0] - structPnt[0], bndPnt[1] - structPnt[1]};
    double d = abs(phiOut/(phiIn - phiOut))*simutils::eucNorm2D(distVec);

    // Finally, compute the interface point location
    X[0] = bndPnt[0] + d*inN[0];
    X[1] = bndPnt[1] + d*inN[1];
}

/**
 * Update the locations and velocities of the spring nodes using Euler's method
*/
void MassSpring2D::eulerSolve(double dt, int elementMode, bool initMode) {
    int pntId1;
    int pntId2;
    massPoint2D pnt1;
    massPoint2D pnt2;

    double ETemp = E;
    double etaTemp = eta;
    
    if (initMode) {
        E = 1.0;
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
            (*qt)[2*i] += (dt/pntList->at(i).mass)*(*f)[2*i];
            (*q)[2*i] += dt*(*qt)[2*i];

            (*qt)[2*i+1] += (dt/pntList->at(i).mass)*(*f)[2*i+1];
            (*q)[2*i+1] += dt*(*qt)[2*i+1];
        } else {
            (*qt)[2*i] = 0.0;
            (*qt)[2*i+1] = 0.0;
        }
    }

    E = ETemp;
    eta = etaTemp;
}

/**
 * Update according to the linearly implicit time stepping method using BICGSTAB
 * algorithm to solve the symmetric system of equations
*/
void MassSpring2D::linearImplicitSolve(double dt, int elementMode, bool initMode) {
    int pntId1;
    int pntId2;
    massPoint2D pnt1;
    massPoint2D pnt2;

    double ETemp = E;
    double etaTemp = eta;
    
    if (initMode) {
        E = 1.0;
        eta = 0.0;
    }

    // simutils::copyVals(2*pntList->size(), qt, qprev);
    *qprev = *qt;

    // 0 out the Hessian
    for (int r = 0; r < 2*pntList->size(); r++) {
        for (int i = matrix->rowBegin(r); i < matrix->rowEndPlusOne(r); i++) {
            matrix->aValue(i) = 0.0;
        }
    }

    // Initialize the linear system RHS
    simutils::set_constant(2*pntList->size(), 0.0, pcgRHS);

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
    double accum = 0;
    double count = 0;

    for (int r = 0; r < 2*pntList->size(); r++) {
        for (int i = matrix->rowBegin(r); i < matrix->rowEndPlusOne(r); i++) {
            colIndex = matrix->getColIndex(i);

            if (colIndex == r) {
                // matrix->aValue(i) = 1.0;
                matrix->aValue(i) = 1 + coeff*matrix->aValue(i);
            } else {
                // matrix->aValue(i) = -coeff*matrix->aValue(i);
                matrix->aValue(i) = coeff*matrix->aValue(i);
                // matrix->aValue(i) = 0.0;
            }

            accum += matrix->aValue(i);
            count ++;
        }
    }

    /* Setup the RHS vector (use Euler step as initial guess) */
    double mass = pntList->at(0).mass;
    for (int i = 0; i < 2*pntList->size(); i++) {
        matrix->bValue(i) = (*qt)[i] + (dt/mass)*(*f)[i] + 1e-16;
        pcgRHS[i] = matrix->bValue(i);
    }

    /* Solve the linear system */
    try {
        int niter=0; // The number of iterations done in the solve

        // Perform the ILU factorization
        if (params->drop_ilu == 0) {
            matrix->sfac(*params);
        }

        simutils::set_constant(2*pntList->size(), 0.0, pcgTol);

        matrix->set_toler(this->pcgTol); // Set the update tolerance

        // Solve the linear system. Solution is stored in provided vector
        matrix->solve(*this->params, pcgRHS, niter);

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
    simutils::copyVals(2*pntList->size(), pcgRHS, qt->data());

    // Do the Euler step to get the new node positions
    *q += dt*(*qt);

    E = ETemp;
    eta = etaTemp;
}

/**
 * Update according to the verlet method
*/
void MassSpring2D::verletSolve(double dt, int elementMode, bool initMode) {
    int pntId1;
    int pntId2;
    massPoint2D pnt1;
    massPoint2D pnt2;

    double ETemp = E;
    double etaTemp = eta;

    if (initMode) {
        E = 5.0;
        eta = 0.1;
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
    double temp1, temp2;
    for (int i = 0; i < pntList->size(); i++) {
        if (iterCount == 0) {
            // Backup current position
            temp1 = (*q)[2*i];
            temp2 = (*q)[2*i+1];

            // Update the current position
            if (!(initMode && pntList->at(i).boundaryPnt))  {
                (*q)[2*i] += dt*((*qt)[2*i]) + 0.5*simutils::square(dt)*((*f)[2*i]/pntList->at(i).mass);
                (*q)[2*i+1] += dt*((*qt)[2*i+1]) + 0.5*simutils::square(dt)*((*f)[2*i+1]/pntList->at(i).mass);
            }

            // Update the solution history.
            (*qprev)[2*i] = temp1;
            (*qprev)[2*i+1] = temp2;

            // Update the velocities
            (*qt)[2*i] = ((*q)[2*i] - (*qprev)[2*i])/dt;
            (*qt)[2*i+1] = ((*q)[2*i+1] - (*qprev)[2*i+1])/dt;

        } else {
            temp1 = (*q)[2*i];
            temp2 = (*q)[2*i+1];

            // Update the current positions
            if (!(initMode && pntList->at(i).boundaryPnt)) {
                (*q)[2*i] = 2.0*(*q)[2*i] - (*qprev)[2*i] + (simutils::square(dt)/pntList->at(i).mass)*(*f)[2*i];
                (*q)[2*i+1] = 2.0*(*q)[2*i+1] - (*qprev)[2*i+1] + (simutils::square(dt)/pntList->at(i).mass)*(*f)[2*i+1];
                // (*q)[2*i] += dt*((*qt)[2*i]) + 0.5*simutils::square(dt)*((*f)[2*i]/pntList->at(i).mass);
                // (*q)[2*i+1] += dt*((*qt)[2*i+1]) + 0.5*simutils::square(dt)*((*f)[2*i+1]/pntList->at(i).mass);
                // (*q)[2*i] = (*q)[2*i] + (*qt)[2*i] + (simutils::square(dt)/pntList->at(i).mass)*(*f)[2*i];
                // (*q)[2*i+1] = (*q)[2*i+1] + (*qt)[2*i+1] + (simutils::square(dt)/pntList->at(i).mass)*(*f)[2*i+1];
            }

            // Update the velocities
            (*qt)[2*i] = ((*q)[2*i] - (*qprev)[2*i])/(2.0*dt);
            (*qt)[2*i+1] = ((*q)[2*i+1] - (*qprev)[2*i+1])/(2.0*dt); 

            // Update the solution history
            (*qprev)[2*i] = temp1;
            (*qprev)[2*i+1] = temp2;
        }
    }

    E = ETemp;
    eta = etaTemp;
}


/**
 * Update the state of the mass-spring system using the information provided by the fluid
 * 
 * NOTE: stress does not have any ghost cells in the first axis, stress = [PoolU, PoolV]
 *       inherited from the pool class.
*/
void MassSpring2D::updateSolidVels(double dt, Pool2D &pool, double ***stress, double fNet[2],
                                int ng, bool initMode) {
    
    // In the case of a static object, do nothing.
    if (!initMode && objType == SolidObject::ObjectType::STATIC) {
        return;
    }

    this->initMode = initMode;
    int temp;

    if (this->initMode) {
        temp = this->updateMode;
        this->updateMode = 0;
    }

    // this->dt = dt;
    if (updateMode == 2) {
        this->dt = this->dtFix;
    } else {
        this->dt = dt;
    }

    // Initialize the force vector to 0.
    f->setZero();

    // Set previous solution to current
    *qprev = *q;

    // Loop through all of the points. For each boundary point, add to the force vector.
    if (!initMode) {
        applyBoundaryForces(pool, stress, ng, fNet);
    }

    // Loop through all of the edges, using the potential energy to compute the displacement of the
    // nodes.
    if (updateMode == 0) {
        verletSolve(dt, elementMode, initMode);
    } else if (updateMode == 1) {
        linearImplicitSolve(dt, elementMode, initMode);
    } else if (updateMode == 2) {
        admmSolver->step(500, 1e-7);
    }

    // Update the x, y vels of each point in the point list using the
    // updated q. (position should not be updated for each point)
    for (int i = 0; i < pntList->size(); i++) {
        pntList->at(i).u = (*qt)[2*i];
        pntList->at(i).v = (*qt)[2*i+1];
    }

    // Update the backups
    *qBackup     = *q;
    *qtBackup    = *qt;
    *qprevBackup = *qprev;

    iterCount++;

    if (this->initMode) {
        this->updateMode = temp;
    }
}

/**
 * Reset the function back to the previous time step
*/
void MassSpring2D::reset() {
    *q = *qBackup;
    *qt = *qtBackup;
    *qprev = *qprevBackup;

    for (int i = 0; i < pntList->size(); i++) {
        pntList->at(i).u = (*qt)[2*i];
        pntList->at(i).v = (*qt)[2*i+1];

        pntList->at(i).x = (*q)[2*i];
        pntList->at(i).y = (*q)[2*i+1];
    }
}

void MassSpring2D::interpBoundary(Pool2D &pool, bool resetRestingLengths) {

    double x, y, phi;
    double phiGrad[2];
    for (auto pnt = boundaryNodeIdList->begin(); pnt != boundaryNodeIdList->end(); ++pnt) {
        x = pntList->at(*pnt).x;
        y = pntList->at(*pnt).y;

        pool.interpolatePhiGrad(x, y, phiGrad);
        simutils::normalize2D(phiGrad);

        phi = pool.interpolatePhi(x, y);

        pntList->at(*pnt).x -= phi*phiGrad[0];
        pntList->at(*pnt).y -= phi*phiGrad[1];

        // Update the node locations.
        (*q)(2*(*pnt))       = pntList->at(*pnt).x;
        (*q)(2*(*pnt)+1)     = pntList->at(*pnt).y;

        (*qprev)(2*(*pnt))   = pntList->at(*pnt).x;
        (*qprev)(2*(*pnt)+1) = pntList->at(*pnt).y;
    }

    // Set the resting lengths on the edges to what they are with the updated boundary positions
    if (resetRestingLengths) {
        double pnt1[2];
        double pnt2[2];
        for (auto edge = edgeList->begin(); edge != edgeList->end(); ++edge) {
            pnt1[0] = pntList->at(edge->pntIds[0]).x;
            pnt1[1] = pntList->at(edge->pntIds[0]).y;

            pnt2[0] = pntList->at(edge->pntIds[1]).x;
            pnt2[1] = pntList->at(edge->pntIds[1]).y;

            edge->l0 = simutils::eucDiff2D(pnt1, pnt2);
        }
    }
}

/**
 * Assuming that updateSolidVels has already been called, update the solid locations based on these
 * velocities
*/
void MassSpring2D::updateSolidLocs(Pool2D &pool, bool interp) {
    // cout << "updating solids locs" << endl;
    // for (int i = 0; i < pntList->size(); i++) {
    //     pntList->at(i).x = (*q)[2*i];
    //     pntList->at(i).y = (*q)[2*i+1];
    // }
    // return;
    // if (!initMode && objType == SolidObject::ObjectType::STATIC) {
    //     return;
    // }


    // Let the boundary points interpolate the pool.
    if (interp) {

        // 0 out the boundary forces
        // simutils::set_constant(2*pntList->size(), 0.0, f);

        // simutils::copyVals(2*pntList->size(), qprev, qt);
        *qt = *qprev;

        // Add the interpolation forces for each of the boundary points
        // double pnt[2];
        // double interpStress[2];
        // for (auto bPnt = boundaryNodeIdList->begin(); bPnt != boundaryNodeIdList->end(); ++bPnt) {
        //     structToLoc(pntList->at(*bPnt), pnt);

        //     // Compute the stress due to the interpolation constraint on the MSS
        //     computeInterpolationStress(pool, pnt, 1, interpStress);

        //     // Add the interpolation stress to the force vector
        //     (*f)[2*(*bPnt)]   -= interpStress[0];
        //     (*f)[2*(*bPnt)+1] -= interpStress[1];
        // }

        // Perform solve
        linearImplicitSolve(dt, elementMode, initMode);

        for (int i = 0; i < 2*pntList->size(); i++) {
            (*q)[i] += dt*(*qt)[i];
        }

        for (int i = 0; i < pntList->size(); i++) {
            pntList->at(i).u = (*qt)[2*i];
            pntList->at(i).v = (*qt)[2*i+1];

            pntList->at(i).x = (*q)[2*i];
            pntList->at(i).y = (*q)[2*i+1];
        }
    } else {


        
        // for (int i = 0; i < 2*pntList->size(); i++) {
        //     (*q)[i] += dt*(*qt)[i];
        // }

        for (int i = 0; i < pntList->size(); i++) {
            pntList->at(i).x = (*q)[2*i];
            pntList->at(i).y = (*q)[2*i+1];
        }
        // assert(false);
    }
}

/**
 * Output the nodes of the MSS with no explicit knowledge of the connectivity between them.
*/
void MassSpring2D::outputNodes(const char* fname) {

    ofstream outFile;
    outFile.open(fname);

    for (auto pnt = pntList->begin(); pnt != pntList->end(); ++pnt) {
        // if (!pnt->boundaryPnt) {
        //     continue;
        // }
        outFile << pnt->x << ", " << pnt->y << ", " << pnt->sigU << ", " << pnt->sigV << endl;
    }

    outFile.close();
}

/**
 * Output the node locations and velocities.
*/
void MassSpring2D::outputNodeVels(const char* fname) {

    ofstream outFile;
    outFile.open(fname);

    for (auto pnt = pntList->begin(); pnt != pntList->end(); ++pnt) {
        outFile << pnt->x << ", " << pnt->y << ", " << pnt->u << ", " << pnt->v << endl;
    }

    outFile.close();
}

/**
 * Output the pairs of nodes.
*/
void MassSpring2D::outputEdges(const char* fname) {
    ofstream outFile;
    outFile.open(fname);

    massPoint2D pnt1;
    massPoint2D pnt2;

    for (auto edge = edgeList->begin(); edge != edgeList->end(); ++edge) {
        // if (!edge->boundaryEdge) {
        //     continue;
        // }
        pnt1 = pntList->at(edge->pntIds[0]);
        pnt2 = pntList->at(edge->pntIds[1]);

        outFile << pnt1.x << ", " << pnt1.y << endl;
        outFile << pnt2.x << ", " << pnt2.y << endl;
    }

    outFile.close();

}

void MassSpring2D::prox(double dt, Eigen::VectorXd &x, Eigen::VectorXd &DXpU, Eigen::VectorXd &z) {
    // First, make z the projection onto the constraints
    z = DXpU;

    // Note: z_i is interpreted as a spring length D_i x_i in the typical MSS way

    // For each edge, project the length constraint : TODO: make more general constraint projection... may want to bake this into the objects,
    // specifiying constraint objects which act on the input
    for (int i = 0; i < edgeList->size(); i++) {
        // Normalize the sub-vector to the resting length of the spring
        z.segment(i*d, d) *= edgeList->at(i).l0/(z.segment(i*d, d).norm());
    }

    z = (this->E*z + simutils::square(this->w)*DXpU)/(this->E + simutils::square(this->w));
}

void MassSpring2D::updateAfterStep(double dt, Eigen::VectorXd &xPrev, Eigen::VectorXd &x) {
    // Update the node positions
    // cout << "Diff after step " << (x - *q).norm() << endl;
    // assert(false);
    *q = x;
    
    // Update the node velocities
    *qt = (x - xPrev)/dt;
}

void MassSpring2D::copyX(Eigen::VectorXd &tar) {
    tar = *q;
}

void MassSpring2D::predictX(double dt, Eigen::VectorXd &xPrev, Eigen::VectorXd &x, Eigen::VectorXd &xBar) {
    // xBar = //2*x - xPrev + (simutils::square(dt)/pntMass)*(*f);
    xBar = x + dt*(*qt) + (simutils::square(dt)/pntMass)*(*f);
}


MassSpring2D::~MassSpring2D() {
    // cout << "in MSS destructor" << endl;
    delete edgeList;
    delete pntList;
    delete boundaryEdgeIdList;
    delete boundaryNodeIdList;

    // delete[] f;
    // delete[] q;
    // delete[] qt;
    // delete[] qprev;
    delete f;
    delete q;
    delete qt;
    delete qprev;
    delete qBackup;
    delete qtBackup;
    delete qprevBackup;

    if (updateMode == 1) {
        delete[] pcgRHS;
        delete[] pcgTol;
        delete matrix;
        delete params;
    } else if (updateMode == 2) {
        delete admmSolver;
    }

    delete nodeCols;
    // cout << "FINISHED in MSS destructor" << endl;
}