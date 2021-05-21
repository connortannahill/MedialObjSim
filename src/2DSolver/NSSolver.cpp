#include "MomentumSolver2D.h"
#include "Boundary.h"
#include "SolidObject.h"
#include "NSSolver.h"
#include "../Utils/SimUtilities.h"
#include "../Utils/Discretizations.h"
#include <cstdlib>
#include <cassert>
#include <vector>
#include "SimParams.h"

using namespace std;

/**
 * Class constructor. In this case we can simply call the parent class.
*/
NSSolver::NSSolver(Boundary &boundary,
                    vector<SolidObject> &solidObjects,
                    SimParams &params,
                    void (*initialConditions)(int,int,int,double*,double*,double**,double**),
                    void (*boundaryConditions)(int,int,double**,double**))
                    : MomentumSolver2D(boundary, solidObjects, params, initialConditions, boundaryConditions)
{
}

/**
 * Function to take a time step. Simply calls the step method in parent class
 * MomentumSolver2D
*/
double NSSolver::step(double tEnd, double safetyFactor) {
    return MomentumSolver2D::step(tEnd, safetyFactor);
}

/**
 * Function which takes explicit time step for the velocity terms in the
 * momentum equations.
 * Note: assuming no body forces
*/
void NSSolver::updateF(Pool2D *pool) {
    int i, j;

    // Extract the Reynolds number
    double Re = params->Re;

    int mo = methodOrd;

    // Compute Fu at the internal points.
    for (j = 0; j < this->ny; j++) {
        for (i = 0; i < this->nx-1; i++) {
            int yi = j + mo;
            int xi = i + mo;
            if (pool->objAtIndex(i, j) == objects::FLUID_C && pool->isUpdateableU(i, j)) {
                this->FU[yi][xi] = this->u[yi][xi]
                    + this->dt*( (1.0/Re)*(discs::firstOrder_lap_uxx(xi, yi, this->dx, this->u) 
                            + discs::firstOrder_lap_uyy(xi, yi, this->dy, this->u))
                            - discs::firstOrder_conv_usqx(xi, yi, this->dx, this->u)
                            - discs::firstOrder_conv_uvy(xi, yi, this->dy, this->u, this->v));
            }
        }
    }

    for (j = 0; j < this->ny-1; j++) {
        for (i = 0; i < this->nx; i++) {
            int yi = j + mo;
            int xi = i + mo;
            if (pool->objAtIndex(i, j) == objects::FLUID_C && pool->isUpdateableV(i, j)) {
                this->FV[yi][xi] = this->v[yi][xi]
                    + this->dt*( (1.0/Re)*(discs::firstOrder_lap_vxx(xi, yi, this->dx, this->v)
                            + discs::firstOrder_lap_vyy(xi, yi, this->dy, this->v))
                            - discs::firstOrder_conv_uvx(xi, yi, this->dx, this->u, this->v)
                            - discs::firstOrder_conv_vsqy(xi, yi, this->dy, this->v) );
            }
        }
    }

    // Apply the boundary conditions (domain boundary as well as the interfaces) using the pool
    // enumeration
    objects::FSIObject obj;
    int xi, yi;

    for (j = 0; j < ny; j++) {
        for (i = 0; i < nx; i++) {
            xi = methodOrd+i;
            yi = methodOrd+j;
            obj = pool->objAtIndex(i, j);
            if (pool->isInterface(obj)) {
                if (pool->hasStructInDir(obj, objects::NORTH)) {
                    FV[yi-1][xi] = v[yi-1][xi];
                }

                if (pool->hasStructInDir(obj, objects::SOUTH)) {
                    FV[yi][xi] = v[yi][xi];
                }

                if (pool->hasStructInDir(obj, objects::EAST)) {
                    FU[yi][xi-1] = u[yi][xi-1];
                }

                if (pool->hasStructInDir(obj, objects::WEST)) {
                    FU[yi][xi] = u[yi][xi];
                }
            }
        }
    }

    // Set the boundary values for F required in the pressure solver.
    // TODO: this should probably be moved into the pressure solver
    // so that it can be most consistent. TODO: find a way more efficient way of
    // implementing these checks (ofc after we have a prototype)
    for (j = 1; j <= this->ny; j++) {
        this->FU[j][0] = this->u[j][0];
        this->FU[j][this->nx] = this->u[j][this->nx];
    }

    for (i = 1; i <= this->nx; i++) {
        this->FV[0][i] = this->v[0][i];
        this->FV[this->ny][i] = this->v[this->ny][i];
    }
}

double NSSolver::getKinematicBoundaryLC(int i, int j, double velObj, double velNeigh,
                                        double cPnt[2], double stagPnt[2], int nDir[2]) {
    double hx = this->x[1] - this->x[0];
    double hy = this->y[1] - this->y[0];
    double nPnt[2];
    double inPnt[2];
    double distVec[2];
    double gamPnt[2];
    double phi_out, phi_in, dist;
    double t;

    // Take the unit outward normal
    double outN[2] = {(double) nDir[0], (double) nDir[1]};
    simutils::normalize2D(outN);

    nPnt[0] = stagPnt[0] + hx*(double)nDir[0];
    nPnt[1] = stagPnt[1] + hy*(double)nDir[1];

    inPnt[0] = stagPnt[0] - hx*(double)nDir[0];
    inPnt[1] = stagPnt[1] - hy*(double)nDir[1];

    // Evaluate the interpolant at the staggered point and the internal point
    phi_out = pool->interpolatePhi(stagPnt[0], stagPnt[1]);
    phi_in = pool->interpolatePhi(inPnt[0], inPnt[1]);

    if (phi_out > 0 && phi_in > 0) {
        // return fluid velocity
        return velNeigh;
    } else if (phi_out < 0 && phi_in < 0){
        return velObj;
    } else {
        // Compute distance to interface (from internal)
        distVec[0] = inPnt[0] - stagPnt[0];
        distVec[1] = inPnt[1] - stagPnt[1];
        dist = abs((phi_in/(phi_in - phi_out))*simutils::eucNorm2D(distVec));

        // Compute location of interface
        gamPnt[0] = inPnt[0] + dist*outN[0];
        gamPnt[1] = inPnt[1] + dist*outN[1];

        if (nDir[0] != 0) {
            t = (gamPnt[0] - inPnt[0])/(nPnt[0] - inPnt[0]);
        } else if (nDir[1] != 0) {
            t = (gamPnt[1] - inPnt[1])/(nPnt[1] - inPnt[1]);
        } else {
            t = 0;
        }

        if (!(t >= 0 && t <= 1)) {
            cout << endl;
            cout << "ERROR IN T" << endl;
            cout << "Phi_in " << phi_in << " Phi_out " << phi_out << endl;
            cout << "dist = " << dist << endl;
            cout << "t = " << t << endl;
            cout << "Normal: nx=" << nDir[0] << ", ny= " << nDir[1] << endl;
            cout << "cPnt: cx=" << cPnt[0] << ", cy= " << cPnt[1] << endl;
            cout << "sPnt: sx=" << stagPnt[0] << ", sy= " << stagPnt[1] << endl;
            cout << "nPnt: nx=" << nPnt[0] << ", ny= " << nPnt[1] << endl;
            cout << "nPnt: inx=" << inPnt[0] << ", iny= " << inPnt[1] << endl;
            cout << "gPnt: sx=" << gamPnt[0] << ", sy= " << gamPnt[1] << endl;
            cout << endl;
            cout << endl;
        }

        assert(t>=0 && t <= 1);

        return t*velObj + (1-t)*velNeigh;
    }
}

/**
 * Apply conditions on the velocity for the structure interfaces
 * 
 * No slip condition for now but can be generalized.
 * 
 * Note: very much assuming no opposing boundary conditions are possible, and no singleton
 *       boundary conditions are possible. Additionally, 1-thick interfaces will fail (which may be OK).
 * 
 * TODO: this should not be virtual, going to be the same for all codes.
*/
void NSSolver::applyInterfaceBCs() {
    this->applyFluidBCs(this->nx, this->ny, this->u, this->v);

    int i, j, xi, yi;
    objects::FSIObject obj;

    int mo = this->methodOrd;

    double hx = this->x[1] - this->x[0];
    double hy = this->y[1] - this->y[0];
    int nDir[2];

    for (j = 1; j < ny-1; j++ ) {
        for (i = 1; i < nx-1; i++) {
            xi = mo+i;
            yi = mo+j;
            obj = pool->objAtIndex(i, j);

            if (pool->isInterface(obj)) {
                // Indicate if the current cell has cells structure cells in the cardinal directions
                bool structUpDown = pool->hasStructInDir(obj, objects::NORTH) || pool->hasStructInDir(obj, objects::SOUTH);
                bool structLeftRight = pool->hasStructInDir(obj, objects::EAST) || pool->hasStructInDir(obj, objects::WEST);
                int neighbourDirY = (pool->hasStructInDir(obj, objects::NORTH))? -1 : 1;
                int neighbourDirX = (pool->hasStructInDir(obj, objects::EAST))? -1 : 1;
                
                double curPnt[2] = {simutils::midpoint(x[i], x[i+1]), simutils::midpoint(y[j], y[j+1])}; // Current cell center
                double stagPnt[2] = {0, 0}; // Staggered grid value

                pool->getNormalDir(obj, nDir);

                double velObj, velNeigh;

                if (structUpDown) {

                    if (neighbourDirY == -1) {
                        stagPnt[0] = curPnt[0];
                        stagPnt[1] = curPnt[1] - hy/2.0;

                        velObj = pool->getObjV(i, j);
                        velNeigh = v[yi-2][xi+nDir[0]];
                        
                        v[yi-1][xi] = this->getKinematicBoundaryLC(i, j, velObj, velNeigh, curPnt, stagPnt, nDir);
                    } else {
                        stagPnt[0] = curPnt[0];
                        stagPnt[1] = curPnt[1] + hy/2.0;

                        velObj = pool->getObjV(i, j);
                        velNeigh = v[yi-1][xi+nDir[0]];
                        v[yi][xi] = this->getKinematicBoundaryLC(i, j, velObj, velNeigh, curPnt, stagPnt, nDir);
                        
                    }

                    // If the cells to either west or east are not structures, need to assign at least the right
                    // u value. If the first case holds, then only assign the west u value if it is a fluid cell.
                    if (!structLeftRight) {
                        if (pool->objAtIndex(i+1, j) == objects::FLUID_C) {
                            stagPnt[0] = curPnt[0] + hx/2.0;
                            stagPnt[1] = curPnt[1];

                            velObj = pool->getObjU(i, j);
                            velNeigh = u[yi+neighbourDirY][xi+nDir[0]];

                            u[yi][xi] = this->getKinematicBoundaryLC(i, j, velObj, velNeigh, curPnt, stagPnt, nDir);
                        } else {
                            u[yi][xi] = 2.0*pool->getObjU(i, j) - u[yi+neighbourDirY][xi];
                        }

                        // If the western cell is a fluid cell, assign the u value
                        if (pool->objAtIndex(i-1, j) == objects::FLUID_C) {
                            stagPnt[0] = curPnt[0] - hx/2.0;
                            stagPnt[1] = curPnt[1];

                            velObj = pool->getObjU(i, j);
                            velNeigh = u[yi+neighbourDirY][xi-1+nDir[0]];

                            u[yi][xi-1] = this->getKinematicBoundaryLC(i, j, velObj, velNeigh, curPnt, stagPnt, nDir);
                        }
                    }
                }

                if (structLeftRight) {

                    if (neighbourDirX == -1) {
                        stagPnt[0] = curPnt[0] - hx/2.0;
                        stagPnt[1] = curPnt[1];

                        velObj = pool->getObjU(i, j);
                        velNeigh = u[yi+nDir[1]][xi-2];
                        
                        u[yi][xi-1] = this->getKinematicBoundaryLC(i, j, velObj, velNeigh, curPnt, stagPnt, nDir);
                    } else {
                        stagPnt[0] = curPnt[0] + hx/2.0;
                        stagPnt[1] = curPnt[1];

                        velObj = pool->getObjU(i, j);
                        velNeigh = u[yi+nDir[1]][xi+1];
                        u[yi][xi] = this->getKinematicBoundaryLC(i, j, velObj, velNeigh, curPnt, stagPnt, nDir);
                    }

                    // If the cells to either west or east are not structures, need to assign at least the right
                    // u value. If the first case holds, then only assign the west u value if it is a fluid cell.
                    if (!structUpDown) {
                        if (pool->objAtIndex(i, j+1) == objects::FLUID_C) {
                            stagPnt[0] = curPnt[0];
                            stagPnt[1] = curPnt[1]+hy/2.0;

                            velObj = pool->getObjV(i, j);
                            velNeigh = v[yi+nDir[1]][xi+neighbourDirX];

                            v[yi][xi] = this->getKinematicBoundaryLC(i, j, velObj, velNeigh, curPnt, stagPnt, nDir);

                        } else {
                            v[yi][xi] = 2.0*pool->getObjV(i, j) - v[yi][xi+neighbourDirX];
                        }

                        // If the western cell is a fluid cell, assign the u value
                        if (pool->objAtIndex(i, j-1) == objects::FLUID_C) {
                            stagPnt[0] = curPnt[0];
                            stagPnt[1] = curPnt[1]-hy/2.0;

                            velObj = pool->getObjV(i, j);
                            velNeigh = v[yi-1+nDir[1]][xi+neighbourDirX];

                            v[yi-1][xi] = this->getKinematicBoundaryLC(i, j, velObj, velNeigh, curPnt, stagPnt, nDir);
                        }
                    }
                }
            } 
        }
    }
}

/**
 * Method to return u
*/
double** NSSolver::getU() {
    return this->iu;
}

/**
 * Method to return v
*/
double** NSSolver::getV() {
    return this->iv;
}

/**
 * Test method to get full U
*/
double** NSSolver::testGetU() {
    return this->u;
}

/**
 * Test method to get full v
*/
double** NSSolver::testGetV() {
    return this->v;
}

/**
 * Function which decides the time step for the current step
*/
double NSSolver::getDt() {
    double umax = abs(this->u[1][1]);
    double vmax = abs(this->v[1][1]);
    int i;
    int j;
    double temp;
    double Re = params->Re;
    int mo = this->methodOrd;

    objects::FSIObject obj;

    for (j = 0; j < ny; j++) {
        for (i = 0; i < nx; i++){
            obj = pool->objAtIndex(i, j);

            if (obj == objects::FLUID_C) {
                temp = abs(this->u[mo+j][mo+i]);
                umax = (temp > umax) ? temp : umax;

                temp = abs(this->v[mo+j][mo+i]);
                vmax = (temp > vmax) ? temp : vmax;
            }
        }
    }

    return min( (Re/2.0)*(1.0/(simutils::square(1.0/this->dx) + simutils::square(1.0/this->dy))),
                    min(this->dx/umax, this->dy/vmax) );
}

void NSSolver::writeToFile(const char *fname) {
    MomentumSolver2D::writeToFile(fname);
}

void NSSolver::writePoolToFile(const char *poolName, const char *poolVelName) {
    MomentumSolver2D::writePoolToFile(poolName, poolVelName);
}

void NSSolver::outputMedialAxis(const char *fname) {
    pool->outputMedialAxisApprox(fname);
}

void NSSolver::outputTracers(const char *fname) {
    pool->outputTracers(fname);
}

void NSSolver::outputStructure(int structNum, const char *fname) {
    pool->outputStructure(structNum, fname);
}

void NSSolver::outputAllStructures(const char *fname) {
    pool->outputAllStructures(fname);
}

void NSSolver::outputStructureNodes(int structNum, const char *fname) {
    pool->outputStructureNodes(structNum, fname);
}

void NSSolver::outputAllStructureNodes(const char *fname) {
    pool->outputAllStructureNodes(fname);
}

void NSSolver::outputStructureVels(int structNum, const char *fname) {
    pool->outputStructureVels(structNum, fname);
}

void NSSolver::outputAllStructureVels(const char *fname) {
    pool->outputAllStructureVels(fname);
}