#include "MomentumSolver3D.h"
#include "Boundary3D.h"
#include "SolidObject3D.h"
#include "NSSolver3D.h"
#include "../Utils/SimUtilities.h"
#include "../Utils/Discretizations.h"
#include <cstdlib>
#include <cassert>
#include "SimParams3D.h"

using namespace std;

/**
 * Class constructor. In this case we can simply call the parent class.
*/
NSSolver3D::NSSolver3D(Boundary3D &boundary, 
                    vector<SolidObject3D> &solidObjects,
                    SimParams3D &params,
                    std::function<void (int,int,int,int,double*,double*,double*,double***,double***,double***)> initialConditions,
                    void (*boundaryConditions)(int,int,int,double***,double***,double***))
                    : MomentumSolver3D(boundary, solidObjects, params, initialConditions, boundaryConditions)
{
}

/**
 * Function to take a time step. Simply calls the step method in parent class
 * MomentumSolver3D
*/
double NSSolver3D::step(double tEnd, double safetyFactor) {
    return MomentumSolver3D::step(tEnd, safetyFactor);
}

/**
 * Function which takes explicit time step for the velocity terms in the
 * momentum equations.
 * Note: assuming no body forces
*/
void NSSolver3D::updateF(Pool3D *pool) {
    int i, j, k;
    int xi, yi, zi;

    // Extract the Reynolds number
    double Re = params->Re;

    int mo = methodOrd;

    double laplacian;
    double convective;
    for (k = 0; k < this->nz; k++) {
        for (j = 0; j < this->ny; j++) {
            for (i = 0; i < this->nx-1; i++) {
                xi = i + mo;
                yi = j + mo;
                zi = k + mo;

                if (pool->objAtIndex(i, j, k) == objects::FLUID_C && pool->isUpdateableU(i, j, k)) {
                    laplacian = (discs::firstOrder3D_lap_uxx(xi, yi, zi, this->dx, this->u) 
                                + discs::firstOrder3D_lap_uyy(xi, yi, zi, this->dy, this->u)
                                + discs::firstOrder3D_lap_uzz(xi, yi, zi, this->dz, this->u));

                    if (params->useEno) {
                        double gam = max(abs((u[zi][yi][xi]*this->dt)/(this->dx)), 
                            max(abs(v[zi][yi][xi]*this->dt)/(this->dy), abs(w[zi][yi][xi]*this->dt)/(this->dz)));
                        convective = discs::firstOrder3D_conv_usqx_stab(xi, yi, zi, this->dx, this->u, gam)
                                + discs::firstOrder3D_conv_uvy_stab(xi, yi, zi, this->dy, this->u, this->v, gam)
                                + discs::firstOrder3D_conv_uwz_stab(xi, yi, zi, this->dz, u, w, gam);
                    } else {
                        convective = discs::firstOrder3D_conv_usqx(xi, yi, zi, this->dx, this->u)
                                + discs::firstOrder3D_conv_uvy(xi, yi, zi, this->dy, this->u, this->v)
                                + discs::firstOrder3D_conv_uwz(xi, yi, zi, this->dz, u, w);
                    }
                    this->FU[zi][yi][xi] = this->u[zi][yi][xi] + this->dt*( (1.0/Re)*laplacian - convective + gx );
                    // this->FU[zi][yi][xi] = this->u[zi][yi][xi] + this->dt*( (1.0/Re)*laplacian - convective);
                }
            }
        }
    }

    for (k = 0; k < this->nz; k++) {
        for (j = 0; j < this->ny-1; j++) {
            for (i = 0; i < this->nx; i++) {
                xi = i + mo;
                yi = j + mo;
                zi = k + mo;

                if (pool->objAtIndex(i, j, k) == objects::FLUID_C && pool->isUpdateableV(i, j, k)) {
                    laplacian = discs::firstOrder3D_lap_vxx(xi, yi, zi, this->dx, this->v)
                                + discs::firstOrder3D_lap_vyy(xi, yi, zi, this->dy, this->v)
                                + discs::firstOrder3D_lap_vzz(xi, yi, zi, this->dz, this->v);
                    if (params->useEno) {
                        double gam = max(abs((u[zi][yi][xi]*this->dt)/(this->dx)), 
                            max(abs(v[zi][yi][xi]*this->dt)/(this->dy), abs(w[zi][yi][xi]*this->dt)/(this->dz)));
                        convective = discs::firstOrder3D_conv_vux_stab(xi, yi, zi, this->dx, this->v, this->u, gam)
                                + discs::firstOrder3D_conv_vsqy_stab(xi, yi, zi, this->dy, this->v, gam)
                                + discs::firstOrder3D_conv_vwz_stab(xi, yi, zi, this->dz, this->v, this->w, gam);
                    } else {
                        convective = discs::firstOrder3D_conv_vux(xi, yi, zi, this->dx, this->v, this->u)
                                + discs::firstOrder3D_conv_vsqy(xi, yi, zi, this->dy, this->v)
                                + discs::firstOrder3D_conv_vwz(xi, yi, zi, this->dz, this->v, this->w);
                    }
                    this->FV[zi][yi][xi] = this->v[zi][yi][xi] + this->dt*( (1.0/Re)*laplacian - convective + gy);
                }
            }
        }
    }

    for (k = 0; k < this->nz-1; k++) {
        for (j = 0; j < this->ny; j++) {
            for (i = 0; i < this->nx; i++) {
                xi = i + mo;
                yi = j + mo;
                zi = k + mo;

                if (pool->objAtIndex(i, j, k) == objects::FLUID_C && pool->isUpdateableW(i, j, k)) {
                    laplacian = discs::firstOrder3D_lap_wxx(xi, yi, zi, this->dx, this->w)
                                + discs::firstOrder3D_lap_wyy(xi, yi, zi, this->dy, this->w)
                                + discs::firstOrder3D_lap_wzz(xi, yi, zi, this->dz, this->w);
                    if (params->useEno) {
                        double gam = max(abs((u[zi][yi][xi]*this->dt)/(this->dx)), 
                            max(abs(v[zi][yi][xi]*this->dt)/(this->dy), abs(w[zi][yi][xi]*this->dt)/(this->dz)));
                        convective = discs::firstOrder3D_conv_wux_stab(xi, yi, zi, this->dx, this->w, this->u, gam)
                                + discs::firstOrder3D_conv_wvy_stab(xi, yi, zi, this->dy, this->w, this->v, gam)
                                + discs::firstOrder3D_conv_wsqz_stab(xi, yi, zi, this->dz, this->w, gam);
                    } else {
                        convective = discs::firstOrder3D_conv_wux(xi, yi, zi, this->dx, this->w, this->u)
                                + discs::firstOrder3D_conv_wvy(xi, yi, zi, this->dy, this->w, this->v)
                                + discs::firstOrder3D_conv_wsqz(xi, yi, zi, this->dz, this->w);

                    }
                    this->FW[zi][yi][xi] = this->w[zi][yi][xi] + this->dt*( (1.0/Re)*laplacian - convective + gz);
                    // cout << "making FW = " << this->FW[zi][yi][xi] << endl;
                    // cout << "w = " << this->w[zi][yi][xi] << endl;
                }
            }
        }
    }
    // assert(false);

    // Set the boundary values for F required in the pressure solver.
    // so that it can be most consistent. TODO: find a way more efficient way of
    // implementing these checks (ofc after we have a prototype)

    // The boundary conditions on the interface
    objects::FSIObject obj;

    for (k = 0; k < nz; k++) {
        for (j = 0; j < ny; j++) {
            for (i = 0; i < nx; i++) {
                xi = methodOrd+i;
                yi = methodOrd+j;
                zi = methodOrd+k;
                obj = pool->objAtIndex(i, j, k);
                if (pool->isInterface(obj)) {
                    if (pool->hasStructInDir(obj, objects::NORTH)) {
                        FV[zi][yi-1][xi] = v[zi][yi-1][xi];
                    }

                    if (pool->hasStructInDir(obj, objects::SOUTH)) {
                        FV[zi][yi][xi] = v[zi][yi][xi];
                    }

                    if (pool->hasStructInDir(obj, objects::EAST)) {
                        FU[zi][yi][xi-1] = u[zi][yi][xi-1];
                    }

                    if (pool->hasStructInDir(obj, objects::WEST)) {
                        FU[zi][yi][xi] = u[zi][yi][xi];
                    }

                    if (pool->hasStructInDir(obj, objects::UP)) {
                        FW[zi-1][yi][xi] = w[zi-1][yi][xi];
                    }

                    if (pool->hasStructInDir(obj, objects::DOWN)) {
                        FW[zi][yi][xi] = w[zi][yi][xi];
                    }
                }
            }
        }
    }


    // The boundary values on the bounding rectangle of the fluid domain
    for (j = 1; j <= this->ny; j++) {
        for (i = 1; i <= this->nx; i++) {
            this->FW[0][j][i] = this->w[0][j][i];
            this->FW[this->nz][j][i] = this->w[this->nz][j][i];
        }
    }

    for (k = 1; k <= this->nz; k++) {
        for (j = 1; j <= this->ny; j++) {
            this->FU[k][j][0] = this->u[k][j][0];
            this->FU[k][j][this->nx] = this->u[k][j][this->nx];
        }

        for (i = 1; i <= this->nx; i++) {
            this->FV[k][0][i] = this->v[k][0][i];
            this->FV[k][this->ny][i] = this->v[k][this->ny][i];
        }
    }
}

/**
 * Compute the no-penetration boundary condition using trilinear interpolation.
*/
double NSSolver3D::getKinematicBoundaryLC(int i, int j, int k, double velObj, double velNeigh,
                                        double cPnt[3], double stagPnt[3], int nDir[3]) {
    double hx = this->x[1] - this->x[0];
    double hy = this->y[1] - this->y[0];
    double hz = this->z[1] - this->z[0];
    double nPnt[3];
    double inPnt[3];
    double distVec[3];
    double gamPnt[3];
    double phi_out, phi_in, dist;
    double t;

    // Take the unit outward normal
    double outN[3] = {(double) nDir[0], (double) nDir[1], (double) nDir[2]};
    simutils::normalize3D(outN);

    nPnt[0] = stagPnt[0] + hx*(double)nDir[0];
    nPnt[1] = stagPnt[1] + hy*(double)nDir[1];
    nPnt[2] = stagPnt[2] + hz*(double)nDir[2];

    inPnt[0] = stagPnt[0] - hx*(double)nDir[0];
    inPnt[1] = stagPnt[1] - hy*(double)nDir[1];
    inPnt[2] = stagPnt[2] - hz*(double)nDir[2];

    phi_out = pool->interpolatePhi(stagPnt[0], stagPnt[1], stagPnt[2]);
    phi_in = pool->interpolatePhi(inPnt[0], inPnt[1], inPnt[2]);

    if (phi_out > 0 && phi_in > 0) {
        // return fluid velocity
        return velNeigh;
    } else if (phi_out < 0 && phi_in < 0){
        return velObj;
    } else {
        // Compute distance to interface (from internal)
        distVec[0] = inPnt[0] - stagPnt[0];
        distVec[1] = inPnt[1] - stagPnt[1];
        distVec[2] = inPnt[2] - stagPnt[2];
        dist = abs((phi_in/(phi_in - phi_out))*simutils::eucNorm3D(distVec));

        // Compute location of interface
        gamPnt[0] = inPnt[0] + dist*outN[0];
        gamPnt[1] = inPnt[1] + dist*outN[1];
        gamPnt[2] = inPnt[2] + dist*outN[2];

        // Compute relative val
        if (nDir[0] != 0) {
            t = (gamPnt[0] - inPnt[0])/(nPnt[0] - inPnt[0]);
        } else if (nDir[1] != 0) {
            t = (gamPnt[1] - inPnt[1])/(nPnt[1] - inPnt[1]);
        } else {
            t = (gamPnt[2] - inPnt[2])/(nPnt[2] - inPnt[2]);
        }

        assert(t>=0 && t <= 1);

        return t*velObj + (1-t)*velNeigh;
    }
}

/**
 * Apply conditions at the fluid-structure interface
 * 
*/
void NSSolver3D::applyInterfaceBCs() {
    this->resetBoundaryConditions();
    this->applyFluidBCs(this->nx,this->ny, this->nz, this->u, this->v, this->w);

    int i, j, k, xi, yi, zi;
    objects::FSIObject obj;

    int mo = this->methodOrd;

    double hx = this->x[1] - this->x[0];
    double hy = this->y[1] - this->y[0];
    double hz = this->z[1] - this->z[0];

    int nDir[3];

    for (k = 1; k < nz-1; k++) {
        for (j = 1; j < ny-1; j++ ) {
            for (i = 1; i < nx-1; i++) {
                xi = mo+i;
                yi = mo+j;
                zi = mo+k;
                obj = pool->objAtIndex(i, j, k);

                if (pool->isInterface(obj)) {
                    // Indicate if the current cell has cells structure cells in the cardinal directions
                    bool structY = pool->hasStructInDir(obj, objects::NORTH) || pool->hasStructInDir(obj, objects::SOUTH);
                    bool structX = pool->hasStructInDir(obj, objects::EAST) || pool->hasStructInDir(obj, objects::WEST);
                    bool structZ = pool->hasStructInDir(obj, objects::UP) || pool->hasStructInDir(obj, objects::DOWN);

                    int neighbourDirY = (pool->hasStructInDir(obj, objects::NORTH)) ? -1 : 1;
                    int neighbourDirX = (pool->hasStructInDir(obj, objects::EAST)) ? -1 : 1;
                    int neighbourDirZ = (pool->hasStructInDir(obj, objects::UP)) ? -1 : 1;
                    
                    double curPnt[3] = {simutils::midpoint(x[i], x[i+1]),
                                            simutils::midpoint(y[j], y[j+1]),
                                            simutils::midpoint(z[k], z[k+1])}; // Current cell center
                    double stagPnt[3] = {0, 0, 0}; // Staggered grid value

                    pool->getNormalDir(obj, nDir);

                    double velObj, velNeigh;

                    // Handle the y-direction edge variables
                    if (structY) {

                        if (neighbourDirY == -1) {
                            stagPnt[0] = curPnt[0];
                            stagPnt[1] = curPnt[1] - hy/2.0;
                            stagPnt[2] = curPnt[2];

                            velObj = pool->getObjV(i, j, k);
                            velNeigh = v[zi+nDir[2]][yi-2][xi+nDir[0]];
                            
                            v[zi][yi-1][xi] = this->getKinematicBoundaryLC(i, j, k, velObj, velNeigh, curPnt, stagPnt, nDir);


                        } else {
                            stagPnt[0] = curPnt[0];
                            stagPnt[1] = curPnt[1] + hy/2.0;
                            stagPnt[2] = curPnt[2];

                            velObj = pool->getObjV(i, j, k);
                            velNeigh = v[zi+nDir[2]][yi+1][xi+nDir[0]];
                            v[zi][yi][xi] = this->getKinematicBoundaryLC(i, j, k, velObj, velNeigh, curPnt, stagPnt, nDir);
                            
                        }

                        // If the cells to either west or east are not structures, need to assign at least the right
                        // u value. If the first case holds, then only assign the west u value if it is a fluid cell.
                        if (!structX) {
                            if (pool->objAtIndex(i+1, j, k) == objects::FLUID_C) {
                                stagPnt[0] = curPnt[0] + hx/2.0;
                                stagPnt[1] = curPnt[1];
                                stagPnt[2] = curPnt[2];

                                velObj = pool->getObjU(i, j, k);
                                velNeigh = u[zi+nDir[2]][yi+neighbourDirY][xi+nDir[0]];

                                u[zi][yi][xi] = this->getKinematicBoundaryLC(i, j, k, velObj, velNeigh, curPnt, stagPnt, nDir);
                            } else {
                                u[zi][yi][xi] = 2.0*pool->getObjU(i, j, k) - u[zi][yi+neighbourDirY][xi];
                            }

                            // If the western cell is a fluid cell, assign the u value
                            if (pool->objAtIndex(i-1, j, k) == objects::FLUID_C) {
                                stagPnt[0] = curPnt[0] - hx/2.0;
                                stagPnt[1] = curPnt[1];
                                stagPnt[2] = curPnt[2];

                                velObj = pool->getObjU(i, j, k);
                                velNeigh = u[zi+nDir[2]][yi+neighbourDirY][xi-1+nDir[0]];

                                u[zi][yi][xi-1] = this->getKinematicBoundaryLC(i, j, k, velObj, velNeigh, curPnt, stagPnt, nDir);
                            }
                        }

                        if (!structZ) {
                            if (pool->objAtIndex(i, j, k+1) == objects::FLUID_C) {
                                stagPnt[0] = curPnt[0];
                                stagPnt[1] = curPnt[1];
                                stagPnt[2] = curPnt[2] + hz/2.0;

                                velObj = pool->getObjW(i, j, k);
                                velNeigh = w[zi+nDir[2]][yi+neighbourDirY][xi+nDir[0]];

                                w[zi][yi][xi] = this->getKinematicBoundaryLC(i, j, k, velObj, velNeigh, curPnt, stagPnt, nDir);
                            } else {
                                w[zi][yi][xi] = 2.0*pool->getObjW(i, j, k) - w[zi][yi+neighbourDirY][xi];
                            }

                            if (pool->objAtIndex(i, j, k-1) == objects::FLUID_C) {
                                stagPnt[0] = curPnt[0];
                                stagPnt[1] = curPnt[1];
                                stagPnt[2] = curPnt[2] - hz/2.0;;

                                velObj = pool->getObjW(i, j, k);
                                velNeigh = w[zi+nDir[2]][yi+neighbourDirY][xi-1+nDir[0]];

                                w[zi-1][yi][xi] = this->getKinematicBoundaryLC(i, j, k, velObj, velNeigh, curPnt, stagPnt, nDir);
                            }
                        }
                    }

                    // Handle the x-direction edge variables
                    if (structX) {

                        if (neighbourDirX == -1) {
                            stagPnt[0] = curPnt[0] - hx/2.0;
                            stagPnt[1] = curPnt[1];
                            stagPnt[2] = curPnt[2];

                            velObj = pool->getObjU(i, j, k);
                            velNeigh = u[zi+nDir[2]][yi+nDir[1]][xi-2];
                            
                            u[zi][yi][xi-1] = this->getKinematicBoundaryLC(i, j, k, velObj, velNeigh, curPnt, stagPnt, nDir);
                        } else {
                            stagPnt[0] = curPnt[0] + hx/2.0;
                            stagPnt[1] = curPnt[1];
                            stagPnt[2] = curPnt[2];

                            velObj = pool->getObjU(i, j, k);
                            velNeigh = u[zi+nDir[2]][yi+nDir[1]][xi+1];
                            u[zi][yi][xi] = this->getKinematicBoundaryLC(i, j, k, velObj, velNeigh, curPnt, stagPnt, nDir);
                        }

                        // If the cells to either west or east are not structures, need to assign at least the right
                        // u value. If the first case holds, then only assign the west u value if it is a fluid cell.
                        if (!structY) {
                            if (pool->objAtIndex(i, j+1, k) == objects::FLUID_C) {
                                stagPnt[0] = curPnt[0];
                                stagPnt[1] = curPnt[1]+hy/2.0;
                                stagPnt[2] = curPnt[2];

                                velObj = pool->getObjV(i, j, k);
                                velNeigh = v[zi+nDir[2]][yi+nDir[1]][xi+neighbourDirX];

                                v[zi][yi][xi] = this->getKinematicBoundaryLC(i, j, k, velObj, velNeigh, curPnt, stagPnt, nDir);

                            } else {
                                v[zi][yi][xi] = 2.0*pool->getObjV(i, j, k) - v[zi][yi][xi+neighbourDirX];
                            }

                            // If the western cell is a fluid cell, assign the u value
                            if (pool->objAtIndex(i, j-1, k) == objects::FLUID_C) {
                                stagPnt[0] = curPnt[0];
                                stagPnt[1] = curPnt[1]-hy/2.0;
                                stagPnt[2] = curPnt[2];

                                velObj = pool->getObjV(i, j, k);
                                velNeigh = v[zi+nDir[2]][yi-1+nDir[1]][xi+neighbourDirX];

                                v[zi][yi-1][xi] = this->getKinematicBoundaryLC(i, j, k, velObj, velNeigh, curPnt, stagPnt, nDir);
                            }
                        }

                        if (!structZ) {
                            if (pool->objAtIndex(i, j, k+1) == objects::FLUID_C) {
                                stagPnt[0] = curPnt[0];
                                stagPnt[1] = curPnt[1];
                                stagPnt[2] = curPnt[2]+hz/2.0;

                                velObj = pool->getObjW(i, j, k);
                                velNeigh = w[zi+nDir[2]][yi+nDir[1]][xi+neighbourDirX];

                                w[zi][yi][xi] = this->getKinematicBoundaryLC(i, j, k, velObj, velNeigh, curPnt, stagPnt, nDir);

                            } else {
                                w[zi][yi][xi] = 2.0*pool->getObjW(i, j, k) - w[zi][yi][xi+neighbourDirX];
                            }

                            // If the western cell is a fluid cell, assign the u value
                            if (pool->objAtIndex(i, j, k-1) == objects::FLUID_C) {
                                stagPnt[0] = curPnt[0];
                                stagPnt[1] = curPnt[1];
                                stagPnt[2] = curPnt[2]-hz/2.0;

                                velObj = pool->getObjW(i, j, k);
                                velNeigh = w[zi+nDir[2]][yi-1+nDir[1]][xi+neighbourDirX];

                                w[zi][yi-1][xi] = this->getKinematicBoundaryLC(i, j, k, velObj, velNeigh, curPnt, stagPnt, nDir);
                            }
                        }
                    }

                    // Handle the x-direction edge variables
                    if (structZ) {

                        if (neighbourDirZ == -1) {
                            stagPnt[0] = curPnt[0];
                            stagPnt[1] = curPnt[1];
                            stagPnt[2] = curPnt[2] - hz/2.0;

                            velObj = pool->getObjW(i, j, k);
                            velNeigh = w[zi-2][yi+nDir[1]][xi+nDir[0]];
                            
                            w[zi-1][yi][xi] = this->getKinematicBoundaryLC(i, j, k, velObj, velNeigh, curPnt, stagPnt, nDir);
                        } else {
                            stagPnt[0] = curPnt[0];
                            stagPnt[1] = curPnt[1];
                            stagPnt[2] = curPnt[2] + hz/2.0;

                            velObj = pool->getObjW(i, j, k);
                            velNeigh = w[zi+1][yi+nDir[1]][xi+nDir[0]];
                            w[zi][yi][xi] = this->getKinematicBoundaryLC(i, j, k, velObj, velNeigh, curPnt, stagPnt, nDir);
                        }

                        // If the cells to either west or east are not structures, need to assign at least the right
                        // u value. If the first case holds, then only assign the west u value if it is a fluid cell.
                        if (!structX) {
                            if (pool->objAtIndex(i+1, j, k) == objects::FLUID_C) {
                                stagPnt[0] = curPnt[0] + hx/2.0;
                                stagPnt[1] = curPnt[1];
                                stagPnt[2] = curPnt[2];

                                velObj = pool->getObjU(i, j, k);
                                velNeigh = u[zi+neighbourDirZ][yi+nDir[1]][xi+nDir[0]];

                                u[zi][yi][xi] = this->getKinematicBoundaryLC(i, j, k, velObj, velNeigh, curPnt, stagPnt, nDir);
                            } else {
                                u[zi][yi][xi] = 2.0*pool->getObjU(i, j, k) - u[zi+neighbourDirZ][yi+nDir[1]][xi+nDir[0]];
                            }

                            // If the western cell is a fluid cell, assign the u value
                            if (pool->objAtIndex(i-1, j, k) == objects::FLUID_C) {
                                stagPnt[0] = curPnt[0] - hx/2.0;
                                stagPnt[1] = curPnt[1];
                                stagPnt[2] = curPnt[2];

                                velObj = pool->getObjU(i, j, k);
                                velNeigh = u[zi+neighbourDirZ][yi+nDir[1]][xi-1+nDir[0]];

                                u[zi][yi][xi-1] = this->getKinematicBoundaryLC(i, j, k, velObj, velNeigh, curPnt, stagPnt, nDir);
                            }
                        }

                        if (!structY) {
                            if (pool->objAtIndex(i, j+1, k) == objects::FLUID_C) {
                                stagPnt[0] = curPnt[0];
                                stagPnt[1] = curPnt[1] + hy/2.0;
                                stagPnt[2] = curPnt[2];

                                velObj = pool->getObjV(i, j, k);
                                velNeigh = v[zi+neighbourDirZ][yi+nDir[1]][xi+nDir[0]];

                                v[zi][yi][xi] = this->getKinematicBoundaryLC(i, j, k, velObj, velNeigh, curPnt, stagPnt, nDir);

                            } else {
                                v[zi][yi][xi] = 2.0*pool->getObjV(i, j, k) - v[zi][yi][xi+neighbourDirX];
                            }

                            // If the western cell is a fluid cell, assign the u value
                            if (pool->objAtIndex(i, j-1, k) == objects::FLUID_C) {
                                stagPnt[0] = curPnt[0];
                                stagPnt[1] = curPnt[1]-hy/2.0;
                                stagPnt[2] = curPnt[2];

                                velObj = pool->getObjV(i, j, k);
                                velNeigh = v[zi+neighbourDirZ][yi-1+nDir[1]][xi+nDir[0]];

                                v[zi][yi-1][xi] = this->getKinematicBoundaryLC(i, j, k, velObj, velNeigh, curPnt, stagPnt, nDir);
                            }
                        }
                    }
                } 
            }
        }
    }
}

/**
 * Apply the boundary conditions
*/
void NSSolver3D::resetBoundaryConditions() {
    int i, j, k;

    // Set all the boundary conditions to 0.
    for (j = 1; j <= this->ny; j++) {
        for (i = 1; i <= this->nx; i++) {
            this->w[0][j][i] = 0.0;
            this->w[this->nz][j][i] = 0.0;

            this->u[0][j][i] = -this->u[1][j][i];
            this->u[this->nz+1][j][i] = -this->u[this->nz][j][i];

            this->v[0][j][i] = -this->v[1][j][i];
            this->v[this->nz+1][j][i] = -this->v[this->nz][j][i];
        }
    }

    for (k = 1; k < nz+1; k++) {
        for (j = 1; j <= this->ny; j++) {
            this->u[k][j][0] = 0.0;
            this->u[k][j][this->nx] = 0.0;

            this->v[k][j][0] = -this->v[k][j][1];
            this->w[k][j][0] = -this->w[k][j][1];

            this->v[k][j][this->nx+1] = -this->v[k][j][this->nx];
            this->w[k][j][this->nx+1] = -this->w[k][j][this->nx];
        }

        for (i = 1; i <= this->nx; i++) {
            this->u[k][0][i] = -this->u[k][1][i];
            this->u[k][this->ny+1][i] = -this->u[k][this->ny][i];

            this->v[k][0][i] = 0.0;
            this->v[k][this->ny][i] = 0.0;

            this->w[k][0][i] = -this->w[k][1][i];
            this->w[k][this->ny+1][i] = -this->w[k][this->ny][i];
        }
    }
}

/**
 * Method to return u
*/
double*** NSSolver3D::getU() {
    return this->iu;
}

/**
 * Method to return v
*/
double*** NSSolver3D::getV() {
    return this->iv;
}

/**
 * Method to return v
*/
double*** NSSolver3D::getW() {
    return this->iw;
}

/**
 * Test method to get full U
*/
double*** NSSolver3D::testGetU() {
    return this->u;
}

/**
 * Test method to get full v
*/
double*** NSSolver3D::testGetV() {
    return this->v;
}

/**
 * Test method to get full v
*/
double*** NSSolver3D::testGetW() {
    return this->w;
}

/**
 * Function which decides the time step for the current step
*/
double NSSolver3D::getDt() {
    double umax = abs(this->u[1][1][1]);
    double vmax = abs(this->v[1][1][1]);
    double wmax = abs(this->w[1][1][1]);
    int i, j, k;
    double temp;
    double Re = params->Re;
    int mo = this->methodOrd;

    objects::FSIObject obj;

    for (k = 0; k < nz; k++) {
        for (j = 0; j < ny; j++) {
            for (i = 0; i < nx; i++){
                obj = pool->objAtIndex(i, j, k);

                if (obj == objects::FLUID_C) {
                    temp = abs(this->u[mo+k][mo+j][mo+i]);
                    umax = (temp > umax) ? temp : umax;

                    temp = abs(this->v[mo+k][mo+j][mo+i]);
                    vmax = (temp > vmax) ? temp : vmax;

                    temp = abs(this->w[mo+k][mo+j][mo+i]);
                    wmax = (temp > wmax) ? temp : wmax;
                }
            }
        }

    }

    return min( (Re/2.0)*(
                    1.0/(simutils::square(1.0/this->dx)
                        + simutils::square(1.0/this->dy)
                        + simutils::square(1.0/this->dz))),
                    min(min(this->dx/umax, this->dy/vmax), this->dz/wmax) );
}

void NSSolver3D::writeToFile(const char *fname) {
    MomentumSolver3D::writeToFile(fname);
}

void NSSolver3D::writePoolToFile(const char *poolName, const char *poolVelName) {
    MomentumSolver3D::writePoolToFile(poolName, poolVelName);
}

void NSSolver3D::outputStructure(int structNum, const char *fname) {
    pool->outputStructure(structNum, fname);
}

void NSSolver3D::outputAllStructures(const char *fname) {
    pool->outputAllStructures(fname);
}

void NSSolver3D::outputStructureNodes(int structNum, const char *fname) {
    pool->outputStructureNodes(structNum, fname);
}

void NSSolver3D::outputAllStructureNodes(const char *fname) {
    pool->outputAllStructureNodes(fname);
}

void NSSolver3D::outputStructureVels(int structNum, const char *fname) {
    pool->outputStructureVels(structNum, fname);
}

void NSSolver3D::outputAllStructureVels(const char *fname) {
    pool->outputAllStructureVels(fname);
}

void NSSolver3D::outputTracers(const char *fname) {
    pool->outputTracers(fname);
}

void NSSolver3D::outputMedialAxis(const char *fname) {
    pool->outputMedialAxis(fname);
}