#include "Boundary.h"
#include "SolidObject.h"
#include "NSSolver.h"
#include "../Utils/SimUtilities.h"
#include "../Utils/Discretizations.h"
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <unordered_map>
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
{
    // Check that the parameters object is set
    if (!params.checkParams()) {
        cout << "ERROR: 2D Sim params not set correctly!" << endl;
        assert(false);
    }

    // Copy basic info
    this->nx = params.nx;
    this->ny = params.nx;
    this->methodOrd = params.methodOrd;

    // Create the uniform 1D meshes using the boundary object
    this->x = boundary.generateXMesh(this->nx);
    this->y = boundary.generateYMesh(this->ny);

    // The mesh spacing
    this->dx = this->x[1] - this->x[0];
    this->dy = this->y[1] - this->y[0];

    // Copy the input parameters using the default copy ctor
    this->params = &params;

    // Create the output arrays with required ghost cells for the staggered grid
    // configuration.
    // TODO: make sure that this is correct
    int hg = this->ny + 2*this->methodOrd;
    int wg = this->nx + 2*this->methodOrd;

    this->u = simutils::new_constant(hg, wg-1, 0.0);
    this->iu = simutils::new_constant(this->ny, this->nx, 0.0);
    this->FU = simutils::new_constant(hg, wg-1, 0.0);
    this->v = simutils::new_constant(hg-1, wg, 0.0);
    this->iv = simutils::new_constant(this->ny, this->nx, 0.0);
    this->FV = simutils::new_constant(hg-1, wg, 0.0);
    this->p = simutils::new_constant(hg, wg, 0.0);

    // Create pressure solver and its parameters. Parameters are default for now
    // TODO: parameters should be passed as an argument.
    this->pressureSolver = new PressureSolver(this->nx, this->x, this->ny, this->y);
    this->pressureParams = new ParamIter();

    // Set the default params
    const int ILU_LEVEL = 5;

    (this->pressureParams)->order = 0; // = 0 natural;
                    //  =1, RCM

    (this->pressureParams)->level = ILU_LEVEL; // level of ilu

    (this->pressureParams)->drop_ilu = 0; // = 0 level ilu
                        // = 1 drop tol

    (this->pressureParams)->iscal = 0; // = 0 no scaling by inverse of diag
                    // = 1 scale by inverse of diag

    (this->pressureParams)->nitmax = 10000; // max number of iterations

    (this->pressureParams)->ipiv = 0; // 0 no pivoting
                                    // 1 pivoting
                                    // only if drop_ilu = 1

    (this->pressureParams)->resid_reduc = 1.e-6; // residual reduction toleranc
                                // iterations stop if any of the following true
                                //   l^2 residual reduced by ctol
                                //   number of itns = nitmax
                                //   update in each variable x[i]<toler[i]

    (this->pressureParams)->info = 0; // = 0 no iteration info written to "output"
                    // = 1 write iteration info

    (this->pressureParams)->drop_tol = 1.e-3; // drop tolerance used if drop_ilu = 1

    (this->pressureParams)->new_rhat = 0; // = 0 use r^0 as rhat
                        // = 1 use (LU)^{-1^0 for rhat

    // (this->pressureParams)->iaccel = -1; // = 0 cgstab
    (this->pressureParams)->iaccel = 0; // = 0 cgstab
                        // = 1 orthomin
                        // = -1 conj gradiend (SPD only!)

    (this->pressureParams)->north = 10; // number of orthogs for orthomin

    // Use initial condition callback function to assign u, v
    // TODO: make initial conditions more general and perhaps more
    //       able to handle irregular boundaries.
    initialConditions(this->nx, this->ny, this->methodOrd, this->x, this->y, this->u, this->v);

    applyFluidBCs = boundaryConditions;

    // Create the pool object for this solver
    this->nStructs = solidObjects.size();
    this->pool = new Pool2D(boundary, solidObjects, params);

    // Pre-compute the Barycentric weights
    // TODO: ensure that this is correct for the FV methods used later
    this->baryWeights = new double[this->methodOrd+1];
    simutils::barycentricInterp(this->methodOrd, this->baryWeights);

    this->stepTaken = false;
}

/**
 * Set the boundary conditions within the structures to ensure things are OK on the next time step
*/
void NSSolver::test_setInternal() {
    int i, j;
    int mo = methodOrd;
    for (j = 0; j < this->ny; j++) {
        for (i = 0; i < this->nx-1; i++) {
            int xi = i + mo;
            int yi = j + mo;

            // if (pool->isInterface(pool->objAtIndex(i, j))) {
            if (pool->objAtIndex(i, j)!=objects::FLUID_C) {
                this->u[yi][xi] = (this->pool)->getObjU(i, j);
                this->u[yi][xi-1] = (this->pool)->getObjU(i, j);
            }
        }
    }

    for (j = 0; j < this->ny-1; j++) {
        for (i = 0; i < this->nx; i++) {
            int xi = i + mo;
            int yi = j + mo;

            // if (pool->isInterface(pool->objAtIndex(i, j))) {
            if (pool->objAtIndex(i, j)!=objects::FLUID_C) {
                this->v[yi][xi] = (this->pool)->getObjV(i, j);
                this->v[yi-1][xi] = (this->pool)->getObjV(i, j);
            }
        }
    }

}

/**
 * Method to update the Pressure equation using a second order method.
 * The linear system is solved using ILU preconditioned PCG.
*/
void NSSolver::updateP() {
    (this->pressureSolver)->solvePCG(this->dt, this->pool, this->FU, this->FV,
                                     this->p, true, *(this->pressureParams));
}

/**
 * Method to update U using F and P. Note only the internal points are updated
 * as we only ever apply a relatively simple Dirichlet condition at the boundaries
 * for the fluids.
 * 
 * // TODO: generalize how the pressure derivatives are being evaluated to at least second order.
 */
void NSSolver::updateU() {
    int i, j;
    int mo = this->methodOrd;

    // Update u
    for (j = 0; j < this->ny; j++) {
        for (i = 0; i < this->nx-1; i++) {
            int xi = i + mo;
            int yi = j + mo;
            if (pool->objAtIndex(i, j) == objects::FLUID_C && pool->isUpdateableU(i, j)) {
                this->u[yi][xi] = this->FU[yi][xi] - (this->dt/this->dx)*(this->p[yi][xi+1] - this->p[yi][xi]);
            }
        }
    }

    // Update v
    for (j = 0; j < this->ny-1; j++) {
        for (i = 0; i < this->nx; i++) {
            int xi = i + mo;
            int yi = j + mo;
            if (pool->objAtIndex(i, j) == objects::FLUID_C && pool->isUpdateableV(i, j)) {
                this->v[yi][xi] = this->FV[yi][xi] - (this->dt/this->dy)*(this->p[yi+1][xi] - this->p[yi][xi]);
            }
        }
    }
}

/**
 * Method which interpolate fluid velocities to the cell centers.
 * 
 * TODO: make sure that this implementation makes sense
*/
void NSSolver::interpolateVelocities() {
    int i, j;
    double x, y;

    double uvals[2];
    double vvals[2];

    // Interpolate u
    for (j = 0; j < this->ny; j++) {
        y = simutils::midpoint(this->y[j], this->y[j+1]);
        for (i = 0; i < this->nx; i++) {
            x = simutils::midpoint(this->x[i], this->x[i+1]);

            // The values being interpolated
            uvals[0] = this->u[this->methodOrd+j][i];
            uvals[1] = this->u[this->methodOrd+j][i+1];

            vvals[0] = this->v[j][this->methodOrd+i];
            vvals[1] = this->v[j+1][this->methodOrd+i];

            // Interpolate the value at the cell centers
            this->iu[j][i] = simutils::barycentricInterv(x, this->x[i], this->x[i+1],
                                this->methodOrd, this->baryWeights, uvals);
            this->iv[j][i] = simutils::barycentricInterv(y, this->y[j], this->y[j+1],
                                this->methodOrd, this->baryWeights, vvals);
        }
    }
}


/**
 * Function to take a time step. General implementation of the explicit-implicit scheme.
*/
double NSSolver::step(double tEnd, double safetyFactor) {
    // Compute the time step
    // this->dt = safetyFactor*this->getDt();
    if (params->dtFixSet) {
        cout << "SETTING POOL DTFIX" << endl;
        this->dt = params->dtFix;
        cout << "SETTING POOL DTFIX" << endl;
    } else {
        this->dt = safetyFactor*this->getDt();

    }
    if (!stepTaken) {
        // pool->dtFix = this->dt;
        this->dtPrev = this->dt;
    }

    cout << "dt = " << dt << endl;

    // If the time step would exceed tEnd, truncate it
    if (this->t + this->dt > tEnd) {
        this->dt = tEnd - this->t;
    }

    this->applyInterfaceBCs();

    // Explicit step: advance the velocity-dependent terms using explicit time
    // discretization.
    this->updateF(pool);

    // Implicit step: update the pressure by solving Poisson equation.
    this->updateP();

    // Combine these together to update the fluid velocities.
    this->updateU();

    // Interpolate the velocities to the cell centers
    this->interpolateVelocities();

    if (this->nStructs > 0) {
        // Update the location of the interfaces
        // cout << "calling updatePool, methodOrd = " << methodOrd << endl;
        (this->pool)->updatePool(dt, iu, iv, p, methodOrd, true);

        // Apply object velocities to boundary points
        this->test_setInternal();
    }

    // Update the point in time
    this->t += this->dt;

    // Save the current time step
    this->dtPrev = this->dt;

    // The fist step has now been taken
    this->stepTaken = true;

    return this->t;
}

/**
 * Helper to get the upwind direction\
 * 
 * Cases:
 *  1 - positive upwind (x_{i+1})
 * -1 - negative upwind (x_{i-1})
 *  0 - no upwind
*/
int NSSolver::eno_get_upwind(bool sten[7], double ySten[7]) {
    int c = 3;

    int upDir;
    if (sten[c-1] && sten[c+1]) {
        if (ySten[c-1] > 0 && ySten[c+1] > 0) {
            upDir = -1;
        } else if (ySten[c-1] < 0 && ySten[c+1] < 0) {
            upDir = 1;
        } else {
            upDir = (ySten[c-1] > ySten[c+1]) ? -1 : 1;
        }
    } else if (sten[c-1] || sten[c+1]) {
        upDir = (sten[c-1]) ? -1 : 1;
    } else {
        upDir = 0;
    }

    return upDir;
}

/**
 * Eno scheme for the NS equations
 * 
 * d/dx u^2
*/
double NSSolver::eno_usqx(int i, int j) {
    // First build the stencil
    bool sten[7] = {false, false, false, false, false, false, false};
    pool->getPossibleStencil(i, j, 0, 3, sten);

    int c = 3;
    int xi = methodOrd + i;
    int yi = methodOrd + j;

    // Build the value and point set on the stencil
    double xSten[7] = {0};
    double ySten[7] = {0};

    if (!sten[c]) {
        return discs::firstOrder_conv_usqx(xi, yi, this->dx, this->u);
    }

    for (int l = 0; l < 7; l++) {
        if (sten[l]) {
            xSten[l] = x[i-c+l+1];
            ySten[l] = simutils::square(u[yi][xi-c+l]);
        }
    }

    // Choose the upwind direction properly (or return centered result)
    int upDir = eno_get_upwind(sten, ySten);

    if (upDir == 0) {
        // In case where ENO stencil fails, use centered approximation
        return discs::firstOrder_conv_usqx(xi, yi, this->dx, this->u);
    }

    return discs::thirdOrdENO(x[i+1], xSten, ySten, upDir, sten);
}

/**
 * Eno scheme for the NS equations
 * 
 * d/dy v^2
*/
double NSSolver::eno_vsqy(int i, int j) {
    // First build the stencil
    bool sten[7] = {false, false, false, false, false, false, false};
    pool->getPossibleStencil(i, j, 1, 3, sten);

    int c = 3;
    int xi = methodOrd + i;
    int yi = methodOrd + j;

    // Build the value and point set on the stencil
    double xSten[7] = {0};
    double ySten[7] = {0};
    for (int l = 0; l < 7; l++) {
        if (sten[l]) {
            xSten[l] = y[j-c+l+1];
            ySten[l] = simutils::square(v[yi-c+l][xi]);
        }
    }

    if (!sten[c]) {
        return discs::firstOrder_conv_vsqy(xi, yi, this->dy, this->v);
    }

    // Choose the upwind direction properly (or return centered result)
    int upDir = eno_get_upwind(sten, ySten);
    if (upDir == 0) {
        // In case where ENO stencil fails, use centered approximation
        return discs::firstOrder_conv_vsqy(xi, yi, this->dy, this->v);
    }

    return discs::thirdOrdENO(y[j+1], xSten, ySten, upDir, sten);
}

/**
 * Eno scheme for the NS equations
 * 
 * d/dx (uv)^2
*/
double NSSolver::eno_uvx(int i, int j) {
    int c = 3;
    int xi = methodOrd + i;
    int yi = methodOrd + j;

    // If center point is not available for interpolation, return the first order scheme
    bool check = pool->isUsableU(i, j) && pool->isUsableU(i, j+1) 
        && pool->isUsableU(i-1, j) && pool->isUsableU(i-1, j+1);
    
    if (!check) {
        return discs::firstOrder_conv_uvx(xi, yi, this->dx, this->u, this->v);
    }

    // First build the stencil
    int check_pnts[7] = {i-2, i-1, i, i, i, i+1, i+2};
    bool sten[7];
    for (int l = 0; l < 7; l++) {
        sten[l] = pool->isUsableV(check_pnts[l], j);
    }

    double xSten[7] = {0};
    double uVals[7] = {0};
    double vVals[7] = {0};

    // Compute center point
    if (!sten[c]) {
        return discs::firstOrder_conv_uvx(xi, yi, this->dx, this->u, this->v);
    }

    xSten[c] = simutils::midpoint(x[i], x[i+1]);
    uVals[c] = 0.5 * (0.5*(u[yi][xi-1] + u[yi+1][xi-1]) + 0.5*(u[yi][xi] + u[yi+1][xi]));
    vVals[c] = v[yi][xi];


    // Build the value and point set on the stencil
    double UC[2];

    for (int l = 0; l < c; l++) {
        if (sten[l]) {
            xSten[l] = x[i-c+l];
            uVals[l] = 0.5*(u[yi][xi-c+l] + u[yi+1][xi-c+l]);
            vVals[l] = 0.5*(v[yi][xi-c+l-1] + v[yi][xi-c+l]);
        }
    }

    for (int l = c+1; l < 7; l++) {
        if (sten[l]) {
            xSten[l+1] = x[i-c+l-1];
            uVals[l] = 0.5*(u[yi][xi-c+l-2] + u[yi+1][xi-c+l-2]);
            vVals[l] = 0.5*(v[yi][xi-c+l-2] + v[yi][xi-c+l-1]);
        }
    }

    // Choose the upwind direction properly (or return centered result)
    int upDir = eno_get_upwind(sten, uVals);
    if (upDir == 0) {
        // In case where ENO stencil fails, use centered approximation
        return discs::firstOrder_conv_uvx(xi, yi, this->dx, this->u, this->v);
    }

    // Using upwind direction, compute ENO
    double uv[7];
    for (int l = 0; l < 7; l++) {
        uv[l] = uVals[l] * vVals[l];
    }

    return discs::thirdOrdENO(simutils::midpoint(x[i], x[i+1]), xSten, uv, upDir, sten);
}

/**
 * Eno scheme for the NS equations
 * 
 * d/dy (uv)^2
*/
double NSSolver::eno_uvy(int i, int j) {
    int c = 3;
    int xi = methodOrd + i;
    int yi = methodOrd + j;

    // First build the stencil
    int check_pnts[7] = {j-2, j-1, j, j, j, j+1, j+2};

    // If center point is not available for interpolation, return the first order scheme
    bool check = pool->isUsableV(i, j) && pool->isUsableV(i+1, j) 
        && pool->isUsableV(i, j-1) && pool->isUsableV(i+1, j-1);
    
    if (!check) {
        return discs::firstOrder_conv_uvy(xi, yi, this->dy, this->u, this->v);
    }

    bool sten[7];
    for (int l = 0; l < 7; l++) {
        sten[l] = pool->isUsableU(i, check_pnts[l]);
    }

    double xSten[7] = {0};
    double uVals[7] = {0};
    double vVals[7] = {0};

    // Compute center point
    if (!sten[c]) {
        return discs::firstOrder_conv_uvy(xi, yi, this->dy, this->u, this->v);
    }
    xSten[c] = simutils::midpoint(y[j], y[j+1]);
    uVals[c] = u[yi][xi];
    vVals[c] = 0.5 * (0.5*(v[yi][xi] + v[yi][xi+1]) + 0.5*(v[yi-1][xi] + v[yi-1][xi+1]));


    // Build the value and point set on the stencil
    double UC[2];
    for (int l = 0; l < c; l++) {
        if (sten[l]) {
            xSten[l] = y[j-c+l]+1;
            uVals[l] = 0.5*(u[yi-c+l][xi] + u[yi-c+l][xi]);
            vVals[l] = 0.5*(v[yi-c+l][xi] + v[yi-c+l-1][xi+1]);
        }
    }

    for (int l = c+1; l < 7; l++) {
        if (sten[l]) {
            xSten[l+1] = y[j-c+l-1];
            uVals[l+1] = 0.5*(u[yi-c+l-2][xi] + u[yi-c+l-1][xi]);
            vVals[l+1] = 0.5*(v[yi-c+l-2][xi] + v[yi-c+l-2][xi+1]);
        }
    }

    // Choose the upwind direction properly (or return centered result)
    int upDir = eno_get_upwind(sten, vVals);
    if (upDir == 0) {
        // In case where ENO stencil fails, use centered approximation
        return discs::firstOrder_conv_uvy(xi, yi, this->dy, this->u, this->v);
    }

    // Using upwind direction, compute ENO
    double uv[7];
    for (int l = 0; l < 7; l++) {
        uv[l] = uVals[l] * vVals[l];
    }

    return discs::thirdOrdENO(simutils::midpoint(y[j], y[j+1]), xSten, uv, upDir, sten);
}


/**
 * Method which interpolate fluid velocities to the cell centers.
 * 
 * TODO: make sure that this implementation makes sense
*/
void NSSolver::interpolateCellCenterVelocities(int i, int j, double outU[2]) {
    double x, y;

    double uvals[2];
    double vvals[2];

    y = simutils::midpoint(this->y[j], this->y[j+1]);
    x = simutils::midpoint(this->x[i], this->x[i+1]);

    // The values being interpolated
    uvals[0] = this->u[this->methodOrd+j][i];
    uvals[1] = this->u[this->methodOrd+j][i+1];

    vvals[0] = this->v[j][this->methodOrd+i];
    vvals[1] = this->v[j+1][this->methodOrd+i];

    // Interpolate the value at the cell centers
    outU[0] = (uvals[0] + uvals[1])/2.0;
    outU[1] = (vvals[0] + vvals[1])/2.0;;
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
    double laplacian;
    double convective;
    for (j = 0; j < this->ny; j++) {
        for (i = 0; i < this->nx-1; i++) {
            int yi = j + mo;
            int xi = i + mo;
            if (pool->objAtIndex(i, j) == objects::FLUID_C && pool->isUpdateableU(i, j)) {
                laplacian = discs::firstOrder_lap_uxx(xi, yi, this->dx, this->u) 
                            + discs::firstOrder_lap_uyy(xi, yi, this->dy, this->u);

                if (params->useEno) {
                    convective = eno_usqx(i, j) + eno_uvy(i, j);
                } else {
                    convective = discs::firstOrder_conv_usqx(xi, yi, this->dx, this->u)
                        + discs::firstOrder_conv_uvy(xi, yi, this->dy, this->u, this->v);
                }
                this->FU[yi][xi] = this->u[yi][xi] + this->dt*( (1.0/Re)*laplacian - convective );
            }
        }
    }

    for (j = 0; j < this->ny-1; j++) {
        for (i = 0; i < this->nx; i++) {
            int yi = j + mo;
            int xi = i + mo;
            if (pool->objAtIndex(i, j) == objects::FLUID_C && pool->isUpdateableV(i, j)) {
                laplacian = discs::firstOrder_lap_vxx(xi, yi, this->dx, this->v)
                            + discs::firstOrder_lap_vyy(xi, yi, this->dy, this->v);

                if (params->useEno) {
                    convective = eno_uvx(i, j) + eno_vsqy(i, j);
                } else {
                    convective = discs::firstOrder_conv_uvx(xi, yi, this->dx, this->u, this->v)
                        + discs::firstOrder_conv_vsqy(xi, yi, this->dy, this->v);
                }
                this->FV[yi][xi] = this->v[yi][xi] + this->dt*( (1.0/Re)*laplacian - convective );
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
    double x, y;
    int mo = methodOrd;

    // Open the output file
    ofstream outFile;
    outFile.open(fname);

    for (int j = 0; j < this->ny; j++) {
        for (int i = 0; i < this->nx; i++) {
            int xi = i + mo;
            int yi = j + mo;
            // if (pool->objAtIndex(i, j) == objects::FLUID_C) {
            x = simutils::midpoint(this->x[i], this->x[i+1]);
            y = simutils::midpoint(this->y[j], this->y[j+1]);
            outFile << x << ", ";
            outFile << y << ", ";
            outFile << this->iu[j][i] << ", ";
            outFile << this->iv[j][i] << ", ";
            outFile << this->p[yi][xi] << "\n";
            // }
        }
    }

    outFile.close();
}

void NSSolver::writePoolToFile(const char *poolName, const char *poolVelName) {
    pool->outputPool(poolName);
    pool->outputPoolVelocity(poolVelName);
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

NSSolver::~NSSolver() {
    delete this->pool;
    delete[] this->x;
    delete[] this->y;
    delete[] this->baryWeights;

    delete pressureSolver;
    delete pressureParams;

    // Free multidimensional arrays using helper
    int hg = this->ny + 2*this->methodOrd;
    int wg = this->nx + 2*this->methodOrd;

    // TODO: how necissary are the derivative arrays? My guess is not very? Will keep
    //       them around for ease of implementation
    simutils::free_double(hg, wg-1, this->u);
    simutils::free_double(this->ny, this->nx, this->iu);
    simutils::free_double(hg, wg-1, this->FU);
    simutils::free_double(hg-1, wg, this->v);
    simutils::free_double(this->ny, this->nx, this->iv);
    simutils::free_double(hg-1, wg, this->FV);
    simutils::free_double(hg, wg, this->p);
}
