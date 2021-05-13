#include "PressureSolver.h"
#include "../Utils/SimUtilities.h"
#include <math.h>
#include "../LASolver/SparseItObj.h"
#include "Boundary.h"
#include <cassert>

// TODO: consider the use of Dirichlet at the boundary strip

/** 
 * Method which creates a new matrix object.
 * 
 * 
*/
void PressureSolver::setUpMatrixStruct(int nx, int ny) {
    int n = nx*ny;
    int i, ii, j, off, offp1;

    // Set up the system matrix. Matrix is set up to have a 5-stencil at every point
    // except the sides & corners to permit changes in the stencils
    MatrixStruc *matrixBuilder = NULL;
    try {
        matrixBuilder = new MatrixStruc(n /*Number of unknowns*/, 0 /*Add non-zeros to diagonals*/);

        /* Create the data structure */
        {
            // SW
            matrixBuilder->set_entry(0, 1); // <i+1, j>
            matrixBuilder->set_entry(0, nx); // <i, j+1>

            // S
            for (int i = 1; i < nx-1; i++) {
                matrixBuilder->set_entry(i, i-1); // <i-1, j>
                matrixBuilder->set_entry(i, i+1); // <i+1, j>
                
                matrixBuilder->set_entry(i, i+nx); // <i, j+1>
            }

            // SE
            matrixBuilder->set_entry(nx-1, nx-2); // <i-1, j>
            matrixBuilder->set_entry(nx-1, (nx-1)+nx); // <i, j+1>

            for (j = 1; j < ny-1; j++) {
                off = nx*j;

                // E
                matrixBuilder->set_entry(off, off+1); // <i+1, j>
                matrixBuilder->set_entry(off, off-nx); // <i, j-1>
                matrixBuilder->set_entry(off, off+nx); // <i, j+1>

                // C
                for (i = 1; i < nx-1; i++) {
                    ii = off+i; // The current row number.
                    matrixBuilder->set_entry(ii, ii-1); // <i-1, j>
                    matrixBuilder->set_entry(ii, ii+1); // <i+1, j>

                    matrixBuilder->set_entry(ii, ii-nx); // <i, j-1>
                    matrixBuilder->set_entry(ii, ii+nx); // <i, j+1>
                }

                // W
                offp1 = nx*(j+1)-1;

                matrixBuilder->set_entry(offp1, offp1-1); // <i-1, j>
                matrixBuilder->set_entry(offp1, offp1-nx); // <i, j-1>
                matrixBuilder->set_entry(offp1, offp1+nx); // <i, j+1>
            }

            // NW
            off = (ny-1)*nx;
            matrixBuilder->set_entry(off, off+1); // <i+1, j>
            matrixBuilder->set_entry(off, off-nx); // <i, j-1>

            // N
            for (i = off+1; i < ny*nx-1; i++) {
                matrixBuilder->set_entry(i, i-nx); // <i, j-1>
                matrixBuilder->set_entry(i, i+1); // <i+1, j>
                matrixBuilder->set_entry(i, i-1); // <i-1, j>
            }

            // W
            off = ny*nx - 1;
            matrixBuilder->set_entry(off, off-1); // <i-1, j>
            matrixBuilder->set_entry(off, off-nx); // <i, j-1>
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
    this->matrix = NULL;

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
}

/**
 * Base constructor. Initialize the PCG pressure solver.
*/
PressureSolver::PressureSolver(int nx, double *x, int ny, double *y) {
    // Counters and increments
    int i;

    this->nx = nx;
    this->ny = ny;

    this->n = nx*ny;

    this->x = new double[nx+1];
    this->y = new double[ny+1];

    for (i = 0 ; i < nx+1; i++)
        this->x[i] = x[i];

    for (i = 0 ; i < ny+1; i++)
        this->y[i] = y[i];
    
    this->pFlat = new double[this->n];
    this->tol = new double[this->n];

    // Set the update tolerances to 0, forces use of risidual reduction criteria
    for (i = 0; i < n; i++) {
        this->tol[i] = 0.0;
    }

    this->hx = this->x[1] - this->x[0];
    this->hy = this->y[1] - this->y[0];

    // Set up the memory for Justin's solver
    this->setUpMatrixStruct(nx, ny);
}

/**
 * Set up the pressure matrix.
 * 
 * TODO: include boundary information. Using the helper functions with appropriately set boundary conditions
 * should make it very easy.
*/
void PressureSolver::setUpPressureMatrix(Pool2D *pool) {
    int i, j, off, offp;
    int row = 0;

    // SW
    if (pool->objAtIndex(0,0) == objects::FLUID_C)
        this->discSW(row, this->matrix);
    else
        this->discSolid(row, this->matrix);
    row++;

    // S
    for (i = 1; i < nx-1; i++) {
        if (pool->objAtIndex(i, 0) == objects::FLUID_C)
            this->discS(row, this->matrix);
        else
            this->discSolid(row, this->matrix);
        row++;
    }

    // SE
    if (pool->objAtIndex(nx-1,0) == objects::FLUID_C)
        this->discSE(row, this->matrix);
    else
        this->discSolid(row, this->matrix);
    row++;

    for (j = 1; j < ny-1; j++) {
        off = nx*j;
        
        // W
        if (pool->objAtIndex(0,j) == objects::FLUID_C)
            this->discW(row, matrix);
        else
            this->discSolid(row, this->matrix);
        row++;

        // C
        for (i = 1; i < nx-1; i++) {
            this->poolDisc(row, i, j, pool, matrix);
            row++;
        }

        // E
        offp = nx*(j+1)-1;
        this->discE(row, matrix);

        if (pool->objAtIndex(nx-1,j) == objects::FLUID_C)
            this->discE(row, matrix);
        else
            this->discSolid(row, this->matrix);
        row++;
    }

    // NW
    off = nx*(ny-1);
    if (pool->objAtIndex(0,ny-1) == objects::FLUID_C)
        this->discNW(row, matrix);
    else
        this->discSolid(row, this->matrix);
    row++;

    // N
    for (i = 1; i < nx-1; i++) {
        if (pool->objAtIndex(i,ny-1) == objects::FLUID_C)
            this->discN(row, matrix);
        else
            this->discSolid(row, this->matrix);
        row++;
    }

    // NE
    if (pool->objAtIndex(nx-1,ny-1) == objects::FLUID_C)
        this->discNE(row, matrix);
    else
        this->discSolid(row, this->matrix);
    row++;
}

/**
 * Set up the pressure matrix.
 * 
 * We only need to compute this for the fluid cells.
 * 
 * TODO: generalize indexing into FU, FV based on method order. Going to have more or less ghost cells.
 * TODO: generalize the computation of FU to second order using sided stencils when possible
*/
void PressureSolver::setUpPressureRHS(double dt, double **FU, double **FV, double **p, Pool2D *pool) {
    int i, j;

    int mo = 1;

    double accum = 0.0;
    int tot = 0;
    int row = 0;

    int nDir[2];
    objects::FSIObject obj;

    double velVec[2];
    double nDirNorm[2];
    double rhs;
    double h;


    try {
        row = 0;
        for (j = 0; j < this->ny; j++) {
            for (i = 0; i < this->nx; i++) {
                obj = pool->objAtIndex(i, j);

                if (obj == objects::FLUID_C) {
                    (this->matrix)->bValue(row)
                        = ((FU[j+mo][i+mo] - FU[j+mo][i])/this->hx + (FV[j+mo][i+mo]-FV[j][i+mo])/this->hy)/dt;

                    accum += ((FU[j+mo][i+mo] - FU[j+mo][i])/this->hx + (FV[j+mo][i+mo]-FV[j][i+mo])/this->hy);
                    tot += 1;
                } else if (pool->isInterface(obj)) {
                    // Get the normal direction
                    pool->getNormalDir(obj, nDir);
                    nDirNorm[0] = nDir[0];
                    nDirNorm[1] = nDir[1];
                    simutils::normalize2D(nDirNorm);
                    velVec[0] = FU[j+mo+nDir[1]][i+mo+nDir[0]] - pool->getObjU(i, j);
                    velVec[1] = FV[j+mo+nDir[1]][i+mo+nDir[0]] - pool->getObjV(i, j);
                    rhs = simutils::ddot2d(velVec, nDirNorm);
                    h = sqrt(simutils::square(hx*nDir[0]) + simutils::square(hy*nDir[1]));
                    // (this->matrix)->bValue(row) = -p[j+mo+nDir[1]][i+mo+nDir[0]];
                    (this->matrix)->bValue(row) = p[j+mo+nDir[1]][i+mo+nDir[0]] - (h/dt)*rhs;
                } else {
                    (this->matrix)->bValue(row) = 0;
                }
                row++;
            }
        }
    } catch(bad_alloc) {
        exit(1);
    } catch(General_Exception excep) {
        cout << excep.p << endl;
        exit(1);
    }
}

/**
 * Solve the pressure equation using PCG.
*/
void PressureSolver::solvePCG(double dt, Pool2D *pool, double **FU, double **FV,
                                double **p, bool matrixRecompute, ParamIter &params) {
    // Compute the matrix. Right now the iluRecompute flag
    // is being used to indicate setting up the pressure matrix at all
    if (matrixRecompute) {
        this->setUpPressureMatrix(pool);

        // Perform the ILU factorization
        if (params.drop_ilu == 0) {
            matrix->sfac(params);
        }
    }

    // Set up the rhs of the linear system
    this->setUpPressureRHS(dt, FU, FV, p, pool);

    // Copy the pressure to the flat matrix
    // TODO: should the pressure contain the boundary squares? Matters more in 3D.
    int off;
    for (int j = 0; j < ny; j++) {
        off = nx*j;
        for (int i = 0; i < nx; i++) {
            this->pFlat[off+i] = p[j+1][i+1];
        }
    }

    try {
        for (int i = 0; i < n; i++) {
            this->tol[i] = 0; // set update tolerance to 0. Forces the use of risidual reduction criteria.
                              // TODO: find out what this means.
                              // TODO: do this in the constructor.
                              // TODO: set the structure tolerances to something large.
        }
        matrix->set_toler(this->tol); // Set the update tolerance

        int niter; // The number of iterations done in the solve

        // Solve the linear system
        matrix->solve(params, this->pFlat, niter);

        cout << "Solved pressure in " << niter << endl;

        assert(niter >= 0);
    }
    catch(bad_alloc) {
        exit(1);
    }
    catch(General_Exception excep) {
        exit(1);
    }

    // Copy the solution to the pressure matrix
    for (int j = 0; j < ny; j++) {
        off = nx*j;
        for (int i = 0; i < nx; i++) {
            p[j+1][i+1] = this->pFlat[off+i];
        }
    }

    // Take the mean
    cout << "Mean pressure (sorta) " << simutils::mean(ny+2, nx+2, p);
}

PressureSolver::~PressureSolver() {
    delete[] this->x;
    delete[] this->y;
    delete[] this->pFlat;
    delete[] this->tol;
    delete matrix;
}

/* Many discretization routines. */

/**
 * Choose discretization based on the kind of fluid structure this point is.
 * 
 * The names for the discretization schemes chosen to be consistant with what is used for the
 * default fluids case.
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
void PressureSolver::poolDisc(int row, int xInd, int yInd, Pool2D *pool, MatrixIter *matrix) {
    objects::FSIObject type = pool->objAtIndex(xInd, yInd);

    // int intType = objects::NORTH*objects::SOUTH*objects::EAST*objects::WEST;

    if (type == objects::STRUCTURE) {
        discSolid(row, matrix);
    } else if (pool->isInterface(type)) {
        discSolid(row, matrix);
    } else {
        discC(row, matrix);
    }

    // If the current cell does not contain a fluid, trivial discretization.
    // Else we choose based on the relative location of the interface points for this
    // fluid cell.
    // if (type != objects::FLUID_C) {
    //     discSolid(row, matrix);
    // } else {
    //     // Check in each cardinal direction if there is an interface. If so, modify intType. If
    //     // there is no interface, it will keep its value 1 label and we relabel it to FLUID_C. Otherwise
    //     // use the appropriate directional discretization.
    //     intType = (pool->objAtIndex(xInd, yInd+1) != objects::FLUID_C) ? (int) intType/objects::NORTH : intType;
    //     intType = (pool->objAtIndex(xInd, yInd-1) != objects::FLUID_C) ? (int) intType/objects::SOUTH : intType;
    //     intType = (pool->objAtIndex(xInd+1, yInd) != objects::FLUID_C) ? (int) intType/objects::EAST : intType;
    //     intType = (pool->objAtIndex(xInd-1, yInd) != objects::FLUID_C) ? (int) intType/objects::WEST : intType;
        
    //     // No interface - set to FLUID_C
    //     if (intType == objects::NORTH*objects::SOUTH*objects::EAST*objects::WEST) {
    //         intType = -1;
    //     }

    //     switch ((objects::FSIObject) intType) {
    //         case objects::FLUID_SW:  discSW(row, matrix);
    //                                  break;

    //         case objects::FLUID_S:   discS(row, matrix);
    //                                  break;

    //         case objects::FLUID_SE:  discSE(row, matrix);
    //                                  break;

    //         case objects::FLUID_W:   discW(row, matrix);
    //                                  break;

    //         case objects::FLUID_C:   discC(row, matrix);
    //                                  break;

    //         case objects::FLUID_E:   discE(row, matrix);
    //                                  break;

    //         case objects::FLUID_NW:  discNW(row, matrix);
    //                                  break;

    //         case objects::FLUID_N:   discN(row, matrix);
    //                                  break;

    //         case objects::FLUID_NE:  discNE(row, matrix);
    //                                  break;

    //         default:                 assert(false);
    //                                  break;
    //     }
    // }
}

void PressureSolver::discSW(int row, MatrixIter *matrix) {
    int colIndex;
    double hxsq = simutils::square(this->hx);
    double hysq = simutils::square(this->hy);

    for (int i = matrix->rowBegin(row); i < matrix->rowEndPlusOne(row); i++) {
        colIndex = matrix->getColIndex(i);

        if (colIndex == row) {
            matrix->aValue(i) = -(1/hxsq + 1/hysq);
        } else if (colIndex == row+1) {
            matrix->aValue(i) = 1/hxsq;
        } else if (colIndex == row + this->nx) {
            matrix->aValue(i) = 1/hysq;
        } else {
            matrix->aValue(i) = 0;
        }
    }
}

void PressureSolver::discS(int row, MatrixIter *matrix) {
    int colIndex;
    double hxsq = simutils::square(this->hx);
    double hysq = simutils::square(this->hy);

    for (int i = matrix->rowBegin(row); i < matrix->rowEndPlusOne(row); i++) {
        colIndex = matrix->getColIndex(i);

        if (colIndex == row) {
            matrix->aValue(i) = -(2/hxsq + 1/hysq);
        } else if (colIndex == row + 1 || colIndex == row - 1) {
            matrix->aValue(i) = 1/hxsq;
        } else if (colIndex == row + this->nx) {
            matrix->aValue(i) = 1/hysq;
        } else {
            matrix->aValue(i) = 0;
        }
    }
}

void PressureSolver::discSE(int row, MatrixIter *matrix) {
    int colIndex;
    double hxsq = simutils::square(this->hx);
    double hysq = simutils::square(this->hy);

    for (int i = matrix->rowBegin(row); i < matrix->rowEndPlusOne(row); i++) {
        colIndex = matrix->getColIndex(i);

        if (colIndex == row) {
            matrix->aValue(i) = -(1/hxsq + 1/hysq);
        } else if (colIndex == row - 1) {
            matrix->aValue(i) = 1/hxsq;
        } else if (colIndex == row + this->nx) {
            matrix->aValue(i) = 1/hysq;
        } else {
            matrix->aValue(i) = 0;
        }
    }

}

void PressureSolver::discW(int row, MatrixIter *matrix) {
    int colIndex;
    double hxsq = simutils::square(this->hx);
    double hysq = simutils::square(this->hy);
    for (int i = matrix->rowBegin(row); i < matrix->rowEndPlusOne(row); i++) {
        colIndex = matrix->getColIndex(i);

        if (colIndex == row) {
            matrix->aValue(i) = -(1/hxsq + 2/hysq);
        } else if (colIndex == row + 1) {
            matrix->aValue(i) = 1/hxsq;
        } else if (colIndex == row - this->nx || colIndex == row + this->nx) {
            matrix->aValue(i) = 1/hysq;
        } else {
            matrix->aValue(i) = 0;
        }
    }
}

void PressureSolver::discC(int row, MatrixIter *matrix) {
    int colIndex;
    double hxsq = simutils::square(this->hx);
    double hysq = simutils::square(this->hy);
    for (int i = matrix->rowBegin(row); i < matrix->rowEndPlusOne(row); i++) {
        colIndex = matrix->getColIndex(i);

        if (colIndex == row) {
            matrix->aValue(i) = -2*(1/hxsq + 1/hysq);
        } else if (colIndex == row-1 || colIndex == row+1) {
            matrix->aValue(i) = 1/hxsq;
        } else if (colIndex == row-this->nx || colIndex == row+this->nx) {
            matrix->aValue(i) = 1/hysq;
        } else {
            matrix->aValue(i) = 0;
        }
    }
}

void PressureSolver::discE(int row, MatrixIter *matrix) {
    int colIndex;
    double hxsq = simutils::square(this->hx);
    double hysq = simutils::square(this->hy);
    for (int i = matrix->rowBegin(row); i < matrix->rowEndPlusOne(row); i++) {
        colIndex = matrix->getColIndex(i);

        if (colIndex == row) {
            matrix->aValue(i) = -(1/hxsq + 2/hysq);
        } else if (colIndex == row - 1) {
            matrix->aValue(i) = 1/hxsq;
        } else if (colIndex == row - this->nx || colIndex == row + this->nx) {
            matrix->aValue(i) = 1/hysq;
        } else {
            matrix->aValue(i) = 0;
        }
    }
}

void PressureSolver::discNW(int row, MatrixIter *matrix) {
    int colIndex;
    double hxsq = simutils::square(this->hx);
    double hysq = simutils::square(this->hy);
    for (int i = matrix->rowBegin(row); i < matrix->rowEndPlusOne(row); i++) {
        colIndex = matrix->getColIndex(i);

        if (colIndex == row) {
            matrix->aValue(i) = -(1/hxsq + 1/hysq);
        } else if (colIndex == row + 1) {
            matrix->aValue(i) = 1/hxsq;
        } else if (colIndex == row - this->nx) {
            matrix->aValue(i) = 1/hysq;
        } else {
            matrix->aValue(i) = 0;
        }
    }
}

void PressureSolver::discN(int row, MatrixIter *matrix) {
    int colIndex;
    double hxsq = simutils::square(this->hx);
    double hysq = simutils::square(this->hy);
    for (int i = matrix->rowBegin(row); i < matrix->rowEndPlusOne(row); i++) {
        colIndex = matrix->getColIndex(i);

        if (colIndex == row) {
            matrix->aValue(i) = -(2/hxsq + 1/hysq);
        } else if (colIndex == row - 1 || colIndex == row+1) {
            matrix->aValue(i) = 1/hxsq;
        } else if (colIndex == row - nx) {
            matrix->aValue(i) = 1/hysq;
        } else {
            matrix->aValue(i) = 0;
        }
    }
}

void PressureSolver::discNE(int row, MatrixIter *matrix) {
    int colIndex;
    double hxsq = simutils::square(this->hx);
    double hysq = simutils::square(this->hy);
    for (int i = matrix->rowBegin(row); i < matrix->rowEndPlusOne(row); i++) {
        colIndex = matrix->getColIndex(i);

        if (colIndex == row) {
            matrix->aValue(i) = -(1/hxsq + 1/hysq);
        } else if (colIndex == row - 1) {
            matrix->aValue(i) = 1/hxsq;
        } else if (colIndex == row - nx) {
            matrix->aValue(i) = 1/hysq;
        } else {
            matrix->aValue(i) = 0;
        }
    }
}

void PressureSolver::discSolid(int row, MatrixIter *matrix) {
    int colIndex;
    for (int i = matrix->rowBegin(row); i < matrix->rowEndPlusOne(row); i++) {
        colIndex = matrix->getColIndex(i);

        // Just avoiding a singular matrix, everything else is set to 0.
        if (colIndex == row) {
            matrix->aValue(i) = 1;
        } else {
            matrix->aValue(i) = 0;
        }
    }
}