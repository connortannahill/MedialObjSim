
#include "PressureSolver.h"
#include "SimUtilities.h"
#include <math.h>
#include "SparseItObj.h"
#include "Boundary.h"

/**
 * Base constructor. Initialize the PCG pressure solver.
*/
PressureSolver::PressureSolver(int nx, double *x, int ny, double *y) {
    // Counters and increments
    int i, ii, off;
    // std::cout << "Enterting ctor" << std::endl;

    this->nx = nx;
    this->ny = ny;

    this->n = (nx+2)*(ny+2);

    // Create the recon matrix
    this->reconMatrix = simutils::new_constant(n, n, 0);

    // Create mesh
    this->x = new double[nx+1];
    this->y = new double[ny+1];

    for (i = 0 ; i < nx+1; i++)
        this->x[i] = x[i];

    for (i = 0 ; i < ny+1; i++)
        this->y[i] = y[i];
    
    this->pFlat = simutils::new_constant(this->n, 0);

    // Set the update tolerances to 0, forces use of risidual reduction criteria
    this->tol = simutils::new_constant(this->n, 0);

    this->hx = this->x[1] - this->x[0];
    this->hy = this->y[1] - this->y[0];

    // Set up the system matrix. Matrix is set up to have a 5-stencil at every point
    // except the sides & corners to permit future generality.

    // Using Ghost cells which are pre-labelled

    MatrixStruc *matrixBuilder = NULL;
    try {
        matrixBuilder = new MatrixStruc(n /*Number of unknowns*/, 0 /*Add non-zeros to diagonals*/);

        /* Create the data structure */
        {
            // BOT BC
            /* Added by default */

            for (int j = 1; j <= ny; j++) {
                off = (nx+2)*j;

                // Internal loop excludes the Dirichlet boundary conditions, which will be added automatically.
                for (int i = 1; i <= nx; i++) {
                    ii = off+i; // The current row number.
                    matrixBuilder->set_entry(ii, ii-1);
                    matrixBuilder->set_entry(ii, ii+1);

                    matrixBuilder->set_entry(ii, ii-(nx+2));
                    matrixBuilder->set_entry(ii, ii+(nx+2));
                }
            }
    // std::cout << "Hi";

            // TOP BC
            /* Covered by the constructor */
        }

        // Pack the data constructor: convert the LinkedList array to csr format
        matrixBuilder->pack();

    } catch(std::bad_alloc) {
        delete matrixBuilder; matrixBuilder = NULL;
        exit(1);
    } catch(General_Exception excep) {
        std::cout << "General exception" << std::endl;
        std::cout << excep.p << std::endl;
        delete matrixBuilder; matrixBuilder = NULL;
        exit(1);
    }

    // Create the actual matrix
    this->matrix = NULL;

    try {
        // std::cout << "Making the matrix" << std::endl;
        matrix = new MatrixIter(*matrixBuilder);
        // std::cout << "Made the matrix" << std::endl;
        assert(matrix != NULL);
        delete matrixBuilder; 
        matrixBuilder = NULL;
    } catch(std::bad_alloc) {
        // std::cout << "Hit excep 1" << std::endl;
        delete matrixBuilder; matrixBuilder = NULL;
        exit(1);
    } catch(General_Exception excep) {
        // std::cout << "Hit excep 2" << std::endl;
        delete matrixBuilder; matrixBuilder= NULL;
        exit(1);
    } 
    // std::cout << "Exiting ctor" << std::endl;
}

/**
 * Set up the pressure matrix.
 * 
 * TODO: include boundary information. Using the helper functions with appropriately set boundary conditions
 * should make it very easy.
*/
void PressureSolver::setUpPressureMatrix(Pool2D &pool) {
    int i, j, ii, off, offp;

    // std::cout << this->matrix << std::endl;
    // Set the matrix to 0 initially

    // BOTTOM BC
    for (i = 0; i < nx+2; i++) {
        matrix->aValue(matrix->rowBegin(i)) = 1.0;
    }
    // // SW
    // this->discSW(0, 0, 0, this->matrix);

    // // S
    // for (i = 1; i < nx-1; i++) {
    //     this->discS(i, 0, i, this->matrix);
    // }

    // // SE
    // this->discSE(nx-1, 0, nx-1, this->matrix);
    // matrix->printMat(n);

    // Internal
    for (j = 1; j <= ny; j++) {
        off = (nx+2)*j;

        // LEFT
        matrix->aValue(matrix->rowBegin(off)) = 1.0;

        // Internal
        for (int i = 1; i <= nx; i++) {
            ii = off+i;
            this->discC(i, j, ii, matrix);
        }

        offp = (nx+2)*(j+1)-1;

        // RIGHT
        matrix->aValue(matrix->rowBegin(offp)) = 1.0;
    }
    // matrix->printMat(n);

    // TOP
    off = (nx+2)*(ny+1);
    for (i = 0; i < nx+2; i++) {
        matrix->aValue(matrix->rowBegin(off+i)) = 1.0;

    }
    // matrix->printMat(n);


    // // NW
    // off = (nx+2)*(ny+1);
    // this->discNW(0, ny-1, off, matrix);

    // // N
    // for (i = 1; i < nx-1; i++) {
    //     this->discN(i, ny-1, off+i, matrix);
    // }

    // // NE
    // this->discNE(nx-1, ny-1, nx*ny-1, matrix);
}

/**
 * Set up the pressure matrix.
 * 
 * Note: at the edges of the domain we are including the boundary information
*/
void PressureSolver::setUpPressureRHS(double dt, double **FU, double **FV, double **p, Pool2D &pool) {
    int i, j, off;

    try {
        // std::cout << "Setting up the RHS" << std::endl;
        // BOT:
        for (i = 0; i < nx+2; i++) {
            // pFlat[i] = p[1][i];
            (this->matrix)->bValue(i) = p[1][i];
        }

        // Internal
        for (j = 1; j <= this->ny; j++) {
            off = (nx+2)*j;

            // LHS
            // this->pFlat[off] = p[j][1];
            (this->matrix)->bValue(off) = p[j][1];

            // Internal
            for (i = 1; i <= this->nx; i++) {
                // (this->matrix)->bValue(off+i) = ((FU[j+1][i+1] - FU[j+1][i])/this->hx + (FV[j+1][i+1]-FV[j][i+1])/this->hy)/dt;
                (this->matrix)->bValue(off+i) = ((FU[j][i] - FU[j][i-1])/this->hx + (FV[j][i]-FV[j-1][i])/this->hy)/dt;
            }

            // RHS
            // this->pFlat[nx*(j+1)-1] = p[j][nx];
            (this->matrix)->bValue(nx*(j+1)-1) = p[j][nx];
        }

        // TOP:
        off = (nx+2)*(ny+1);
        for (i = 0; i < nx+2; i++) {
            // pFlat[off+i] = p[ny][i];
            (this->matrix)->bValue(off+i) = p[ny][i];
        }
    } catch(std::bad_alloc) {
        exit(1);
    } catch(General_Exception excep) {
        std::cout << excep.p << std::endl;
        exit(1);
    }
}

/**
 * Solve the pressure equation using PCG.
 * 
 * // TODO: make the solver aware of the boundaries
 * // TODO: add the ability for ILU recomputation
*/
void PressureSolver::solvePCG(double dt, Pool2D &pool, double **FU, double **FV,
                                double **p, bool iluRecompute, ParamIter &params) {
    // Compute the matrix. Right now the iluRecompute flag
    // is being used to indicate setting up the pressure matrix at all
    if (iluRecompute) {
        std::cout << "Setting up the pressure matrix" << std::endl;
        this->setUpPressureMatrix(pool);
    }


    // assert(false);

    // Set up the rhs of the linear system
    this->setUpPressureRHS(dt, FU, FV, p, pool);

    matrix->printMat(n);
    matrix->printRHS(n);

    // Copy the pressure to the flat matrix
    // TODO: should the pressure contain the boundary squares? Matters more in 3D.
    int off;
    std::cout << "Copying RHS from pressure" << std::endl;
    for (int j = 0; j < ny+2; j++) {
        off = (nx+2)*j;
        for (int i = 0; i < nx+2; i++) {
            this->pFlat[off+i] = p[j][i];
            std::cout << this->pFlat[off+i] << std::endl;
        }
        // std::cout << std::endl;
    }

    // Perform the ILU factorization
    if (params.drop_ilu == 0) {
        matrix->sfac(params);
    }

    try {
        for (int i = 0; i < n; i++) {
            this->tol[i] = 0; // set update tolerance to 0. Forces the use of risidual reduction criteria.
                              // TODO: find out what this means.
        }
        matrix->set_toler(this->tol); // Set the update tolerance

        int niter; // The number of iterations done in the solve

        // Solve the linear system
        matrix->solve(params, this->pFlat, niter);

        std::cout << "Solved pressure in " << niter << std::endl;

        assert(niter >= 0);
    }
    catch(std::bad_alloc) {
        exit(1);
    }
    catch(General_Exception excep) {
        exit(1);
    }

    // Copy the solution to the pressure matrix
    for (int j = 0; j < ny+2; j++) {
        off = (nx+2)*j;
        for (int i = 0; i < nx+2; i++) {
            p[j][i] = this->pFlat[off+i];
        }
    }

    // simutils::printMat(this->ny+2, this->nx+2, p);
    // std::cout << "Printing the recon matrix" << std::endl;
    // simutils::printMat(this->n, this->n, this->reconMatrix);
    std::cout << "Mean pressure " << simutils::mean(ny+2, nx+2, p);
    assert(false);

    // Take the mean
}

PressureSolver::~PressureSolver() {
    delete[] this->x;
    delete[] this->y;
    delete[] this->pFlat;
    delete[] this->tol;
    simutils::free_double(this->n, this->n, this->reconMatrix);
    delete matrix;
}

/* Many discretization routines. TODO: generalize this all using ENUM cases */
void PressureSolver::discSW(int xInd, int yInd, int row, MatrixIter *matrix) {
    int colIndex;
    double hxsq = simutils::square(this->hx);
    double hysq = simutils::square(this->hy);

    for (int i = matrix->rowBegin(row); i < matrix->rowEndPlusOne(row); i++) {
        colIndex = matrix->getColIndex(i);

        if (colIndex == row) {
            matrix->aValue(i) = -(1/hxsq + 1/hysq);
            this->reconMatrix[row][colIndex] = -2;

        } else if (colIndex == row+1) {
            matrix->aValue(i) = 1/hxsq;
            this->reconMatrix[row][colIndex] = 1; 

        } else if (colIndex == row + this->nx) {
            matrix->aValue(i) = 1/hysq;
            this->reconMatrix[row][colIndex] = 1; 
        } else {
            assert(false);
        }
    }
}

void PressureSolver::discS(int xInd, int yInd, int row, MatrixIter *matrix) {
    int colIndex;
    double hxsq = simutils::square(this->hx);
    double hysq = simutils::square(this->hy);

    for (int i = matrix->rowBegin(row); i < matrix->rowEndPlusOne(row); i++) {
        colIndex = matrix->getColIndex(i);

        if (colIndex == row) {
            matrix->aValue(i) = -(2/hxsq + 1/hysq);
            this->reconMatrix[row][colIndex] = -3;
        } else if (colIndex == row + 1 || colIndex == row - 1) {
            matrix->aValue(i) = 1/hxsq;
            this->reconMatrix[row][colIndex] = 1;
        } else if (colIndex == row + this->nx || colIndex == row - this->nx) {
            matrix->aValue(i) = 1/hysq;
            this->reconMatrix[row][colIndex] = 1;
        } else {
            assert(false);
        }
    }

}

void PressureSolver::discSE(int xInd, int yInd, int row, MatrixIter *matrix) {
    int colIndex;
    double hxsq = simutils::square(this->hx);
    double hysq = simutils::square(this->hy);

    for (int i = matrix->rowBegin(row); i < matrix->rowEndPlusOne(row); i++) {
        colIndex = matrix->getColIndex(i);

        if (colIndex == row) {
            matrix->aValue(i) = -(1/hxsq + 1/hysq);
            this->reconMatrix[row][colIndex] = -2;
        } else if (colIndex == row - 1) {
            matrix->aValue(i) = 1/hxsq;
            this->reconMatrix[row][colIndex] = 1;
        } else if (colIndex == row + this->nx) {
            matrix->aValue(i) = 1/hysq;
            this->reconMatrix[row][colIndex] = 1;
        } else {
            assert(false);
        }
    }

}

void PressureSolver::discW(int xInd, int yInd, int row, MatrixIter *matrix) {
    int colIndex;
    double hxsq = simutils::square(this->hx);
    double hysq = simutils::square(this->hy);
    for (int i = matrix->rowBegin(row); i < matrix->rowEndPlusOne(row); i++) {
        colIndex = matrix->getColIndex(i);

        if (colIndex == row) {
            matrix->aValue(i) = -(1/hxsq + 2/hysq);
            this->reconMatrix[row][colIndex] =-3;
        } else if (colIndex == row + 1) {
            matrix->aValue(i) = 1/hxsq;
            this->reconMatrix[row][colIndex] =1;
        } else if (colIndex == row - this->nx || colIndex == row + this->nx) {
            matrix->aValue(i) = 1/hysq;
            this->reconMatrix[row][colIndex] =1;
        } else {
            assert(false);
        }
    }
}

void PressureSolver::discC(int xInd, int yInd, int row, MatrixIter *matrix) {
    int colIndex;
    double hxsq = simutils::square(this->hx);
    double hysq = simutils::square(this->hy);
    for (int i = matrix->rowBegin(row); i < matrix->rowEndPlusOne(row); i++) {
        colIndex = matrix->getColIndex(i);

        if (colIndex == row) {
            matrix->aValue(i) = -2*(1/hxsq + 1/hysq);
            this->reconMatrix[row][colIndex] = -4;
        } else if (colIndex == row-1 || colIndex == row+1) {
            matrix->aValue(i) = 1/hxsq;
            this->reconMatrix[row][colIndex] = 1;
        } else if (colIndex == row-(this->nx+2) || colIndex == row+(this->nx+2)) {
            matrix->aValue(i) = 1/hysq;
            this->reconMatrix[row][colIndex] = 1;
        } else {
            assert(false);
        }
    }
}

void PressureSolver::discE(int xInd, int yInd, int row, MatrixIter *matrix) {
    int colIndex;
    double hxsq = simutils::square(this->hx);
    double hysq = simutils::square(this->hy);
    for (int i = matrix->rowBegin(row); i < matrix->rowEndPlusOne(row); i++) {
        colIndex = matrix->getColIndex(i);

        if (colIndex == row) {
            matrix->aValue(i) = -(1/hxsq + 2/hysq);
            this->reconMatrix[row][colIndex]= -3;
        } else if (colIndex == row - 1) {
            matrix->aValue(i) = 1/hxsq;
            this->reconMatrix[row][colIndex] =1;
        } else if (colIndex == row - this->nx || colIndex == row + this->nx) {
            matrix->aValue(i) = 1/hysq;
            this->reconMatrix[row][colIndex] =1;
        } else {
            assert(false);
        }
    }
}

void PressureSolver::discNW(int xInd, int yInd, int row, MatrixIter *matrix) {
    int colIndex;
    double hxsq = simutils::square(this->hx);
    double hysq = simutils::square(this->hy);
    for (int i = matrix->rowBegin(row); i < matrix->rowEndPlusOne(row); i++) {
        colIndex = matrix->getColIndex(i);

        if (colIndex == row) {
            matrix->aValue(i) = -(1/hxsq + 1/hysq);
            this->reconMatrix[row][colIndex] = -2;
        } else if (colIndex == row + 1) {
            matrix->aValue(i) = 1/hxsq;
            this->reconMatrix[row][colIndex] = 1;
        } else if (colIndex == row - this->nx) {
            matrix->aValue(i) = 1/hysq;
            this->reconMatrix[row][colIndex] = 1;
        } else {
            assert(false);
        }
    }
}

void PressureSolver::discN(int xInd, int yInd, int row, MatrixIter *matrix) {
    int colIndex;
    double hxsq = simutils::square(this->hx);
    double hysq = simutils::square(this->hy);
    for (int i = matrix->rowBegin(row); i < matrix->rowEndPlusOne(row); i++) {
        colIndex = matrix->getColIndex(i);

        if (colIndex == row) {
            matrix->aValue(i) = -(2/hxsq + 1/hysq);
            this->reconMatrix[row][colIndex] = -3;
        } else if (colIndex == row - 1 || colIndex == row+1) {
            matrix->aValue(i) = 1/hxsq;
            this->reconMatrix[row][colIndex] = 1;
        } else if (colIndex == row - nx) {
            matrix->aValue(i) = 1/hysq;
            this->reconMatrix[row][colIndex] = 1;
        } else {
            assert(false);
        }
    }
}

void PressureSolver::discNE(int xInd, int yInd, int row, MatrixIter *matrix) {
    int colIndex;
    double hxsq = simutils::square(this->hx);
    double hysq = simutils::square(this->hy);
    for (int i = matrix->rowBegin(row); i < matrix->rowEndPlusOne(row); i++) {
        colIndex = matrix->getColIndex(i);

        if (colIndex == row) {
            matrix->aValue(i) = -(1/hxsq + 1/hysq);
            this->reconMatrix[row][colIndex] = -2;
        } else if (colIndex == row - 1) {
            matrix->aValue(i) = 1/hxsq;
            this->reconMatrix[row][colIndex] = 1;
        } else if (colIndex == row - nx) {
            matrix->aValue(i) = 1/hysq;
            this->reconMatrix[row][colIndex] =1;
        } else {
            assert(false);
        }
    }
}

/**
 * Keeping for testing and simplicity.
*/
// void PressureSolver::solveSOR(int niter, double omega, double dt, int nx, double *x,
//                                 int ny, double *y, double **FU, double **FV, double **p) {
//     int i, j, k;

//     double dx = x[1] - x[0];
//     double dy = y[1] - y[0];

//     double dxsq = simutils::square(dx);
//     double dysq = simutils::square(dy);

//     double rhs;

//     double res;
//     double res_ij;
//     double eW, eE, eS, eN;

//     // Apply the boundary conditions to the pressure at this iteration
//     for (j = 1; j <= ny; j++) {
//         p[j][0] = p[j][1];
//         p[j][nx+1] = p[j][nx];
//     }
//     for (i = 1; i <= nx; i++) {
//         p[0][i] = p[1][i];
//         p[ny+1][i] = p[ny][i];
//     }

//     for (k = 0; k < niter; k++) {
//         // Set the residual to 0 initially
//         res = 0.0;

//         // SOR loop
//         for (j = 1; j <= ny; j++) {
//             for (i = 1; i <= nx; i++) {
//                 // The epsilon functions
//                 // eW = (i == 1) ? 0.0 : 1.0;
//                 // eE = (i == nx) ? 0.0 : 1.0;
//                 // eS = (j == 1) ? 0.0 : 1.0;
//                 // eN = (i == ny) ? 0.0 : 1.0;

//                 // Compute residual on the current step
//                 // Very simple derivative calculation... doing everything by the book
//                 // rhs = ((FU[j][i] - FU[j][i-1])/(dx) + (FV[j][i]-FV[j-1][i])/(dy))/dt;
//                 rhs = 0.0;

//                 // res_ij = (eE*(p[j][i+1]-p[j][i]) - eW*(p[j][i] - p[j][i-1]))/(simutils::square(dx))
//                 //         + (eN*(p[j+1][i] - p[j][i]) - eS*(p[j][i] - p[j-1][i]))/(simutils::square(dy)) - rhs;
//                 res_ij = (p[j][i+1] - 2.0*p[j][i] + p[j][i-1])/dxsq
//                         + (p[j+1][i] - 2.0*p[j][i] + p[j-1][i])/dysq - rhs;

//                 // Update the pressure
//                 // p[j][i] = (1 - omega)*p[j][i] + omega/((eE+eW)/(simutils::square(dx)) + (eN+eS)/(simutils::square(dy)))*res_ij;
//                 p[j][i] = (1 - omega)*p[j][i] + (omega/(2.0/dxsq + 2.0/dysq))*res_ij;

//                 // Accumulate the residual for this iteration
//                 res += simutils::square(res_ij);
//             }
//         }

//         // If the norm residual meets the tolerance, break
//         // Not expecting much in terms of accuracy.
//         if (sqrt(res/(nx*ny)) < 1e-6) {
//             break;
//         }

//         // Apply the boundary conditions to the pressure at this iteration
//         for (j = 1; j <= ny; j++) {
//             p[j][0] = p[j][1];
//             p[j][nx+1] = p[j][nx];
//         }
//         for (i = 1; i <= nx; i++) {
//             p[0][i] = p[1][i];
//             p[ny+1][i] = p[ny][i];
//         }
//     }

//     std::cout << "RHS" << std::endl;
//     for (j = 1; j <= ny; j++) {
//         for (i = 1; i <= nx; i++) {
//             std::cout <<  ((FU[j][i] - FU[j][i-1])/(dx) + (FV[j][i]-FV[j-1][i])/(dy))/dt << ", ";
//         }
//         std::cout << std::endl;
//     }

//     // std::cout << "Pressure" << std::endl;
//     std::cout << "SOR finished in " << k << "iterations" << std::endl;
//     std::cout << "Res was " << sqrt(res/(nx*ny)) << std::endl;
//     for (j = 1; j <= ny; j++) {
//         for (i = 1; i <= nx; i++) {
//             std::cout << p[j][i] << ", ";
//         }
//         std::cout << std::endl;
//     }
// }