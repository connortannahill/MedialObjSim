#include "Pool3D.h"
#include "Boundary3D.h"
#include "SolidObject3D.h"
#include "../Utils/SimUtilities.h"
#include <iostream>
#include <math.h>
#include <set>
#include <limits>
#include <fstream>
#include "../Utils/Discretizations.h"
#include <assert.h>
#include <utility>
#include <map>
#include <queue>
#include "MassSpring3D.h"

using namespace mass_spring;
using namespace std;

// Infinity 
const double DINFINITY = numeric_limits<double>::infinity();

// Constants
int DOMAIN_UNDESCOVERED_3D = -2;
int DOMAIN_INTERSECTION_3D = -3;
int DOMAIN_FLUID_3D = -1;

/* PRIVATE METHODS */

const int CLOSE = 0;
const int ACCEPTED = 1;
const int FAR = 2;

/**
 * Check that this grid point is valid
*/
bool Pool3D::indInRange(int i, int j, int k) {
    return !(i < 0 || i > nx-1 || j < 0 || j > ny-1 || k < 0 || k > nz-1);
}

/**
 * Find closest distance to the boundary faces of the MSS
*/
void Pool3D::closestBoundaryPnt(int structNum, double inPnt[3], double outPnt[3]) {
    double baryCoords[3];
    double closestDist = INFINITY;
    double d;

    int closestId = -1;

    double baryCoordsClosest[3];
    const int MAXCANDS = 100;
    vector<pair<size_t,double> > ret_matches;
    std::vector<size_t> ret_index(MAXCANDS);
    std::vector<double> out_dist_sqr(MAXCANDS);

    // Candidate list
    int numFound = kdTree->knnSearch(&inPnt[0],
                    MAXCANDS, &ret_index[0], &out_dist_sqr[0]);
    
    assert(numFound > 0);
    
    for (int i = 0; i < numFound; i++) {
        mass_spring::massPoint3D pnt = *kdPointCloud->points->at(ret_index.at(i));

        // Must correspond to this one
        if (pnt.structNum != structNum) {
            continue;
        }

        // Iterate through faces
        for (auto face = pnt.faceIds.begin(); face != pnt.faceIds.end(); ++face) {
            face3D faceLoc = solids->at(structNum).faceList->at(*face);
            d = solids->at(structNum).projTriangleDist(inPnt, faceLoc.pntIds[0], faceLoc.pntIds[1], faceLoc.pntIds[2], baryCoords);

            if (d < closestDist) {
                // Record the closest distance
                closestDist = d;
                closestId = *face;

                simutils::copyVals(3, baryCoords, baryCoordsClosest);
            }

        }
    }

    // for (auto face = faceList->begin(); face != faceList->end(); ++face) {
    //     // Get the distance of the input point to this triangle (baryCoords are irrelevant)
    //     d = projTriangleDist(inPnt, face->pntIds[0], face->pntIds[1], face->pntIds[2], baryCoords);

    //     if (d < closestDist) {
    //         // Record the closest distance
    //         closestDist = d;

    //         // Record the nearest IDs and coordinates
    //         closestId = face - faceList->begin();

    //         simutils::copyVals(3, baryCoords, baryCoordsClosest);
    //     }
    // }

    // Make the output point the Barycentric reconstruction of the closest point recorded
    for (int j = 0; j < 3; j++) {
        outPnt[j] = 0.0;
    }

    double pnt[3];
    for (int i = 0; i < 3; i++) {
        pnt[0] = solids->at(structNum).pntList->at(solids->at(structNum).faceList->at(closestId).pntIds[i]).x;
        pnt[1] = solids->at(structNum).pntList->at(solids->at(structNum).faceList->at(closestId).pntIds[i]).y;
        pnt[2] = solids->at(structNum).pntList->at(solids->at(structNum).faceList->at(closestId).pntIds[i]).z;

        for (int j = 0; j < 3; j++) {
            outPnt[j] += baryCoordsClosest[i]*pnt[j];
        }
    }
}

/**
 * Find closest distance to the boundary faces of the MSS
*/
double Pool3D::closestBoundaryDist(int structNum, double inPnt[3]) {
    double baryCoords[3];
    double closestDist = INFINITY;
    double d;

    const int MAXCANDS = 10;
    vector<pair<size_t,double> > ret_matches;
    std::vector<size_t> ret_index(MAXCANDS);
    std::vector<double> out_dist_sqr(MAXCANDS);

    // Candidate list
    int numFound = kdTree->knnSearch(&inPnt[0],
                    MAXCANDS, &ret_index[0], &out_dist_sqr[0]);
    
    assert(numFound > 0);
    
    for (int i = 0; i < numFound; i++) {
        mass_spring::massPoint3D pnt = *kdPointCloud->points->at(ret_index.at(i));

        // Must correspond to this one
        if (pnt.structNum != structNum) {
            continue;
        }

        // Iterate through faces
        for (auto face = pnt.faceIds.begin(); face != pnt.faceIds.end(); ++face) {
            face3D faceLoc = solids->at(structNum).faceList->at(*face);
            d = solids->at(structNum).projTriangleDist(inPnt, faceLoc.pntIds[0], faceLoc.pntIds[1], faceLoc.pntIds[2], baryCoords);

            if (d < closestDist) {
                // Record the closest distance
                closestDist = d;
            }

        }
    }


    // for (auto face = faceList->begin(); face != faceList->end(); ++face) {
    //     // Get the distance of the input point to this triangle (baryCoords are irrelevant)
    //     d = projTriangleDist(inPnt, face->pntIds[0], face->pntIds[1], face->pntIds[2], baryCoords);

    //     if (d < closestDist) {
    //         // Record the closest distance
    //         closestDist = d;
    //     }
    // }
    return closestDist;
}

/**
 * TODO: the phi update rules may break if there is a negative sqrt. The wiki page provides what to
 * do in this case. Should be implemneted when we do non-constant velocity extrap.
 * 
*/
void Pool3D::fastMarchSetNVal(int i, int j, int k, bool nExtrap, int mode) {
    double hx = x[1] - x[0];
    double hy = y[1] - y[0];
    double hz = z[1] - z[0];
    int mo = methodOrd;

    // Ensure this point is in bounds
    if (!indInRange(i, j, k)) {
        return;
    }

    // If this point is in far set it to close, else return
    if (fastMarchingState[k][j][i] == FAR) {
        fastMarchingState[k][j][i] = CLOSE;
    } else {
        return;
    }

    // Use upwinding in each axis to determine upwinding direction for this point,
    // as well as the point value to be used in solving the Eikonal equation
    int upwindX = 0;
    int upwindY = 0;
    int upwindZ = 0;

    // x
    if (i+1 > nx-1 || i-1 < 0) {
        // First checks: the point must be in bounds. If one of the points is OOB but the
        // other is in far, we can't use this direction.
        upwindX = (i+1 > nx-1) ? -1 : 1;
        if (fastMarchingState[k][j][i+upwindX] == FAR) {
            upwindX = 0;
        }
    } else if (fastMarchingState[k][j][i+1] == FAR && fastMarchingState[k][j][i-1] == FAR) {
        // Secondly, if both of the points are in FAR, we can't use this direction
        upwindX = 0;
    } else if (fastMarchingState[k][j][i+1] != FAR && fastMarchingState[k][j][i-1] != FAR) {
        // Next, if we have both directions, find the proper upwind direction
        upwindX = (phiReInit[mo+k][mo+j][mo+i-1] < phiReInit[mo+k][mo+j][mo+i+1]) ? -1 : 1;

        // If the upwind direction is an interface, this will always be the upwind direction
        if (oneGridFromInterfaceStructure(i, j, k)) {
            if (isInterface(objAtIndex(i-1, j, k))) {
                upwindX = -1;
            } else if (isInterface(objAtIndex(i+1, j, k))) {
                upwindX = 1;
            }
        }
    } else {
        // Finally, if there is only one direction, and there are no issues with bounds, we choose this one.
        upwindX = (fastMarchingState[k][j][i+1] != FAR) ? 1 : -1;
    }

    // y
    if (j+1 > ny-1 || j-1 < 0) {
        upwindY = (j+1 > ny-1) ? -1 : 1;
        if (fastMarchingState[k][j+upwindY][i] == FAR) {
            upwindY = 0;
        }
    } else if (fastMarchingState[k][j+1][i] == FAR && fastMarchingState[k][j-1][i] == FAR) {
        // Secondly, if both of the points are in FAR, we can't use this direction
        upwindY = 0;
    } else if (fastMarchingState[k][j+1][i] != FAR && fastMarchingState[k][j-1][i] != FAR) {
        // Next, if we have both directions, find the proper upwind direction
        upwindY = (phiReInit[mo+k][mo+j-1][mo+i] < phiReInit[mo+k][mo+j+1][mo+i]) ? -1 : 1;

        if (oneGridFromInterfaceStructure(i, j, k)) {
            if (isInterface(objAtIndex(i, j-1, k))) {
                upwindY = -1;
            } else if (isInterface(objAtIndex(i, j+1, k))) {
                upwindY = 1;
            }
        }
    } else {
        // Finally, if there is only one direction, and there are no issues with bounds, we choose this one.
        upwindY = (fastMarchingState[k][j+1][i] != FAR) ? 1 : -1;
    }

    // z
    if (k+1 > nz-1 || k-1 < 0) {
        upwindZ = (k+1 > nz-1) ? -1 : 1;
        if (fastMarchingState[k+upwindZ][j][i] == FAR) {
            upwindZ = 0;
        }
    } else if (fastMarchingState[k+1][j][i] == FAR && fastMarchingState[k-1][j][i] == FAR) {
        // Secondly, if both of the points are in FAR, we can't use this direction
        upwindZ = 0;
    } else if (fastMarchingState[k+1][j][i] != FAR && fastMarchingState[k-1][j][i] != FAR) {
        // Next, if we have both directions, find the proper upwind direction
        upwindZ = (phiReInit[mo+k-1][mo+j][mo+i] < phiReInit[mo+k+1][mo+j][mo+i]) ? -1 : 1;
        if (oneGridFromInterfaceStructure(i, j, k)) {
            if (isInterface(objAtIndex(i, j, k-1))) {
                upwindZ = -1;
            } else if (isInterface(objAtIndex(i, j, k+1))) {
                upwindZ = 1;
            }
        }
    } else {
        // Finally, if there is only one direction, and there are no issues with bounds, we choose this one.
        upwindZ = (fastMarchingState[k+1][j][i] != FAR) ? 1 : -1;
    }

    // Now we apply the method outlined in Bridson 2015 to approximate the solution to the
    // reinitialization Eikonal equation
    double d;

    // Case that should be impossible
    assert(!(upwindX == 0 && upwindY == 0 && upwindZ == 0));

    // Build arrays with step sizes
    double phiVals[3] = {
        (upwindX == 0) ? DINFINITY : phiReInit[mo+k][mo+j][mo+i+upwindX],
        (upwindY == 0) ? DINFINITY : phiReInit[mo+k][mo+j+upwindY][mo+i],
        (upwindZ == 0) ? DINFINITY : phiReInit[mo+k+upwindZ][mo+j][mo+i]
    };

    double h[3] = {hx, hy, hz};

    double uVals[3];
    double vVals[3];
    double wVals[3];

    if (nExtrap) {
        uVals[0] = (upwindX == 0) ? DINFINITY : poolU[mo+k][mo+j][mo+i+upwindX];
        uVals[1] = (upwindY == 0) ? DINFINITY : poolU[mo+k][mo+j+upwindY][mo+i];
        uVals[2] = (upwindZ == 0) ? DINFINITY : poolU[mo+k+upwindZ][mo+j][mo+i];

        vVals[0] = (upwindX == 0) ? DINFINITY : poolV[mo+k][mo+j][mo+i+upwindX];
        vVals[1] = (upwindY == 0) ? DINFINITY : poolV[mo+k][mo+j+upwindY][mo+i];
        vVals[2] = (upwindZ == 0) ? DINFINITY : poolV[mo+k+upwindZ][mo+j][mo+i];

        wVals[0] = (upwindX == 0) ? DINFINITY : poolW[mo+k][mo+j][mo+i+upwindX];
        wVals[1] = (upwindY == 0) ? DINFINITY : poolW[mo+k][mo+j+upwindY][mo+i];
        wVals[2] = (upwindZ == 0) ? DINFINITY : poolW[mo+k+upwindZ][mo+j][mo+i];
    }

    // Simultaneous selection sort of all of the above arrays in increasing phi vals
    double temp;
    for (int j = 0; j < 2; j++) {
        for (int i = j+1; i < 3; i++) {
            if (phiVals[i] < phiVals[j]) {
                // Swap phi
                temp = phiVals[j];
                phiVals[j] = phiVals[i];
                phiVals[i] = temp;

                // Swap h
                temp = h[j];
                h[j] = h[i];
                h[i] = temp;

                if (nExtrap) {
                    // Swap u vals
                    temp = uVals[j];
                    uVals[j] = uVals[i];
                    uVals[i] = temp;

                    // Swap v vals
                    temp = vVals[j];
                    vVals[j] = vVals[i];
                    vVals[i] = temp;

                    // Swap w vals
                    temp = wVals[j];
                    wVals[j] = wVals[i];
                    wVals[i] = temp;
                }
            }
        }
    }

    assert(phiVals[0] <= phiVals[1] && phiVals[1] <= phiVals[2]
        && phiVals[0] <= phiVals[2]);

    // Now we apply the Bridson algorithm

    // First, try the closest neighbour
    d = phiVals[0] + h[0];

    // Try the two closest neighbours
    double h0sq = simutils::square(h[0]);
    double h1sq = simutils::square(h[1]);
    double h2sq = simutils::square(h[2]);
    double a, b, c;

    if (d > phiVals[1]) {

        d = 1.0/(2.0*(h0sq+h1sq)) * ((2.0*h1sq*phiVals[0]+2.0*h0sq*phiVals[1])
                    + 2.0*sqrt(simutils::square(h1sq*phiVals[0]+h0sq*phiVals[1])
                        - (h0sq+h1sq)*(h1sq*simutils::square(phiVals[0])
                        + h0sq*simutils::square(phiVals[1]) - h1sq*h0sq)));
        
        if (d > phiVals[2]) {
            // Solve the quadratic for the distance
            a = (h1sq*h2sq + h0sq*h2sq + h0sq*h1sq);
            b = -2.0*(h1sq*h2sq*phiVals[0] + h0sq*h2sq*phiVals[1] + h0sq*h1sq*phiVals[2]);
            c = h1sq*h2sq*simutils::square(phiVals[0]) + h0sq*h2sq*simutils::square(phiVals[1])
                + h0sq*h1sq*simutils::square(phiVals[2]) - h0sq*h1sq*h2sq;
            
            d = (-b + sqrt(simutils::square(b) - 4.0*a*c))/(2.0*a);

            // Three point formula for the velocities
            if (nExtrap) {
                double denom = h1sq*h2sq*(d - phiVals[0]) + h0sq*h2sq*(d - phiVals[1]) + h0sq*h1sq*(d - phiVals[2]);

                poolU[mo+k][mo+j][mo+i] = (h1sq*h2sq*uVals[0]*(d - phiVals[0])
                    + h0sq*h2sq*uVals[1]*(d - phiVals[1]) + h0sq*h1sq*uVals[2]*(d - phiVals[2]))/denom;
                poolV[mo+k][mo+j][mo+i] = (h1sq*h2sq*vVals[0]*(d - phiVals[0])
                    + h0sq*h2sq*vVals[1]*(d - phiVals[1]) + h0sq*h1sq*vVals[2]*(d - phiVals[2]))/denom;
                poolW[mo+k][mo+j][mo+i] = (h1sq*h2sq*wVals[0]*(d - phiVals[0])
                    + h0sq*h2sq*wVals[1]*(d - phiVals[1]) + h0sq*h1sq*wVals[2]*(d - phiVals[2]))/denom;
            }
        } else {
            if (nExtrap) {
                // Two point formula for velocities
                poolU[mo+k][mo+j][mo+i] = (h1sq*uVals[0]*(d - phiVals[0]) + h0sq*uVals[1]
                    * (d - phiVals[1]))/(h1sq*(d - phiVals[0]) + h0sq*(d - phiVals[1]));
                poolV[mo+k][mo+j][mo+i] = (h1sq*vVals[0]*(d - phiVals[0]) + h0sq*vVals[1]
                    * (d - phiVals[1]))/(h1sq*(d - phiVals[0]) + h0sq*(d - phiVals[1]));
                poolW[mo+k][mo+j][mo+i] = (h1sq*wVals[0]*(d - phiVals[0]) + h0sq*wVals[1]
                    * (d - phiVals[1]))/(h1sq*(d - phiVals[0]) + h0sq*(d - phiVals[1]));
            }
        }
    } else {
        if (nExtrap) {
            // One point formula for velocities
            poolU[mo+k][mo+j][mo+i] = uVals[0];
            poolV[mo+k][mo+j][mo+i] = vVals[0];
            poolW[mo+k][mo+j][mo+i] = wVals[0];
        }
    }

    // Assign the distance
    phiReInit[mo+k][mo+j][mo+i] = d;

    if (mode == 2) {
        if (d > collisionDist) {
            fastMarchingState[k][j][i] = ACCEPTED;
            return;
        }
    }

    // Add the point to the heap structure
    heap.push(GridVal3D(i, j, k, phiReInit[mo+k][mo+j][mo+i]));
}

/**
 * Assign the marching state of the current node.
*/
void Pool3D::assignDomainMemberships(int i, int j, int k, double val, int mode) {
    int curMembership = domainMembership(i, j, k);
    int neighMembership;

    int iList[6] = {i,   i,   i-1, i+1, i,   i};
    int jList[6] = {j,   j,   j,   j,   j+1, j-1};
    int kList[6] = {k+1, k-1, k,   k,   k,   k};

    double medX = 0;
    double medY = 0;
    double medZ = 0;

    for (int l = 0; l < 6; l++) {
        int ni = iList[l];
        int nj = jList[l];
        int nk = kList[l];

        if (indInRange(ni, nj, nk)) {
            neighMembership = domainMembership(ni, nj, nk);

            if (neighMembership == DOMAIN_FLUID_3D) {
                domainTracker[nk][nj][ni] = curMembership;
            } else if (neighMembership != curMembership && domainTracker[nk][nj][ni] != DOMAIN_INTERSECTION_3D) {

                // Keep track of the number if intersections about this point. We accumulate
                // them and take the average.

                int offX = (ni - i > 0) ? 1 : 0;
                int offY = (nj - j > 0) ? 1 : 0;
                int offZ = (nk - k > 0) ? 1 : 0;

                medX = x[i+offX];
                medY = y[i+offY];
                medZ = z[i+offZ];

                if (mode == 2 && val < collisionDist) {
                    medialAxisPnts->push_back(make_tuple(medX, medY, medZ));
                }

                domainTracker[k+offZ][j+offY][i+offX] = DOMAIN_INTERSECTION_3D;
            }
        }
    }
}

/**
 * Apply the fast marching method from Adalsteinsson and Sethian 1999
 * to extrapolate the velocity field normal to the interface. This strategy
 * generates a velocity field which should preserve signed distances.
 * 
 * Assumes that the current phi approximates a signed distance function. That domain
 * membership has been established and center of mass velocity computer (only rigid body right now)
 * 
 * mode controls what the fast marching algorithm does
 *  mode = 0: the interpolation based reinititialization
 *  mode = 1: use of the previous isocontour
 *  mode = 2: collision plane detectio
*/
void Pool3D::fastMarch(bool nExtrap, int mode) {
    // Identify all of the points on the boundaries and initialize the heap with their signed distance
    // function values.
    int mo = methodOrd;

    if (mode == 2) {
        medialAxisPnts->clear();
    }

    // Initialize the velocity fields and the temp signed distance function to infinity (negative
    // infinity in the interior of the domain.)

    double phiVal = 0.0;
    objects::FSIObject obj;
    double pnt[3];
    double vTemp[3];

    for (int k = 0; k < this->nz; k++) {
        pnt[2] = simutils::midpoint(z[k], z[k+1]);
        for (int j = 0; j < this->ny; j++) {
            pnt[1] = simutils::midpoint(y[j], y[j+1]);
            for (int i = 0; i < this->nx; i++) {
                pnt[0] = simutils::midpoint(x[i], x[i+1]);

                phiVal = INFINITY;
                obj = this->objAtIndex(i, j, k);

                // Set the default values in the fluid region.
                if (!isInterface(obj)) {
                    if (nExtrap) {
                        this->poolU[mo+k][mo+j][mo+i] = DINFINITY;
                        this->poolV[mo+k][mo+j][mo+i] = DINFINITY;
                        this->poolW[mo+k][mo+j][mo+i] = DINFINITY;
                    }
                    this->phiReInit[mo+k][mo+j][mo+i] = DINFINITY;
                } else {
                    // If an interface cell, set phi to distance from nearest interface
                    double phiVal;
                    if (mode == 0) {
                        double sign = simutils::sign(phi[mo+k][mo+j][mo+i]);
                        int structNum = domainMembership(i, j, k);
                        phiVal = sign*closestBoundaryDist(structNum, pnt);
                    } else {
                        phiVal = phi[mo+k][mo+j][mo+i];
                    }
                    this->phiReInit[mo+k][mo+j][mo+i] = phiVal;
                    
                    if (nExtrap) {
                        solids->at(domainMembership(i, j, k)).interpFaceVels(pnt, vTemp);

                        this->poolU[mo+k][mo+j][mo+i] = vTemp[0];
                        this->poolV[mo+k][mo+j][mo+i] = vTemp[1];
                        this->poolW[mo+k][mo+j][mo+i] = vTemp[2];
                    }
                }

                // Label each of the points
                fastMarchingState[k][j][i] = (obj == objects::FLUID_C) ? FAR : ACCEPTED;

                // Add the value of the curent level set (assumed SDF) at the interfaces.
                if (this->isInterface(obj)) {
                    fastMarchingState[k][j][i] = CLOSE;
                    heap.push(GridVal3D(i, j, k, phiVal));
                }
            }
        }
    }

    // Now run the FM algorithm in the fluid region.
    int i, j, k;
    double val;
    while (!heap.empty()) {
        // Get and pop the top value 
        GridVal3D pnt = heap.top();
        heap.pop();

        // Extract the point information
        i = pnt.getX();
        j = pnt.getY();
        k = pnt.getZ();
        val = pnt.getVal();

        // Assign the domain memberships of the neighbouring points
        if (domainTracker[k][j][i] != DOMAIN_INTERSECTION_3D) {
            assignDomainMemberships(i, j, k, val, mode);
        }

        // Update all of the neighbours according to the fast marching algorithm
        // and add them to the heap structure
        fastMarchSetNVal(i+1, j, k, nExtrap, mode);
        fastMarchSetNVal(i-1, j, k, nExtrap, mode);
        fastMarchSetNVal(i, j+1, k, nExtrap, mode);
        fastMarchSetNVal(i, j-1, k, nExtrap, mode);
        fastMarchSetNVal(i, j, k+1, nExtrap, mode);
        fastMarchSetNVal(i, j, k-1, nExtrap, mode);

        fastMarchingState[k][j][i] = ACCEPTED;
    }

    if (mode == 2) {
        cout << "MODE 2 RUNNING... FOUND " << medialAxisPnts->size() << endl;
        return;
    }

    if (mode != 0) {
        for (int k = 0; k < this->nz; k++) {
            for (int j = 0; j < this->ny; j++) {
                for (int i = 0; i < this->nx; i++) {
                    obj = objAtIndex(i, j, k);

                    fastMarchingState[k][j][i] = (isInterface(obj)) ? CLOSE : FAR;
                }
            }
        }
    }

    double velTemp[3];

    // Now, run the algorithm in the structure domain. Note that at this point,
    // all of the fluid and inteface regions are in the "accepted class" at this point.
    for (int k = 0; k < this->nz; k++) {
        pnt[2] = simutils::midpoint(z[k], z[k+1]);
        for (int j = 0; j < this->ny; j++) {
            pnt[1] = simutils::midpoint(y[j], y[j+1]);
            for (int i = 0; i < this->nx; i++) {
                pnt[0] = simutils::midpoint(x[i], x[i+1]);

                // We only act on structure grids here
                obj = this->objAtIndex(i, j, k);

                if (obj == objects::STRUCTURE) {
                    phiVal = INFINITY;
                    assert(phiVal >= 0);

                    if (this->oneGridFromInterfaceStructure(i, j, k)) {
                        if (mode == 0) {
                            int structNum = domainMembership(i, j, k);
                            phiVal = closestBoundaryDist(structNum, pnt);

                            // Normal extrapolation (very approximate)
                            if (nExtrap) {
                                solids->at(domainMembership(i, j, k)).interpFaceVels(pnt, velTemp);
                                this->poolU[mo+k][mo+j][mo+i] = velTemp[0];
                                this->poolV[mo+k][mo+j][mo+i] = velTemp[1];
                                this->poolW[mo+k][mo+j][mo+i] = velTemp[2];
                            }

                            // Put this in the close field and push it onto the heap
                            heap.push(GridVal3D(i, j, k, phiVal));
                        } else {
                            fastMarchSetNVal(i, j, k, nExtrap, phiVal);
                        }
                    } else {
                        if (nExtrap) {
                            // Otherwise, set the velocity to infinity to make errors easy to spot.
                            this->poolU[mo+k][mo+j][mo+i] = DINFINITY;
                            this->poolV[mo+k][mo+j][mo+i] = DINFINITY;
                            this->poolW[mo+k][mo+j][mo+i] = DINFINITY;
                        }
                    }
                }

                if (mode == 0) {
                    // Seperate the logic for assigning the marching state for brain reasons
                    if (obj == objects::STRUCTURE && this->oneGridFromInterfaceStructure(i, j, k)) {
                        fastMarchingState[k][j][i] = CLOSE;
                        phiReInit[mo+k][mo+j][mo+i] = phiVal;
                    } else if (obj == objects::STRUCTURE) {
                        fastMarchingState[k][j][i] = FAR;
                        phiReInit[mo+k][mo+j][mo+i] = DINFINITY;
                    } else {
                        fastMarchingState[k][j][i] = ACCEPTED;
                    }
                }
            }
        }
    }

    if (mode != 0) {
        for (int k = 0; k < this->nz; k++) {
            for (int j = 0; j < this->ny; j++) {
                for (int i = 0; i < this->nx; i++) {
                    obj = objAtIndex(i, j, k);
                    if (obj == objects::STRUCTURE && this->oneGridFromInterfaceStructure(i, j, k)) {
                        fastMarchingState[k][j][i] = CLOSE;
                    } else if (obj == objects::STRUCTURE) {
                        fastMarchingState[k][j][i] = FAR;
                        phiReInit[mo+k][mo+j][mo+i] = DINFINITY;
                    } else {
                        fastMarchingState[k][j][i] = ACCEPTED;
                    }
                }
            }
        }
    }

    // Now, run fast marching on the interior data.
    while (!heap.empty()) {
        GridVal3D pnt = heap.top();
        heap.pop();

        // Extract the point information
        i = pnt.getX();
        j = pnt.getY();
        k = pnt.getZ();
        val = pnt.getVal();

        // Update all of the neighbours according to the fast marching algorithm
        // and add them to the heap structure
        fastMarchSetNVal(i+1, j, k, nExtrap, mode);
        fastMarchSetNVal(i-1, j, k, nExtrap, mode);
        fastMarchSetNVal(i, j+1, k, nExtrap, mode);
        fastMarchSetNVal(i, j-1, k, nExtrap, mode);
        fastMarchSetNVal(i, j, k+1, nExtrap, mode);
        fastMarchSetNVal(i, j, k-1, nExtrap, mode);

        fastMarchingState[k][j][i] = ACCEPTED;
    }

    // Finally, flip the sign of the phi value in all of structure points.
    for (int k = 0; k < this->nz; k++) {
        for (int j = 0; j < this->ny; j++) {
            for (int i = 0; i < this->nx; i++) {
                if (objAtIndex(i, j, k) == objects::STRUCTURE) {
                    phiReInit[k+mo][j+mo][i+mo] *= -1;
                }
            }
        }
    }
}


/**
 * Algorithm to detect collisions and setup data structures for possible force calculations
*/
void Pool3D::detectCollisions() {
    // Wipe the collision memory from the MSS's
    for (auto mss = solids->begin(); mss != solids->end(); ++mss) {

        // For each node, clear the set in the node collisions
        for (auto nodeSet = mss->nodeCols->begin(); nodeSet != mss->nodeCols->end(); ++nodeSet) {
            nodeSet->clear();
        }
    }

    // Now, check the medial axis for any possible collisions using thresholding on phi.
    double medX, medY, medZ, medPhi;
    int xCell, yCell, zCell;
    set<int> allCollisions;
    set<int> collidingIds;
    vector<tuple<double,double,double>> medialAxisCollisionPnts;

    for (auto tup = medialAxisPnts->begin(); tup != medialAxisPnts->end(); ++tup) {
        // Find the cell location of the medial axis point
        medX = get<0>(*tup);
        medY = get<1>(*tup);
        medZ = get<2>(*tup);

        xCell = simutils::findLimInfMeshPoint(medX, this->x, this->nx+1);
        yCell = simutils::findLimInfMeshPoint(medY, this->y, this->ny+1);
        zCell = simutils::findLimInfMeshPoint(medZ, this->z, this->nz+1);

        // Find the value of the level set function in this cell
        medPhi = interpolatePhi(medX, medY, medZ);

        // If this value is beneath the threshold for possible collision, look at which objects are colliding
        if (medPhi < collisionDist) {
            // Build up a list of the neighbours around this cell to find which unique colliding nodes
            for (int k = zCell-2; k <= zCell+2; k++) {
                for (int j = yCell-2; j <= yCell+2; j++) {
                    for (int i = xCell-2; i <= xCell+2; i++) {
                        if (indInRange(i, j, k) && domainMembership(i, j, k) != DOMAIN_INTERSECTION_3D) {
                            collidingIds.insert(domainMembership(i, j, k));
                            allCollisions.insert(domainMembership(i, j, k));
                        }
                    }
                }

            }

            tuple<double,double,double> medLoc(medX, medY, medZ);
            medialAxisCollisionPnts.push_back(medLoc);
        }

        collidingIds.clear();
    }

    // If we are using proper repulsive forces, we now build up the KDTree and
    // all of the node collections
    if (repulseMode == 2 && medialAxisCollisionPnts.size() > 0) {

        // Place all of the colliding MSS's into the KDTree
        // kdPointCloud->resetCloud();
        // for (auto colMss = allCollisions.begin(); colMss != allCollisions.end(); ++colMss) {
        //     kdPointCloud->addMSS(solids->at(*colMss));
        // }

        // // Now rebuild the kdtree for efficient searching.
        // kdTree->buildIndex();
        // vector<pair<uint32_t,double> > ret_matches;
        vector<pair<size_t,double> > ret_matches;
        nanoflann::SearchParams params;
        set<massPoint3D*> allCols; // TESTING: build a set of the detected collision nodes and output them to a file for viz.
        double SCAL_FAC = 1.0;
        int nMatches;
        for (auto tup = medialAxisCollisionPnts.begin(); tup != medialAxisCollisionPnts.end(); ++tup) {
            medX = get<0>(*tup);
            medY = get<1>(*tup);
            medZ = get<2>(*tup);

            double queryPnt[3] = {medX, medY, medZ};

            // Do a radius search around this medial axis cell to discover all of the
            // potentially interacting nodes
            nMatches = kdTree->radiusSearch(&queryPnt[0],
                simutils::square(SCAL_FAC*collisionDist), ret_matches, params);

            for (int i = 0; i < nMatches; i++) {
                allCols.insert(kdPointCloud->points->at(ret_matches.at(i).first));
            }
            
            // For each of the matches, we keep track in each MSS of which of the points it is potentially colliding with
            vector<int> nearestMss;
            vector<int> nearestPointMatch;
            int id = 0;
            for (auto match = ret_matches.begin(); match != ret_matches.end(); ++match) {
                // Get the id within the MSS and the MSS itself
                int mssId = kdPointCloud->points->at(match->first)->structNum;

                // If this is a new Id, since the array is ordered according to proximity to m, we add it
                if (std::find(nearestMss.begin(), nearestMss.end(), mssId) == nearestMss.end()) {
                    nearestMss.push_back(mssId);
                    nearestPointMatch.push_back(id);
                }
                id++;
            }

            // if (nearestMss.size() < 2) {
            //     cout << "ERROR IN FINDING MSS NODES AT COLLISION POINT" << endl;
            //     assert(nearestMss.size() < 2);
            // }

            // For each of the matches, we keep track in each MSS of which of the points
            // it is potentially colliding with
            for (auto match = ret_matches.begin(); match != ret_matches.end(); ++match) {
                // Get the id within the MSS and the MSS itself
                int id = kdPointCloud->points->at(match->first)->nodeId;
                int mssId = kdPointCloud->points->at(match->first)->structNum;

                // For this node, place every collision node nearby into the node tracking data structure
                for (int i = 0; i < nearestMss.size(); i++) {
                    if (nearestMss.at(i) != mssId) {
                        solids->at(mssId).nodeCols->at(id).insert(
                            kdPointCloud->points->at(ret_matches.at(nearestPointMatch.at(i)).first));
                    }
                }
            }

            ret_matches.clear();
        }
    }
}

void Pool3D::create3DPool(Boundary3D &boundary,
                    vector<SolidObject3D> &structures,
                    SimParams3D &params) {

    // Get the fixed time step
    if (params.dtFixSet) {
        this->dtFix = params.dtFix;
    } else {
        this->dtFix = -1;
    }

    // Set the body forces
    this->gx = params.gx;
    this->gy = params.gy;
    this->gz = params.gz;
    
    // Set the method order used for the evolution of the level set equation
    // TODO: generalize this
    this->methodOrd = 3;

    // Set the repulsion mode and distance
    this->collisionStiffness = params.collisionStiffness;
    this->collisionDist = params.collisionDist;
    this->repulseDist = 2.0 * this->collisionDist;
    this->repulseMode = params.repulseMode;

    // Create the kdTree data object and point cloud
    kdPointCloud = new KDPointCloud3D();
    kdTree = new KDTree3D(3, *kdPointCloud, KDTreeSingleIndexAdaptorParams(10));

    // Get the fluid viscosity
    this->mu = params.mu;

    // Create the uniform meshes from the boundary object
    this->nx = params.nx;
    this->ny = params.ny;
    this->nz = params.nz;

    // MSS nx, ny, nz
    int mssNx = params.mssNx;
    int mssNy = params.mssNy;
    int mssNz = params.mssNz;

    // Create the arrays for the externally generated velocity field
    this->poolU = simutils::new_constant(nz+2*methodOrd, ny+2*methodOrd,
                                        nx+2*methodOrd, 0.0);
    this->poolV = simutils::new_constant(nz+2*methodOrd, ny+2*methodOrd,
                                        nx+2*methodOrd, 0.0);
    this->poolW = simutils::new_constant(nz+2*methodOrd, ny+2*methodOrd,
                                        nx+2*methodOrd, 0.0);

    this->nStructs = structures.size();

    // Create the pool labels
    this->pool = new objects::FSIObject**[nz];
    for (int k = 0; k < nz; k++) {
        this->pool[k] = new objects::FSIObject*[ny];
        for (int j = 0; j < ny; j++) {
            this->pool[k][j] = new objects::FSIObject[nx];
        }
    }

    // Set all cells to fluids to start.
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                this->pool[k][j][i] = objects::FLUID_C;
            }
        }
    }

    // Generate the 1D meshes
    this->x = boundary.generateXMesh(this->nx);
    this->y = boundary.generateYMesh(this->ny);
    this->z = boundary.generateZMesh(this->nz);

    this->hx = x[1] - x[0];
    this->hy = y[1] - y[0];
    this->hz = z[1] - z[0];

    // Allocate the phi array and temp array for time stepping
    this->phi = simutils::new_constant(nz+2*methodOrd, ny+2*methodOrd,
                                        nx+2*methodOrd, DINFINITY);
    this->phiReInit = simutils::new_constant(nz+2*methodOrd, ny+2*methodOrd,
                                        nx+2*methodOrd, DINFINITY);
    this->phiRk1 = simutils::new_constant(nz+2*methodOrd, ny+2*methodOrd,
                                        nx+2*methodOrd, 0.0);
    this->phiRk2 = simutils::new_constant(nz+2*methodOrd, ny+2*methodOrd,
                                        nx+2*methodOrd, 0.0);
    
    // Apply the boundaries using the implicit SolidObject functions.
    // Also places the tracer particals
    // cout << "Embedding the shapes" << endl;
    if (nStructs > 0) {
        this->tracers = new simstructs::tracer3D[nStructs];
        if (structures.size() != 0) {
            for (int l = 0; l < nStructs; l++) {
                this->embedShape(structures.at(l), l);
            }
        } else {
            assert(false);
        }
    }
    // cout << "FINISHED Embedding the shapes" << endl;

    // Set the initial pool object velocities based on the information from the particles. Important
    // For assigning initial boundary values around the structure.
    if (nStructs == 1) {
        simutils::set_constant(nz+2*methodOrd, ny+2*methodOrd, nx+2*methodOrd, tracers[0].u, this->poolU);
        simutils::set_constant(nz+2*methodOrd, ny+2*methodOrd, nx+2*methodOrd, tracers[0].v, this->poolV);
        simutils::set_constant(nz+2*methodOrd, ny+2*methodOrd, nx+2*methodOrd, tracers[0].w, this->poolW);
    } else if (nStructs == 0) {
        simutils::set_constant(nz+2*methodOrd, ny+2*methodOrd, nx+2*methodOrd, 0.0, this->poolU);
        simutils::set_constant(nz+2*methodOrd, ny+2*methodOrd, nx+2*methodOrd, 0.0, this->poolV);
        simutils::set_constant(nz+2*methodOrd, ny+2*methodOrd, nx+2*methodOrd, 0.0, this->poolW);
    }

    // Create the mass-spring system for each of the objects

    // Enumerate the Pool object
    this->enumeratePool();
    this->enumChanged = true;

    // Set up the domain membership matrix, first allocating it
    domainTracker = simutils::new_constant(nz, ny, nx, DOMAIN_UNDESCOVERED_3D);
    this->setUpDomainArray();

    // Create the medial axis point vector
    medialAxisPnts = new vector<tuple<double,double,double>>();

    // if (isDeformable) {
    solids = new vector<MassSpring3D>();
    if (nx != mssNx || ny != mssNy || nz != mssNz) {
        assert(false);
        // Create different dimension Pool to represent the solid and copy the current Pool.

        SimParams3D temp = params;
        params.setNx(params.mssNx);
        params.setNy(params.mssNy);
        params.setNz(params.mssNz);

        // Create new temp Pool
        Pool3D *poolTemp = new Pool3D(boundary, structures, temp);
        
        // Copy the Mass Spring objects from the old array to the new one
        for (int i = 0; i < nStructs; i++) {
            // cout << "Adding MSS " << i << endl;
            solids->push_back(MassSpring3D((poolTemp->solids)->at(i)));
            // cout << "FINISHED Adding MSS " << i << endl;
        }

        double x, y, z;
        for (int k = 0; k < nz; k++) {
            z = simutils::midpoint(this->z[k], this->z[k+1]);
            for (int j = 0; j < ny; j++) {
                y = simutils::midpoint(this->y[j], this->y[j+1]);
                for (int i = 0; i < nx; i++) {
                    x = simutils::midpoint(this->x[i], this->x[i+1]);
                    phi[this->methodOrd+k][this->methodOrd+j][this->methodOrd+i] = poolTemp->interpolatePhi(x, y, z);
                }
            }
        }
        this->enumeratePool();
        this->setUpDomainArray();

        // Delete the temp Pool
        delete poolTemp;
    } else {
        // Create the mass-spring system for each of the objects
        solids = new vector<MassSpring3D>();

        if (params.updateMode == 2 && !params.dtFixSet) {
            cout << "Attempting to use ADMM without fixed time step!" << endl;
            assert(false);
        }

        // cout << "nStructs = " << nStructs << endl;
        for (int i = 0; i < nStructs; i++) {
            // cout << "i = " << i << endl;
            solids->push_back(MassSpring3D(*this, i, structures.at(i),
                                params.updateMode, params.elementMode));
            if (params.updateMode == 2) {
                // cout << "i = " << i << endl;
                // cout << "size of solids vector " << solids->size() << endl;
                solids->at(i).setAdmmTol(params.admmTol);
                assert(false);
            }
        }
    }

    // Add all of the MSS's
    kdPointCloud->resetCloud();
    for (int i = 0 ; i < nStructs; i++) {
        kdPointCloud->addMSS(solids->at(i));
    }

    // Now rebuild the kdtree for efficient searching.
    kdTree->buildIndex();

    // State of each grid point within the FMM
    fastMarchingState = simutils::new_constant(nz, ny, nz, FAR);
    
    fastMarch(false, 0);
    simutils::copyVals(nx+2*methodOrd, ny+2*methodOrd, nz+2*methodOrd, phiReInit, phi);

    this->enumeratePool();

    if (repulseMode == 2) {
        detectCollisions();
    }
}

/**
 * Interpolate phi
*/
double Pool3D::interpolatePhi(double x, double y, double z) {
    int mo = methodOrd;

    // Find the "lim inf" mesh points for the current tracer location. Needed for all modes
    int i = simutils::findLimInfMeshPoint(x, this->x, this->nx+1);
    int j = simutils::findLimInfMeshPoint(y, this->y, this->ny+1);
    int k = simutils::findLimInfMeshPoint(z, this->z, this->nz+1);

    // Handle the case where the left-most points are outside of the boundary.
    int il;
    int jl;
    int kl;

    if (i == 0) {
        il = i+1;
    } else if (i == nx-1) {
        il = i-1;
    } else {
        il = i;
    }

    if (j == 0) {
        jl = j+1;
    } else if (j == ny-1) {
        jl = j-1;
    } else {
        jl = j;
    }

    if (k == 0) {
        kl = k+1;
    } else if (k == nz-1) {
        kl = k-1;
    } else {
        kl = k;
    }

    double xim = simutils::midpoint(this->x[il-1], this->x[il]);
    double xi = simutils::midpoint(this->x[il], this->x[il+1]);
    double xip = simutils::midpoint(this->x[il+1], this->x[il+2]);

    double yim = simutils::midpoint(this->y[jl-1], this->y[jl]);
    double yi = simutils::midpoint(this->y[jl], this->y[jl+1]);
    double yip = simutils::midpoint(this->y[jl+1], this->y[jl+2]);

    double zim = simutils::midpoint(this->z[kl-1], this->z[kl]);
    double zi = simutils::midpoint(this->z[kl], this->z[kl+1]);
    double zip = simutils::midpoint(this->z[kl+1], this->z[kl+2]);

    double xMesh[2];
    double yMesh[2];
    double zMesh[2];

    int xl, xr, yl, yr, zl, zr;

    if (xim <= x && x <= xi) {
        xMesh[0] = xim;
        xMesh[1] = xi;

        xl = mo+il-1;
        xr = mo+il;
    } else {
        xMesh[0] = xi;
        xMesh[1] = xip;

        xl = mo+il;
        xr = mo+il+1;
    }

    if (yim <= y && y <= yi) {
        yMesh[0] = yim;
        yMesh[1] = yi;

        yl = mo+jl-1;
        yr = mo+jl;
    } else {
        yMesh[0] = yi;
        yMesh[1] = yip;

        yl = mo+jl;
        yr = mo+jl+1;
    }

    if (zim <= z && z <= zi) {
        zMesh[0] = zim;
        zMesh[1] = zi;

        zl = mo+kl-1;
        zr = mo+kl;
    } else {
        zMesh[0] = zi;
        zMesh[1] = zip;

        zl = mo+kl;
        zr = mo+kl+1;
    }

    double vals[8] = {phi[zl][yl][xl], phi[zl][yl][xr],
                            phi[zl][yr][xl], phi[zl][yr][xr],

                      phi[zr][yl][xl], phi[zr][yl][xr],
                            phi[zr][yr][xl], phi[zr][yr][xr]};
    

    // Calcaulte phi with the trilinear interpolation formula
    return simutils::triLinearInterpolation(x, y, z, xMesh, yMesh, zMesh, vals);
}

/**
 * Interpolate the gradient of phi
*/
void Pool3D::interpolatePhiGrad(double x, double y, double z, double phiGrad[3]) {
    int mo = methodOrd;

    // Find the "lim inf" mesh points for the current tracer location. Needed for all modes
    int i = simutils::findLimInfMeshPoint(x, this->x, this->nx+1);
    int j = simutils::findLimInfMeshPoint(y, this->y, this->ny+1);
    int k = simutils::findLimInfMeshPoint(z, this->z, this->nz+1);

    // Handle the case where the left-most points are outside of the boundary.
    int il;
    int jl;
    int kl;

    if (i == 0) {
        il = i+1;
    } else if (i == nx-1) {
        il = i-1;
    } else {
        il = i;
    }

    if (j == 0) {
        jl = j+1;
    } else if (j == ny-1) {
        jl = j-1;
    } else {
        jl = j;
    }

    if (k == 0) {
        kl = k+1;
    } else if (k == nz-1) {
        kl = k-1;
    } else {
        kl = k;
    }

    double xim = simutils::midpoint(this->x[il-1], this->x[il]);
    double xi = simutils::midpoint(this->x[il], this->x[il+1]);
    double xip = simutils::midpoint(this->x[il+1], this->x[il+2]);

    double yim = simutils::midpoint(this->y[jl-1], this->y[jl]);
    double yi = simutils::midpoint(this->y[jl], this->y[jl+1]);
    double yip = simutils::midpoint(this->y[jl+1], this->y[jl+2]);

    double zim = simutils::midpoint(this->z[kl-1], this->z[kl]);
    double zi = simutils::midpoint(this->z[kl], this->z[kl+1]);
    double zip = simutils::midpoint(this->z[kl+1], this->z[kl+2]);

    double xMesh[2];
    double yMesh[2];
    double zMesh[2];

    int xl, xr, yl, yr, zl, zr;

    if (xim <= x && x <= xi) {
        xMesh[0] = xim;
        xMesh[1] = xi;

        xl = mo+il-1;
        xr = mo+il;
    } else {
        xMesh[0] = xi;
        xMesh[1] = xip;

        xl = mo+il;
        xr = mo+il+1;
    }

    if (yim <= y && y <= yi) {
        yMesh[0] = yim;
        yMesh[1] = yi;

        yl = mo+jl-1;
        yr = mo+jl;
    } else {
        yMesh[0] = yi;
        yMesh[1] = yip;

        yl = mo+jl;
        yr = mo+jl+1;
    }

    if (zim <= z && z <= zi) {
        zMesh[0] = zim;
        zMesh[1] = zi;

        zl = mo+kl-1;
        zr = mo+kl;
    } else {
        zMesh[0] = zi;
        zMesh[1] = zip;

        zl = mo+kl;
        zr = mo+kl+1;
    }

    double vals[8] = {phi[zl][yl][xl], phi[zl][yl][xr],
                            phi[zl][yr][xl], phi[zl][yr][xr],

                      phi[zr][yl][xl], phi[zr][yl][xr],
                            phi[zr][yr][xl], phi[zr][yr][xr]};
    

    // Calcaulte phi with the trilinear interpolation formula
    simutils::triLinearGradient(x, y, z, xMesh, yMesh, zMesh, vals, phiGrad);
}

/**
 * Update the position of the tracer particals within each level set function.
 * 
 * Two modes of operation: 
 *    int mode = -1: Update using gradient descent mode.
 *    int mode = 1: Update using interpolation of the local velocity field.
 *    int mode = 0: Bring the tracer to the interface (useful when the tracer moves outside of the level set)
 * 
 * NOTE: all tracer operations and math assumes that ~
 *         - the boundaries of each structure are always perserved (they can not join together)
 *         - phi approximates a signed distance function
 *         - The structure can not leave through the edge of the pool (so no periodic boundary conditions)
 * 
 * NOTE: Prevents the tracer partical from wandering outside of the domain
*/
void Pool3D::updateTracer(int structNum, double dt, int mode) {
    // Extract information from the tracer for this structure
    double tracerX = this->tracers[structNum].x;
    double tracerY = this->tracers[structNum].y;
    double tracerZ = this->tracers[structNum].z;

    // Order of the numerical method, used for offsets
    int mo = methodOrd;

    // Find the "lim inf" mesh points for the current tracer location. Needed for all modes
    int i = simutils::findLimInfMeshPoint(tracers[structNum].x, this->x, this->nx+1);
    int j = simutils::findLimInfMeshPoint(tracers[structNum].y, this->y, this->ny+1);
    int k = simutils::findLimInfMeshPoint(tracers[structNum].z, this->z, this->nz+1);

    /* Noting that phi is defined at the midpoints of cell grids compute the stencil for the bilinear inteprolation.
        Use the right (up) shifted stencil unless the rightmost (topmost) cell is where the tracer lays, in this case we
        use the left-most stencil */
    
    int il = (i == nx-1) ? i-1 : i;
    int jl = (j == ny-1) ? j-1 : j;
    int kl = (k == nz-1) ? k-1 : k;

    // The horizontal (x) stencil
    double xMesh[2] = {simutils::midpoint(this->x[il], this->x[il+1]),
                           simutils::midpoint(this->x[il+1], this->x[il+2])};

    // The vertical (y) stencil
    double yMesh[2] = {simutils::midpoint(this->y[jl], this->y[jl+1]),
                            simutils::midpoint(this->y[jl+1], this->y[jl+2])};

    // The up/down (z) stencil
    double zMesh[2] = {simutils::midpoint(this->z[kl], this->z[kl+1]), 
                            simutils::midpoint(this->z[kl+1], this->z[kl+2])};
    
    double gradPhi_ijk[3];
    double phi_ijk;

    if (mode == 0) {
        /* Jump to interface mode */
        cout << "Interface jumping mode not yet implemented in 3D! Try again later!" << endl;
        assert(false);
    } else if (mode == 1) {
        /* Velocity field jump mode */

        // Use trilinear interpolation to find the values of the velocity field at the tracer partical

        // u:
        double uVals[8] = {poolU[mo+kl][mo+jl][mo+il], poolU[mo+kl][mo+jl][mo+il+1],
                                poolU[mo+kl][mo+jl+1][mo+il], poolU[mo+kl][mo+jl+1][mo+il+1],

                           poolU[mo+kl+1][mo+jl][mo+il], poolU[mo+kl+1][mo+jl][mo+il+1],
                                poolU[mo+kl+1][mo+jl+1][mo+il], poolU[mo+kl+1][mo+jl+1][mo+il+1]};

        double uCur = simutils::triLinearInterpolation(tracerX, tracerY, tracerZ,
                                                         xMesh, yMesh, zMesh, uVals);

        // v:
        double vVals[8] = {poolV[mo+kl][mo+jl][mo+il], poolV[mo+kl][mo+jl][mo+il+1],
                                poolV[mo+kl][mo+jl+1][mo+il], poolV[mo+kl][mo+jl+1][mo+il+1],

                           poolV[mo+kl+1][mo+jl][mo+il], poolV[mo+kl+1][mo+jl][mo+il+1],
                                poolV[mo+kl+1][mo+jl+1][mo+il], poolV[mo+kl+1][mo+jl+1][mo+il+1]};

        double vCur = simutils::triLinearInterpolation(tracerX, tracerY, tracerZ,
                                                         xMesh, yMesh, zMesh, vVals);

        double wVals[8] = {poolW[mo+kl][mo+jl][mo+il], poolW[mo+kl][mo+jl][mo+il+1],
                                poolW[mo+kl][mo+jl+1][mo+il], poolW[mo+kl][mo+jl+1][mo+il+1],

                           poolW[mo+kl+1][mo+jl][mo+il], poolW[mo+kl+1][mo+jl][mo+il+1],
                                poolW[mo+kl+1][mo+jl+1][mo+il], poolW[mo+kl+1][mo+jl+1][mo+il+1]};

        double wCur = simutils::triLinearInterpolation(tracerX, tracerY, tracerZ,
                                                         xMesh, yMesh, zMesh, wVals);

        // Update the position of the tracer partical
        tracers[structNum].x += dt*uCur;
        tracers[structNum].y += dt*vCur;
        tracers[structNum].w += dt*wCur;
    } else if (mode == 2) {
        /* Use gradient descent with line search and the interior perserving step limit
           to update the position of the tracer partical */
        
        this->interpolatePhiGrad(tracerX, tracerY, tracerZ, gradPhi_ijk);
        phi_ijk = this->interpolatePhi(tracerX, tracerY, tracerZ);
        
        
        // Use the closest point result to compute initial stepsize, ensuring that the updated point will be
        // within the negative part of the domain
        double alpha = 0.1;

        // First gradient update
        double xStep = tracers[structNum].x - alpha*gradPhi_ijk[0];
        double yStep = tracers[structNum].y - alpha*gradPhi_ijk[1];
        double zStep = tracers[structNum].z - alpha*gradPhi_ijk[2];

        const int max_iter = 10;
        int iters = 0;

        while (interpolatePhi(xStep, yStep, zStep) > phi_ijk) {
            alpha /= 2;
            
            xStep = tracerX - alpha*gradPhi_ijk[0];
            yStep = tracerY - alpha*gradPhi_ijk[1];
            zStep = tracerZ - alpha*gradPhi_ijk[2];

            iters++;

            if (iters == max_iter) {
                break;
            }
        }

        if (iters < max_iter) {
            tracers[structNum].x = xStep;
            tracers[structNum].y = yStep;
            tracers[structNum].z = yStep;
        }
    } else {
        tracers[structNum].x += dt * tracers[structNum].u;
        tracers[structNum].y += dt * tracers[structNum].v;
        tracers[structNum].z += dt * tracers[structNum].w;
    }
}

/**
 * Embed the implcitly defined shapes into the pool and initalize the tracer
 * particals.
*/
void Pool3D::embedShape(SolidObject3D &struc, int structNum) {
    int i, j, k;
    int mo = this->methodOrd;
    double minValue = phi[mo][mo][mo];
    int min_x, min_y, min_z;

    // Compute the signed distance function
    if (this->nStructs >= 0) {
        for (k = 0; k < this->nz; k++) {
            for (j = 0; j < this->ny; j++) {
                for (i = 0; i < this->nx; i++) {
                    this->phiReInit[k+mo][j+mo][i+mo] =
                        struc.shapeEval(simutils::midpoint(this->x[i],
                            this->x[i+1]), simutils::midpoint(this->y[j], this->y[j+1]), 
                            simutils::midpoint(this->z[k], this->z[k+1]));
                    
                    if (this->phiReInit[k+mo][j+mo][i+mo] < minValue) {
                        minValue = this->phiReInit[k+mo][j+mo][i+mo];
                        min_x = i;
                        min_y = j;
                        min_z = k;
                    }
                }
            }
        }
    }

    if (this->nStructs >= 0) {
        this->tracers[structNum].x = simutils::midpoint(this->x[min_x-1], this->x[min_x]);
        this->tracers[structNum].y = simutils::midpoint(this->y[min_y-1], this->y[min_y]);
        this->tracers[structNum].z = simutils::midpoint(this->z[min_z-1], this->z[min_z]);
        this->tracers[structNum].mass = struc.getMass();
        this->tracers[structNum].u = struc.getU0(); // Velocity of the stucture, useful for update rules
        this->tracers[structNum].v = struc.getV0();
        this->tracers[structNum].w = struc.getW0();
        this->tracers[structNum].isDeformable = struc.getObjType() == SolidObject3D::ObjectType::DEFORMABLE;
        this->isDeformable = this->isDeformable || struc.getObjType() == SolidObject3D::ObjectType::DEFORMABLE;
        
        cout << "tracer for the first (x, y, z) = (" << this->tracers[structNum].x << ", " << this->tracers[structNum].y << ", " << this->tracers[structNum].z << ")" << endl;
    }

    if (structNum == 0) {
        // If this is the first shape being embedded, simply copy into phi.
        simutils::copyVals(this->nx+2*methodOrd, this->ny+2*methodOrd, this->nz+2*methodOrd, phiReInit, phi);
    } else {
        // If this is not the first structure, take the union of the current structure with the existing
        // level set.
        for (k = 0; k < this->nz; k++) {
            for (j = 0; j < this->ny; j++) {
                for (i = 0; i < this->nx; i++) {
                    this->phi[k+mo][j+mo][i+mo] = min( this->phi[k+mo][j+mo][i+mo], this->phiReInit[k+mo][j+mo][i+mo] );
                }
            }
        }
    }
}

/**
 * Label a fluid interface according to the structure. Note that we only allow
 * interfaces that have at least two neighbouring fluid cells. If it does not satisfy this
 * condition, we label it as a structure.
 * 
 * Note interfaces are defined as follows: (plus UP and DOWN values)
 * 
 *   NW      N      NE
 *   |---------------|
 *   |               |
 *   |               |
 * W |       C       | E
 *   |               |
 *   |---------------|
 *   SW      S       SE
 * 
*/
void Pool3D::labelInterface(int i, int j, int k) {
    int n, e, s, w, u, d;

    // Use the definition of the enum (prime factors for N, E, S, W) to
    // quickly decide the interface
    n = ( this->pool[k][j+1][i] != objects::STRUCTURE) ? objects::NORTH : 1;
    e = ( this->pool[k][j][i+1] != objects::STRUCTURE) ? objects::EAST : 1;
    s = ( this->pool[k][j-1][i] != objects::STRUCTURE) ? objects::SOUTH : 1;
    w = ( this->pool[k][j][i-1] != objects::STRUCTURE) ? objects::WEST : 1;
    u = ( this->pool[k+1][j][i] != objects::STRUCTURE) ? objects::UP : 1;
    d = ( this->pool[k-1][j][i] != objects::STRUCTURE) ? objects::DOWN : 1;

    int label = u*d*n*e*s*w;

    // Get rid of the special case where all of the the 6-connected neighbours
    // of this cell are either fluids or interfaces by just assigning it to a fluid
    // cell
    if (!isStructConnected(i, j, k)) {
        label = objects::FLUID_C;
    }

    this->pool[k][j][i] = (objects::FSIObject)label;
}

bool Pool3D::enumHasChanged() {
    return this->enumChanged;
}

/**
 * Enumerate the pool based on the current phi.
 * 
 * Boundaries of the rectangular domain are labelled in the obvious way.
 * All others are labelled adaptively using the signed distance function/
*/
void Pool3D::enumeratePool() {
    int i, j, k;

    int mo = this->methodOrd;

    // If there is not a signed distance function defined, set all internal points to fluids
    if (this->nStructs != 0) {

        for (k = 1; k < nz-1; k++) {
            for (j = 1; j < ny-1; j++) {
                for (i = 1; i < nx-1; i++) {
                    if (this->phi[k+mo][j+mo][i+mo] < 0) {
                        this->enumChanged = (this->pool[k][j][i] != objects::STRUCTURE);
                        this->pool[k][j][i] = objects::STRUCTURE;
                    } else {
                        this->enumChanged = (this->pool[k][j][i] != objects::FLUID_C);
                        this->pool[k][j][i] = objects::FLUID_C;
                    }
                }
            }
        }

        // Find the interface cells
        for (k = 1; k < nz-1; k++) {
            for (j = 1; j < ny-1; j++) {
                for (i = 1; i < nx-1; i++) {
                    if (this->pool[k][j][i] == objects::FLUID_C) {
                        if (this->isStructConnected(i, j, k)) {
                            this->pool[k][j][i] = objects::UNLABELLED_INTERFACE;
                        }
                    }
                }
            }
        }

        // Go through and label the unlabelled interfaces
        for (k = 1; k < nz-1; k++) {
            for (j = 1; j < ny-1; j++) {
                for (i = 1; i < nx-1; i++) {
                    if (this->pool[k][j][i] == objects::UNLABELLED_INTERFACE) {
                        this->labelInterface(i, j, k);
                    }
                }
            }
        }
    }
}

/** 
 * Checks if there is a fluid connected to this stucture. Used for finding unlabelled interface points
*/
bool Pool3D::isFluidConnected(int i, int j, int k) {
    return this->pool[k][j][i+1] == objects::FLUID_C || this->pool[k][j][i-1] == objects::FLUID_C
            || this->pool[k][j+1][i] == objects::FLUID_C || this->pool[k][j-1][i] == objects::FLUID_C
            || this->pool[k+1][j][i] == objects::FLUID_C || this->pool[k-1][j][i] == objects::FLUID_C;
}

/** 
 * Checks if there is a structure in any of the cardinal directions from this point (only for boundary calls)
*/
bool Pool3D::isStructConnected(int i, int j, int k) {
    return this->pool[k][j][i+1] == objects::STRUCTURE || this->pool[k][j][i-1] == objects::STRUCTURE
            || this->pool[k][j+1][i] == objects::STRUCTURE || this->pool[k][j-1][i] == objects::STRUCTURE
            || this->pool[k+1][j][i] == objects::STRUCTURE || this->pool[k-1][j][i] == objects::STRUCTURE;
}

/**
 * Get normal direction from interface point.
*/
void Pool3D::getNormalDir(objects::FSIObject obj, int nDir[3]) {
    nDir[0] = 0; nDir[1] = 0; nDir[2] = 0;

    if (this->hasStructInDir(obj, objects::EAST)) {
        nDir[0] -= 1;
    } 
    if (this->hasStructInDir(obj, objects::WEST)) {
        nDir[0] += 1;
    }
    if (this->hasStructInDir(obj, objects::NORTH)) {
        nDir[1] -= 1;
    }
    if (this->hasStructInDir(obj, objects::SOUTH)) {
        nDir[1] += 1;
    }
    if (this->hasStructInDir(obj, objects::UP)) {
        nDir[2] -= 1;
    }
    if (this->hasStructInDir(obj, objects::DOWN)) {
        nDir[2] += 1;
    }
}

/**
 * Third order TVD Runge-Kutta method for the various Hamilton-Jacobi equations we need to solve.
 * 
 * Solving the PDE 
 * 
 * S_t = F(S, grad(S), V, S_t0) for some vector field V depending on time, space.
 * 
 * bcType:
 *    0 : Ghost cell method (ENO scheme is unconstrained, bcFun is assumed to set the ghost cells)
 *    1 : Neumann BC without ghost cells: bcFun gives the value of the derivative on outer boundary strip
 *        (in the ghost cell spots). ENO stencil is constrained. 
 *          Note: in this case the bcFun callback is applied in the derivative computation
 *    2 : Dirichlet BC without ghost cells: bcFun gives the values of the output, potentially based on
 *        the values of neighbouring points from the previous time step. ENO stencil is constrained.
*/
void tvdRK3HJ(double dt, double ***arr, Pool3D *pool, int bcType,
                    double (Pool3D::*rhs)(int,int,int,double*,double*,double*,double*,double*,double*),
                    void (Pool3D::*bcFun)(double***)) {
    int i, j, k;
    int start, endx, endy, endz;

    int mo = pool->methodOrd;

    double hx = pool->x[1] - pool->x[0];
    double hy = pool->y[1] - pool->y[0];
    double hz = pool->z[1] - pool->z[0];

    // Start and end of the loop based on the type of BC
    if (bcType == 2) {
        start = 1; endx = pool->nx-2; endy = pool->ny-2; endz = pool->nz-2;
    } else {
        start = 0; endx = pool->nx-1; endy = pool->ny-1; endz = pool->nz-1;
    }

    // Stage 1 - apply BCs to the input array (necissary as boundary conditions
    // will depend on the order of the WENO scheme.
    if (bcType != 1) {
        (pool->*bcFun)(arr);
    }

    // Evaluate the RHS and store the result in the first temp array
    rhsHJENO(hx, hy, hz, arr, pool->phiRk1, bcType, pool, rhs, bcFun);

    // Stage 1 Euler method
    for (k = start; k <= endz; k++) {
        for (j = start; j <= endy; j++) {
            for (i = start; i <= endx; i++) {
                pool->phiRk2[mo+k][mo+j][mo+i]
                    = arr[mo+k][mo+j][mo+i] + dt*pool->phiRk1[mo+k][mo+j][mo+i];
            }
        }
    }

    // Stage 2
    if (bcType != 1) {
        (pool->*bcFun)(pool->phiRk2);
    }

    rhsHJENO(hx, hy, hz, pool->phiRk2, pool->phiRk1, bcType, pool, rhs, bcFun);

    for (k = start; k <= endz; k++) {
        for (j = start; j <= endy; j++) {
            for (i = start; i <= endx; i++) {
                pool->phiRk2[mo+k][mo+j][mo+i] = 0.75*arr[mo+k][mo+j][mo+i]
                    + 0.25*pool->phiRk2[mo+k][mo+j][mo+i]
                    + 0.25*dt*pool->phiRk1[mo+k][mo+j][mo+i];
            }
        }
    }

    // Stage 3 (solution attained)
    if (bcType != 1) {
        (pool->*bcFun)(pool->phiRk2);
    }

    rhsHJENO(hx, hy, hz, pool->phiRk2, pool->phiRk1, bcType, pool, rhs, bcFun);

    for (k = start; k <= endz; k++) {
        for (j = start; j <= endy; j++) {
            for (i = start; i <= endx; i++) {
                arr[mo+k][mo+j][mo+i] = (1.0/3.0)*arr[mo+k][mo+j][mo+i]
                    + (2.0/3.0)*pool->phiRk2[mo+k][mo+j][mo+i]
                    + (2.0/3.0)*dt*pool->phiRk1[mo+k][mo+j][mo+i];
            }
        }
    }

    // With the Dirichlet boundary condition, enforce the condition one more time, to manage the case
    // where it is different on the new time level.
    if (bcType == 2) {
        (pool->*bcFun)(arr);
    }
}

/**
 * RHS function for discretization fo the HJ equation governing movement of the level set
 * function
*/
void rhsHJENO(double hx, double hy, double hz, double ***in, double ***out, int bcType, Pool3D *pool,
            double (Pool3D::*rhs)(int,int,int,double*,double*,double*,double*,double*,double*),
            void (Pool3D::*bcFun)(double***)) {
    int i, j, k, l;
    int start, endx, endy, endz;

    int mo = pool->methodOrd;

    int nvals = 2*mo+1;
    double mxVals[nvals];
    double myVals[nvals];
    double mzVals[nvals];
    double xVals[nvals];
    double yVals[nvals];
    double zVals[nvals];

    // Based on the boundary condition, choose how the RHS is updated on the boundary
    if (bcType == 0) {
        start = 0; endx = pool->nx-1; endy = pool->ny-1; endz = pool->nz-1;
    } else {
        start = 1; endx = pool->nx-2; endy = pool->ny-2; endz = pool->nz-2;
    }

    for (i = 0; i < nvals; i++) {
        mxVals[i] = ((double)i)*hx;
        myVals[i] = ((double)i)*hy;
        mzVals[i] = ((double)i)*hz;
    }

    // Compute the RHS array
    for (k = start; k <= endz; k++) {
        for (j = start; j <= endy; j++) {
            for (i = start; i <= endx; i++) {
                for (l = 0; l < nvals; l++) {
                    yVals[l] = in[mo+k][j+l][mo+i];
                    xVals[l] = in[mo+k][mo+j][i+l];
                    zVals[l] = in[k+l][mo+j][mo+i];
                }
                out[mo+k][mo+j][mo+i] = (pool->*rhs)(i, j, k, mxVals, xVals,
                                            myVals, yVals, mzVals, zVals);
            }
        }
    }

    // If using Neumann BCs, apply the condition on the spatial discretization.
    if (bcType == 1) {
        (pool->*bcFun)(out);
    }
}

/**
 * RHS function for the level set function using the third order ENO scheme
*/
double Pool3D::levelSetRHS_ENO3(int i, int j, int k, double *mxVals, double *xVals,
                                double *myVals, double *yVals, double *mzVals,
                                double *zVals) {
    bool stencil[7] = {true, true, true, true, true, true, true};
    int mo = this->methodOrd;

    int upwindU = (poolU[k+mo][j+mo][i+mo] > 0) ? -1 : 1;
    int upwindV = (poolV[k+mo][j+mo][i+mo] > 0) ? -1 : 1;
    int upwindW = (poolW[k+mo][j+mo][i+mo] > 0) ? -1 : 1;

    double pntX = mxVals[3];
    double pntY = myVals[3];
    double pntZ = mzVals[3];

    return -(poolU[k+mo][j+mo][i+mo]*discs::thirdOrdENO(pntX, mxVals, xVals, upwindU, stencil, 3)
                +   poolV[k+mo][j+mo][i+mo]*discs::thirdOrdENO(pntY, myVals, yVals, upwindV, stencil, 3)
                +   poolW[k+mo][j+mo][i+mo]*discs::thirdOrdENO(pntZ, mzVals, zVals, upwindW, stencil, 3));
}

/**
 * Checks if a point (i, j) is one pixel away from the interface using the methods from Russo
 * and Smereka 2000.
 * 
 * NOTE: assumes that phi contains the original interface
*/
bool Pool3D::oneGridFromInterface(int i, int j, int k) {
    int ox = methodOrd+i;
    int oy = methodOrd+j;
    int oz = methodOrd+k;
    double cVal = phi[oz][oy][ox];
    return phi[oz][oy][ox-1]*cVal < 0 || phi[oz][oy][ox+1]*cVal < 0
        || phi[oz][oy-1][ox]*cVal < 0 || phi[oz][oy+1][ox]*cVal < 0
        || phi[oz+1][oy][ox]*cVal < 0 || phi[oz-1][oy][ox]*cVal < 0; 
}

/**
 * Check if point is one away from an interface in the solid direction
 * 
 * NOTE: assumes that phi contains the original interface
*/
bool Pool3D::oneGridFromInterfaceStructure(int i, int j, int k) {
    int ox = methodOrd+i;
    int oy = methodOrd+j;
    int oz = methodOrd+k;
    return phi[oz][oy][ox] < 0 && (isInterface(objAtIndex(i-1, j, k)) || isInterface(objAtIndex(i+1, j, k))
            || isInterface(objAtIndex(i, j+1, k)) || isInterface(objAtIndex(i, j-1, k))
            || isInterface(objAtIndex(i, j, k+1)) || isInterface(objAtIndex(i, j, k-1)));
}


/**
 * First-order approximation for the distance to the interface. Taken from Russo and Smereka 2000.
 * 
 * NOTE: assumes that phi contains the original interface before the beginning of time integration.
*/
double Pool3D::distanceFromInterface(int i, int j, int k) {
    double hx, hy, hz;
    double eps = 1e-16;
    int ox = methodOrd+i;
    int oy = methodOrd+j;
    int oz = methodOrd+k;

    hx = x[1] - x[0];
    hy = y[1] - y[0];
    hz = z[1] - z[0];

    // Compute centered finite difference discretization of the first derivatives
    double phiGrad[3];
    phiGrad[0] = (phi[oz][oy][ox+1] - phi[oz][oy][ox-1])/(2.0*hx);
    phiGrad[1] = (phi[oz][oy+1][ox] - phi[oz][oy-1][ox])/(2.0*hy);
    phiGrad[2] = (phi[oz+1][oy][ox] - phi[oz-1][oy][ox])/(2.0*hz);

    // Take the Euclidian norm of the gradient
    double phi_norm = simutils::eucNorm3D(phiGrad);

    // Compute the distance approximation, with a small epsilon in the denom to avoid
    // dividing by 0.
    return phi[oz][oy][ox]/(phi_norm+eps);
}

void Pool3D::applyReinitializationSlice(int k, double ***phi) {
    int i, j;

    int mo = this->methodOrd;

    // Sides
    for (j = 0; j <= ny-1; j++) {
        phi[mo+k][mo+j][mo] = phi[mo+k][mo+j][mo+1]
            + this->sign(phi[mo+k][mo+j][mo+1])*abs(phi[mo+k][mo+j][mo+1] - phi[mo+k][mo+j][mo+2]); // LHS
        
        phi[mo+k][mo+j][mo+(nx-1)] = phi[mo+k][mo+j][mo+(nx-1)-1]
            + this->sign(phi[mo+k][mo+j][mo+(nx-1)-1])*abs(phi[mo+k][mo+j][mo+(nx-1)-1] - phi[mo+k][mo+j][mo+(nx-1)-2]); // RHS
    }

    for (i = 0; i <= nx-1; i++) {
        phi[mo+k][mo][mo+i] = phi[mo+k][mo+1][mo+i]
            + this->sign(phi[mo+k][mo+1][mo+i])*abs(phi[mo+k][mo+1][mo+i] - phi[mo+k][mo+2][mo+i]); // BOT

        phi[mo+k][mo+(ny-1)][mo+i] = phi[mo+k][mo+(ny-1)-1][mo+i]
            + this->sign(phi[mo+k][mo+(ny-1)-1][mo+i])*abs(phi[mo+k][mo+(ny-1)-1][mo+i] - phi[mo+k][mo+(ny-1)-2][mo+i]); // TOP
    }
}

/**
 * Apply the reinitialization boundary conditions mentioned in the Fatori thesis.
*/
void Pool3D::applyReinitializationBCs(double ***phi) {
    int mo = this->methodOrd;

    for (int j = 0; j <= ny-1; j++) {
        for (int i = 0; i <= nx-1; i++) {
            // DOWN
            phi[mo][mo+j][mo+i] = phi[mo+1][mo+j][mo+i]
                + this->sign(phi[mo+1][mo+j][mo+i])*abs(phi[mo+1][mo+j][mo+i] - phi[mo+2][mo+j][mo+i]);
            
            // UP
            phi[mo+(nz-1)][mo+j][mo+i] = phi[mo+(nz-1)-1][mo+j][mo+i]
                + this->sign(phi[mo+(nz-1)-1][mo+j][mo+i])*abs(phi[mo+(nz-1)-1][mo+j][mo+i] - phi[mo+(nz-1)-2][mo+j][mo+i]);
        }
    }

    for (int k = 0; k <= nz-1; k++) {
        this->applyReinitializationSlice(k, phi);
    }
}

double Pool3D::signedDistanceReinitializationRHS_CorrectUpwind(int i, int j, int k, double *mxVals,
                                                      double *xVals, double *myVals, 
                                                      double *yVals, double *mzVals, double *zVals) {
    
    // Compute the gradient vector
    int s = 3;
    int l = 2;
    int r = 4;
    double mo = methodOrd;
    int ox = mo+i;
    int oy = mo+j;
    int oz = mo+k;

    double hx, hy, hz, h;

    hx = x[1] - x[0];
    hy = y[1] - y[0];
    hz = z[1] - z[0];
    h = sqrt(simutils::square(hx)+simutils::square(hy)+simutils::square(hz));

    // Compute the gradient using scheme from Russo and Smerka
    double G = 0;
    double a, b, c, d, e, f;
    
    if (this->oneGridFromInterface(i, j, k)) {
        // Use interface-conscious upwind scheme
        return (-this->sign(phi[oz][oy][ox])*abs(xVals[s]) + this->distanceFromInterface(i, j, k))/h;
    } else {
        // Compute the upwind discretization
        if (phi[oz][oy][ox] > 0) {
            a = simutils::dmax((xVals[s] - xVals[l])/hx, 0.0);
            b = simutils::dmin((xVals[r] - xVals[s])/hx, 0.0);
            c = simutils::dmax((yVals[s] - yVals[l])/hy, 0.0);
            d = simutils::dmin((yVals[r] - yVals[s])/hy, 0.0);
            e = simutils::dmax((zVals[s] - zVals[l])/hz, 0.0);
            f = simutils::dmin((zVals[r] - zVals[s])/hz, 0.0);

            G = 1 - sqrt(simutils::dmax(simutils::square(a), simutils::square(b))
                + simutils::dmax(simutils::square(c), simutils::square(d))
                + simutils::dmax(simutils::square(e), simutils::square(f)));
        } else if (phi[oz][oy][ox] < 0) {
            a = simutils::dmin((xVals[s] - xVals[l])/hx, 0.0);
            b = simutils::dmax((xVals[r] - xVals[s])/hx, 0.0);
            c = simutils::dmin((yVals[s] - yVals[l])/hy, 0.0);
            d = simutils::dmax((yVals[r] - yVals[s])/hy, 0.0);
            e = simutils::dmin((zVals[s] - zVals[l])/hz, 0.0);
            f = simutils::dmax((zVals[r] - zVals[s])/hz, 0.0);

            G = 1 - sqrt(simutils::dmax(simutils::square(a), simutils::square(b))
                + simutils::dmax(simutils::square(c), simutils::square(d))
                + simutils::dmax(simutils::square(e), simutils::square(f)));
        } else {
            G = 0;
        }
        return this->sign(phi[oz][oy][ox])*G;
    }
}


/**
 * Compute the gradient of the level set using upwind finite differencing.
 * 
 * TODO: add higher order gradient approximation
*/
void Pool3D::levelSetGradient(int i, int j, int k, double g[3]) {
    double hx, hy, hz;

    int mo = this->methodOrd;
    
    // Mesh spacing in each dimension
    hx = x[1] - x[0];
    hy = y[1] - y[0];
    hz = z[1] - z[0];

    // Indices into phi
    int xi = mo+i;
    int yi = mo+j;
    int zi = mo+k;

    bool ip1Less = phi[zi][yi][xi+1] < phi[zi][yi][xi-1];

    /* Computing the gradient */

    // Take the x partial
    if (ip1Less) {
        g[0] = (phi[zi][yi][xi+1] - phi[zi][yi][xi])/hx;
    } else {
        g[0] = (phi[zi][yi][xi] - phi[zi][yi][xi-1])/hx;
    }

    // Take the y partial
    ip1Less = phi[zi][yi+1][xi] < phi[zi][yi-1][xi];
    if (ip1Less) {
        g[1] = (phi[zi][yi+1][xi] - phi[zi][yi][xi])/hy;
    } else {
        g[1] = (phi[zi][yi][xi] - phi[zi][yi-1][xi])/hy;
    }

    // Take the z partial
    ip1Less = phi[zi+1][yi][xi] < phi[zi-1][yi][xi];
    if (ip1Less) {
        g[2] = (phi[zi+1][yi][xi] - phi[zi][yi][xi])/hz;
    } else {
        g[2] = (phi[zi][yi][xi] - phi[zi-1][yi][xi])/hz;
    }
} 

/** 
 * Compute the outward normal vector of the level set function at grid point (i, j)
 * (assumed to be an interior point)
*/
void Pool3D::levelSetUnitNormal(int i, int j, int k, double n[3]) {

    // Take the gradient
    this->levelSetGradient(i, j, k, n);

    // Normalize the gradient
    simutils::normalize3D(n);
}

/**
 * Apply the periodic BC for slice of level set function
*/
void Pool3D::applyPeriodicLevelSetSlice(int k) {
    int i, j;
    int mo = methodOrd;

    // Apply a periodic boundary condition using the existing pool
    for (j = 0; j < ny; j++) {
        for (i = 0; i < mo; i++) {
            // LEFT
            phi[k+mo][j+mo][i] = phi[k+mo][j+mo][nx+i];

            // RIGHT
            phi[k+mo][j+mo][mo+nx+i] = phi[k+mo][j+mo][mo+i];
        }
    }

    for (i = 0; i < this->nx; i++) {
        for (j = 0; j < mo; j++) {
            // BOT
            phi[k+mo][j][i+mo] = phi[k+mo][ny+j][i+mo];

            // TOP
            phi[k+mo][mo+ny+j][i+mo] = phi[k+mo][mo+j][i+mo];
        }
    }
}

/**
 * Apply 3d periodic boundary conditions.
*/
void Pool3D::applyLevelSetPeriodicBCs(double ***phi) {
    int i, j, k;
    int mo = this->methodOrd;

    // DOWN-MOST
    for (j = 0; j < ny; j++) {
        for (i = 0; i < nx; i++) {
            for (k = 0; k < mo; k++) {
                phi[k][mo+j][mo+i] = phi[nz+k][mo+j][mo+i];
            }
        }
    }
    this->applyPeriodicLevelSetSlice(0);

    // Middle
    for (int k = 1; k < nz-1; k++) {
        this->applyPeriodicLevelSetSlice(k);
    }

    // UP-MOST
    for (j = 0; j < ny; j++) {
        for (i = 0; i < nx; i++) {
            for (k = 0; k < mo; k++) {
                phi[mo+nz+k][mo+j][mo+i] = phi[mo+k][mo+j][mo+i];
            }
        }
    }
    this->applyPeriodicLevelSetSlice(nz-1);
}


/**
 * The numerically smeared sign function for the reinitialization equation.
 * 
 * Takes in a grid spacing h and a scalar value of a signed distance function phi.
*/
double Pool3D::phiSmearSign(double h, double phi, double phi_x, double phi_y,
                            double phi_z) {
    return phi/sqrt(simutils::square(phi)
        + simutils::square(h)*(simutils::square(phi_x)
                                    + simutils::square(phi_y)
                                    + simutils::square(phi_z)));
}

/**
 * The numerically smeared sign function for the reinitialization equation.
 * 
 * Takes in a grid spacing h and a scalar value of a signed distance function phi.
*/
double Pool3D::smearSign(double h, double phi) {
    return phi/sqrt(simutils::square(phi) + simutils::square(h));
}

/**
 * Sign function.
*/
double Pool3D::sign(double phi) {
    return (phi >= 0.0) ? 1.0 : -1.0;
}

/* PUBLIC METHODS */

/**
 * Note that nx, ny, from the perspective of the pool, is the number of subintervals.
*/
Pool3D::Pool3D(Boundary3D &boundary, vector<SolidObject3D> &structures,
               SimParams3D &params) {
    // Create the pool based on the given implementation
    this->create3DPool(boundary, structures, params);
}
/**
 * 
 * Use the two-way velocity extrapolation to extrapolate the velocity into the level set to
 * set constant values along the normal of the level set. Note: two-way formula is more complicated
 * but has the advantage 
 * 
 * Note: this only uses the interior points, i.e., i-1, i+1, j-1, j+1 will never be ghost cells. 
*/
double Pool3D::velocityExtrapolationRHS_ENO3(int i, int j, int k, double *mxVals,
                                             double *xVals, double *myVals, 
                                             double *yVals, double *mzVals,
                                             double *zVals) {

    // Backup plan
    if (objAtIndex(i, j, k) != objects::FLUID_C) {
        return 0.0;
    }

    double sign;
    bool stencil[7] = {true, true, true, true, true, true, true};
    
    int mo = this->methodOrd;

    double hx, hy, hz;
    hx = x[1] - x[0];
    hy = y[1] - y[0];
    hz = z[1] - z[0];

    // Compute the sign function
    sign = this->sign(phi[mo+k][mo+j][mo+i]);

    // Compute the normal vector for the level set at the grid point (i, j)
    double phiGrad[3] = {0.0, 0.0, 0.0};

    // x deriv
    if (i == 0) {
        phiGrad[0] = (-(3.0/2.0)*phi[mo+k][mo+j][mo+i] + 2.0*phi[mo+k][mo+j][mo+i+1] - 0.5*phi[mo+k][mo+j][mo+i+2])/hx;
    } else if (i == nx-1) {
        phiGrad[0] = ((3.0/2.0)*phi[mo+k][mo+j][mo+i] - 2.0*phi[mo+k][mo+j][mo+i-1] + 0.5*phi[mo+k][mo+j][mo+i-2])/hx;
    } else {
        phiGrad[0] = (phi[mo+k][mo+j][mo+i+1] - phi[mo+k][mo+j][mo+i-1])/(2.0*hx);
    }

    
    // y deriv
    if (j == 0) {
        phiGrad[1] = (-(3.0/2.0)*phi[mo+k][mo+j][mo+i] + 2.0*phi[mo+k][mo+j+1][mo+i] - 0.5*phi[mo+k][mo+j+2][mo+i])/hy;
    } else if (j == ny-1) {
        phiGrad[1] = ((3.0/2.0)*phi[mo+k][mo+j][mo+i] - 2.0*phi[mo+k][mo+j-1][mo+i] + 0.5*phi[mo+k][mo+j-2][mo+i])/hy;
    } else {
        phiGrad[1] = (phi[mo+k][mo+j+1][mo+i] - phi[mo+k][mo+j-1][mo+i])/(2.0*hy);
    }

    // z deriv
    if (k == 0) {
        phiGrad[2] = (-(3.0/2.0)*phi[mo+k][mo+j][mo+i] + 2.0*phi[mo+k+1][mo+j][mo+i] - 0.5*phi[mo+k+2][mo+j][mo+i])/hz;
    } else if (k == nz-1) {
        phiGrad[2] = ((3.0/2.0)*phi[mo+k][mo+j][mo+i] - 2.0*phi[mo+k-1][mo+j][mo+i] + 0.5*phi[mo+k-2][mo+j][mo+i])/hz;
    } else {
        phiGrad[2] = (phi[mo+k+1][mo+j][mo+i] - phi[mo+k-1][mo+j][mo+i])/(2.0*hz);
    }

    // Normalize the gradient
    simutils::normalize3D(phiGrad);

    int upwindX = (sign*phiGrad[0] > 0)? -1 : 1;
    int upwindY = (sign*phiGrad[1] > 0)? -1 : 1;
    int upwindZ = (sign*phiGrad[2] > 0)? -1 : 1;

    double pntX = mxVals[3];
    double pntY = myVals[3];
    double pntZ = mzVals[3];

    // Calculate the RHS of the HJ equation using the third order ENO finite difference scheme.
    return -( phiGrad[0]*discs::thirdOrdENO(pntX, mxVals,xVals,upwindX,stencil, 3)
                    + phiGrad[1]*discs::thirdOrdENO(pntY, myVals,yVals,upwindY,stencil, 3)
                    + phiGrad[2]*discs::thirdOrdENO(pntZ, mzVals,zVals,upwindZ,stencil, 3) );
}

/**
 * Calculate the contribution of the current interface cell to the net force
 * for the object it is associated with
*/
void Pool3D::calculateLocalFNet(int i, int j, int k, objects::FSIObject obj, int ng,
                                double ***u, double ***v, double ***w, double ***p) {
    double dA = 0.0;
    double n[3];
    double uGrad[3];
    double vGrad[3];
    double wGrad[3];
    int is, js, is_mv, js_mv, ks, ks_mv;
    double hx = x[1] - x[0];
    double hy = y[1] - y[0];
    double hz = z[1] - z[0];

    double pres;

    // Re-initialize the normal vector to 0
    n[0] = 0.0; n[1] = 0.0; n[2] = 0.0;

    // Set default stencil
    is = i; is_mv = i+1;
    js = j; js_mv = j+1;
    ks = k; ks_mv = k+1;

    if (this->hasStructInDir(obj, objects::EAST)) {
        is = i-1;
        is_mv = i-2;

        n[0] -= 1.0;
    } else if (this->hasStructInDir(obj, objects::WEST)) {
        n[0] += 1.0;
    } 
    
    if (this->hasStructInDir(obj, objects::NORTH))  {
        js = j-1;
        js_mv = j-2;

        n[1] -= 1.0;
    } else if (this->hasStructInDir(obj, objects::SOUTH)) {
        n[1] += 1.0;
    }

    if (this->hasStructInDir(obj, objects::UP)) {
        ks = k-1;
        ks_mv = k-2;

        n[2] -= 1.0;
    } else if (this->hasStructInDir(obj, objects::DOWN)) {
        n[2] += 1.0;
    }

    if (!simutils::eps_equal(n[0], 0.0, 1e-16)) {
        dA = hz*hy;
    } else if (!simutils::eps_equal(n[1], 0.0, 1e-16)) {
        dA = hx*hz;
    } else {
        dA = hx*hy;
    }

    // Compute the information for the hydrodynamic stress tensor.
    // By taking the value the unit normal points to.
    pres = p[ng+k+(int)n[2]][ng+j+(int)n[1]][ng+i+(int)n[0]];

    // Use either backwards or forwards difference formula to compute the derivatives
    // in each direction

    // x derivatives
    if (is < is_mv) {
        // Forward difference
        uGrad[0] = (u[ng+k][ng+j][ng+is_mv] - u[ng+k][ng+j][ng+is])/hx;
        vGrad[0] = (v[ng+k][ng+js][ng+is_mv] - v[ng+k][ng+js][ng+is])/hx;
        wGrad[0] = (w[ng+ks][ng+j][ng+is_mv] - w[ng+ks][ng+j][ng+is])/hx;
    } else {
        // Backward difference
        uGrad[0] = (u[ng+k][ng+j][ng+is] - u[ng+k][ng+j][ng+is_mv])/hx;
        vGrad[0] = (v[ng+k][ng+js][ng+is] - v[ng+k][ng+js][ng+is_mv])/hx;
        wGrad[0] = (w[ng+ks][ng+j][ng+is] - w[ng+ks][ng+j][ng+is_mv])/hx;
    }

    // y derivatives
    if (js < js_mv) {
        uGrad[1] = (u[ng+k][ng+js_mv][ng+is] - u[ng+k][ng+js][ng+is])/hy;
        vGrad[1] = (v[ng+k][ng+js_mv][ng+i] - v[ng+k][ng+js][ng+i])/hy;
        wGrad[1] = (w[ng+ks][ng+js_mv][ng+i] - w[ng+ks][ng+js][ng+i])/hy;
    } else {
        uGrad[1] = (u[ng+k][ng+js][ng+is] - u[ng+k][ng+js_mv][ng+is])/hy;
        vGrad[1] = (v[ng+k][ng+js][ng+i] - v[ng+k][ng+js_mv][ng+i])/hy;
        wGrad[1] = (w[ng+ks][ng+js][ng+i] - w[ng+ks][ng+js_mv][ng+i])/hy;
    }

    // z derivatives
    if (ks < ks_mv) {
        uGrad[2] = (u[ng+ks_mv][ng+j][ng+is] - u[ng+ks][ng+j][ng+is])/hz;
        vGrad[2] = (v[ng+ks_mv][ng+js][ng+i] - v[ng+ks][ng+js][ng+i])/hz;
        wGrad[2] = (w[ng+ks_mv][ng+j][ng+i] - w[ng+ks][ng+j][ng+i])/hz;
    } else {
        uGrad[2] = (u[ng+ks][ng+j][ng+is] - u[ng+ks_mv][ng+j][ng+is])/hz;
        vGrad[2] = (v[ng+ks][ng+js][ng+i] - v[ng+ks_mv][ng+js][ng+i])/hz;
        wGrad[2] = (w[ng+ks][ng+j][ng+i] - w[ng+ks_mv][ng+j][ng+i])/hz;
    }

    // Normalize the unit vector
    simutils::normalize3D(n);

    // Add to the net force
    tracers[domainMembership(i, j, k)].fNetX
        += (-pres*n[0] + this->mu*(2.0*uGrad[0]*n[0]+(uGrad[1]+vGrad[0])*n[1]+(uGrad[2]+wGrad[0])*n[2]))*dA;

    tracers[domainMembership(i, j, k)].fNetY
        += (-pres*n[1] + this->mu*((vGrad[0]+uGrad[1])*n[0] + 2.0*vGrad[1]*n[1] + (vGrad[2]+wGrad[1])*n[2]))*dA;

    tracers[domainMembership(i, j, k)].fNetZ
        += (-pres*n[2] + this->mu*((wGrad[0]+uGrad[2])*n[0] + (wGrad[1]+vGrad[2])*n[1] + 2.0*wGrad[2]*n[2]))*dA;
}

/** 
 * Method which computes the values of the stresses at each of the boundary points.
 * Note the relationship
 * 
 * Note: i, j has to be on the structure domain.
*/
void Pool3D::computeBoundaryStress(int i, int j, int k, objects::FSIObject obj, int ng,
                                double ***u, double ***v, double ***w, double ***p) {
    double n[3] = {0.0, 0.0, 0.0};
    double uGrad[3];
    double vGrad[3];
    double wGrad[3];
    int is, js, is_mv, js_mv, ks, ks_mv;
    double hx = x[1] - x[0];
    double hy = y[1] - y[0];
    double hz = z[1] - z[0];
    int mo = this->methodOrd;

    double pres;

    // Set default stencil
    is = i; is_mv = i+1;
    js = j; js_mv = j+1;
    ks = k; ks_mv = k+1;

    if (this->hasStructInDir(obj, objects::EAST)) {
        is = i-1;
        is_mv = i-2;

        n[0] -= 1.0;
    } else if (this->hasStructInDir(obj, objects::WEST)) {
        n[0] += 1.0;
    } 
    
    if (this->hasStructInDir(obj, objects::NORTH))  {
        js = j-1;
        js_mv = j-2;

        n[1] -= 1.0;
    } else if (this->hasStructInDir(obj, objects::SOUTH)) {
        n[1] += 1.0;
    }

    if (this->hasStructInDir(obj, objects::UP)) {
        ks = k-1;
        ks_mv = k-2;

        n[2] -= 1.0;
    } else if (this->hasStructInDir(obj, objects::DOWN)) {
        n[2] += 1.0;
    }

    assert(abs(n[0]) <= 1 && abs(n[1]) <= 1 && abs(n[2]) <= 2);

    // Compute the information for the hydrodynamic stress tensor.
    // By taking the value the unit normal points to.
    pres = p[ng+k+(int)n[2]][ng+j+(int)n[1]][ng+i+(int)n[0]];

    // Use either backwards or forwards difference formula to compute the derivatives
    // in each direction

    if (is < is_mv) {
        // Forward difference
        uGrad[0] = (u[k][j][is_mv]  - u[k][j][is])/hx;
        vGrad[0] = (v[k][js][is_mv] - v[k][js][is])/hx;
        wGrad[0] = (w[ks][j][is_mv] - w[ks][j][is])/hx;
    } else {
        // Backward difference
        uGrad[0] = (u[k][j][is]  - u[k][j][is_mv])/hx;
        vGrad[0] = (v[k][js][is] - v[k][js][is_mv])/hx;
        wGrad[0] = (w[ks][j][is] - w[ks][j][is_mv])/hx;
    }

    // y derivatives
    if (js < js_mv) {
        uGrad[1] = (u[k][js_mv][is] - u[k][js][is])/hy;
        vGrad[1] = (v[k][js_mv][i]  - v[k][js][i])/hy;
        wGrad[1] = (w[ks][js_mv][i] - w[ks][js][i])/hy;
    } else {
        uGrad[1] = (u[k][js][is] - u[k][js_mv][is])/hy;
        vGrad[1] = (v[k][js][i]  - v[k][js_mv][i])/hy;
        wGrad[1] = (w[ks][js][i] - w[ks][js_mv][i])/hy;
    }

    // z derivatives
    if (ks < ks_mv) {
        uGrad[2] = (u[ks_mv][j][is] - u[ks][j][is])/hz;
        vGrad[2] = (v[ks_mv][js][i] - v[ks][js][i])/hz;
        wGrad[2] = (w[ks_mv][j][i]  - w[ks][j][i])/hz;
    } else {
        uGrad[2] = (u[ks][j][is] - u[ks_mv][j][is])/hz;
        vGrad[2] = (v[ks][js][i] - v[ks_mv][js][i])/hz;
        wGrad[2] = (w[ks][j][i]  - w[ks_mv][j][i])/hz;
    }

    simutils::normalize3D(n);

    // Compute the stresses (assigning them to the pool velocities temporarily for memory saving)
    poolU[mo+k][mo+j][mo+i]
        = (-pres*n[0] + this->mu*(2.0*uGrad[0]*n[0]+(uGrad[1]+vGrad[0])*n[1]+(uGrad[2]+wGrad[0])*n[2]));
    poolV[mo+k][mo+j][mo+i]
        = (-pres*n[1] + this->mu*((vGrad[0]+uGrad[1])*n[0] + 2.0*vGrad[1]*n[1] + (vGrad[2]+wGrad[1])*n[2]));
    poolW[mo+k][mo+j][mo+i]
        = (-pres*n[2] + this->mu*((wGrad[0]+uGrad[2])*n[0] + (wGrad[1]+vGrad[2])*n[1] + 2.0*wGrad[2]*n[2]));
}

/**
 * Set up the Pool velocity field for the case where we have multiple immersed structures.
*/
void Pool3D::multiStructureVelocitySet(double ***u, double ***v, double ***w, int ng) {
    int i, j, k;
    int mo = this->methodOrd;
    objects::FSIObject obj;
    

    for (k = 0; k < this->nz; k++) {
        for (j = 0; j < this->ny; j++) {
            for (i = 0; i < this->nx; i++) {
                obj = objAtIndex(i, j, k);

                if (obj != objects::FLUID_C) {
                    this->poolU[k+mo][j+mo][i+mo] = tracers[domainMembership(i, j, k)].u;
                    this->poolV[k+mo][j+mo][i+mo] = tracers[domainMembership(i, j, k)].v;
                    this->poolW[k+mo][j+mo][i+mo] = tracers[domainMembership(i, j, k)].w;
                } else {
                    this->poolU[k+mo][j+mo][i+mo] = 0.0;
                    this->poolV[k+mo][j+mo][i+mo] = 0.0;
                    this->poolW[k+mo][j+mo][i+mo] = 0.0;
                }
            }
        }
    }

    this->fastMarch(true, 1);
}

/**
 * Function to calculate the velocities of the objects in the pool.
 * Note: has only been implemented in the trivial case of one structure. Must
 *       find a velocity blending function for the more involved cases.
 * 
 * Using the method from Justin's brain paper to calculate force (almost).
 * Assuming that the tracking particals are within the respective level set functions.
 * 
 * NOTE: the use of the midpoint scheme gauranteed second order accuracy for the quadrature.
 * 
 * Note: poolU and poolV are used as temp arrays in this computation
 * 
 * TODO: consider adaptive stencils for the stress tensor. Be careful around the boundaries
*/
void Pool3D::updatePoolVelocities(double dt, double ***u, double ***v, double ***w,
                                    double ***p, int ng) {
    int i, j, k;
    int mo = methodOrd;
    objects::FSIObject obj;

    if (this->nStructs >= 1) {
        // Re-initialize the tracer partical forces to 0
        for (i = 0; i < nStructs; i++) {
            tracers[i].fNetX = 0.0;
            tracers[i].fNetY = 0.0;
            tracers[i].fNetZ = 0.0;
        }

        for (k = 1; k < nz-1; k++) {
            for (j = 1; j < ny-1; j++) {
                for (i = 1; i < nx-1; i++) {
                    // Get the object at the current cell
                    obj = this->objAtIndex(i, j, k);

                    // For each interface point, calculate its contribution to the net force acting on the object.
                    if (isInterface(obj)) {
                        if (!this->isDeformable) {
                            this->calculateLocalFNet(i, j, k, obj, ng, u, v, w, p);
                        } else {
                            this->computeBoundaryStress(i, j, k, obj, ng, u, v, w, p);
                        }
                    } else {
                        this->poolU[mo+k][mo+j][mo+i] = 0.0;
                        this->poolV[mo+k][mo+j][mo+i] = 0.0;
                        this->poolW[mo+k][mo+j][mo+i] = 0.0;
                    }
                }
            }
        }

        // Use the net force being applied on the object to update the tracer velocities using Newton's
        // second law

        // Update the velocities of each of the structures
        if (!this->isDeformable) {
            for (i = 0; i < nStructs; i++) {
                tracers[i].u = tracers[i].u + ((tracers[i].fNetX)/(tracers[i].mass))*dt;
                tracers[i].v = tracers[i].v + ((tracers[i].fNetY)/(tracers[i].mass))*dt;
                tracers[i].w = tracers[i].w + ((tracers[i].fNetZ)/(tracers[i].mass))*dt;

                // Move the mass spring system nodes along with the tracer
                for (int k = 0; k < solids->at(i).pntList->size(); k++) {
                    (*solids->at(i).q)[3*k]   += tracers[i].u*dt;
                    (*solids->at(i).q)[3*k+1] += tracers[i].v*dt;
                    (*solids->at(i).q)[3*k+2] += tracers[i].w*dt;

                    (*solids->at(i).qt)[3*k]   = tracers[i].u;
                    (*solids->at(i).qt)[3*k+1] = tracers[i].v;
                    (*solids->at(i).qt)[3*k+2] = tracers[i].w;

                    solids->at(i).pntList->at(k).x = (*solids->at(i).q)[3*k];
                    solids->at(i).pntList->at(k).y = (*solids->at(i).q)[3*k+1];
                    solids->at(i).pntList->at(k).z = (*solids->at(i).q)[3*k+2];

                    solids->at(i).pntList->at(k).u = (*solids->at(i).qt)[3*k];
                    solids->at(i).pntList->at(k).v = (*solids->at(i).qt)[3*k+1];
                    solids->at(i).pntList->at(k).w = (*solids->at(i).qt)[3*k+2];

                }
            }

            // The the existing pool velocity arrays to 0. The velocity values
            // are already encoded into the tracer particals.
            if (nStructs == 1) {
                simutils::set_constant(nz+2*mo, ny+2*mo, nx+2*mo, tracers[0].u, this->poolU);
                simutils::set_constant(nz+2*mo, ny+2*mo, nx+2*mo, tracers[0].v, this->poolV);
                simutils::set_constant(nz+2*mo, ny+2*mo, nx+2*mo, tracers[0].w, this->poolW);
            } else {
                this->multiStructureVelocitySet(u, v, w, ng);
            }


        } else {
            // For each of the structures, update the deformable solids, compute the net force,
            // and update the positions of the tracer particals.
            if (repulseMode != 0) {
                // TODO: maybe think about a more efficient 
                // cout << "Doing first fast march" << endl;
                this->fastMarch(false, 2);
                // cout << "Finished the fast march" << endl;
                // cout << "Doing detect collisions" << endl;
                this->detectCollisions();
                // cout << "Finished Doing detect collisions" << endl;
                // cout << "setting up domain array" << endl;
                // this->setUpDomainArray();
                // cout << "FINISHED setting up domain array" << endl;
            }

            // cout << "stress" << endl;
            double ***stress[3] = {poolU, poolV, poolW};
            // cout << "done stress" << endl;

            // TODO: don't actually need the net force in principal
            double fNet[3] = {0.0, 0.0, 0.0};

            // cout << "updating tracers" << endl;
            for (int i = 0; i < nStructs; i++) {
                solids->at(i).updateSolidVels(dt, *this, stress, fNet, mo, false);

                // Update the tracer info (not particularly necissary perhaps but will be nice for testing)
                tracers[i].fNetX = fNet[0];
                tracers[i].fNetY = fNet[1];
                tracers[i].fNetZ = fNet[2];

                tracers[i].u = tracers[i].u + ((tracers[i].fNetX)/(tracers[i].mass))*dt;
                tracers[i].v = tracers[i].v + ((tracers[i].fNetY)/(tracers[i].mass))*dt;
                tracers[i].w = tracers[i].w + ((tracers[i].fNetZ)/(tracers[i].mass))*dt;
            }
            // cout << "done updating tracers" << endl;

            // Compute the velocities on the interface points
            double inPnt[3];
            double vel[3];
            // cout << "computing face velocities" << endl;
            for (int k = 0; k < nz; k++) {
                inPnt[2] = simutils::midpoint(z[k], z[k+1]);
                for (int j = 0; j < ny; j++) {
                    inPnt[1] = simutils::midpoint(y[j], y[j+1]);
                    for (int i = 0; i < nx; i++) {
                        inPnt[0] = simutils::midpoint(x[i], x[i+1]);

                        if (this->isInterface(this->objAtIndex(i, j, k))) {
                            try {
                                solids->at(this->domainMembership(i, j, k)).interpFaceVels(inPnt, vel);
                            } catch (const out_of_range& e) {
                                cout << "Caught out of range" << endl;
                                cout << "(i, j, k) = (" << i << ", " << j << ", " << k << " )" << endl;
                                cout << "Domain = " << this->domainMembership(i, j, k) << endl;
                                assert(false);
                            }
                            poolU[mo+k][mo+j][mo+i] = vel[0];
                            poolV[mo+k][mo+j][mo+i] = vel[1];
                            poolW[mo+k][mo+j][mo+i] = vel[2];
                        } else {
                            poolU[mo+k][mo+j][mo+i] = 0.0;
                            poolV[mo+k][mo+j][mo+i] = 0.0;
                            poolW[mo+k][mo+j][mo+i] = 0.0;
                        }
                    }
                }
            }
            // cout << "done computing face velocities" << endl;

            // Now that we have the velocities, we can apply the fast marching algorithm.
            // cout << "fast march the velocities with interp" << endl;
            this->fastMarch(true, 0);
            // cout << "done fast march the velocities" << endl;
        }
    } else {
        cout << "Multiple structures are not yet considered in the pool algorithms" << endl;
        cout << "Whatever the solution is will be far more sophisticated" << endl;
        assert(false);
    }
}

/**
 * Update the pool using the third order ENO scheme. Using periodic
 * BC's for simplicity.
 * 
 * ng := number of ghost cells on the boundary.
*/
void Pool3D::updatePool(double dt, double ***u, double ***v,
                        double ***w, double ***p, int ng, bool reinitialize) {
    // Update the pool using the third order ENO scheme to test level set
    // in the equation. Using periodic boundary condition.

    // If this is called but there are no immersed structures, simply return
    if (this->nStructs < 1) {
        return;
    }

    // Add all of the MSS's
    kdPointCloud->resetCloud();
    for (int i = 0 ; i < nStructs; i++) {
        kdPointCloud->addMSS(solids->at(i));
    }

    // Now rebuild the kdtree for efficient searching.
    kdTree->buildIndex();


    // Update the velocity field
    enumeratePool();
    setUpDomainArray();

    // cout << "updating pool vels" << endl;
    this->updatePoolVelocities(dt, u, v, w, p, ng);
    // cout << "FINISHED updating pool vels" << endl;

    // Reinitilize after the first update. Ensures that we are using a cut-cell approximation;
    if (reinitialize) {
        simutils::copyVals(nz+2*methodOrd, nx+2*methodOrd, ny+2*methodOrd, phiReInit, phi);
    }

    // Update the position of the tracer particals

    // Advance the level set function
    // cout << "updaitng phi" << endl;
    tvdRK3HJ(dt, phi, this, 0, &Pool3D::levelSetRHS_ENO3,
             &Pool3D::applyLevelSetPeriodicBCs);
    // cout << "finished updaitng phi" << endl;

    for (int k = 0; k < nStructs; k++) {
        // cout << "updating tracers" << endl;
        updateTracer(k, dt, 3);
        // cout << "FINSIHED updating tracers" << endl;
    }
    
    // Update the enumeration and domain array
    // cout << "last section" << endl;
    enumeratePool();
    setUpDomainArray();

    for (int i = 0; i < nStructs; i++) {
        solids->at(i).updateSolidLocs(*this, false);
    }

    // If significant drift from the interface detected, correct.
    if (shouldRefitSDF(min(hx, min(hy, hz)))) {
        refitToSolids(ng);
    }
    // cout << "finished last section" << endl;
}

/**
 * Build squeeze field to extrapolate the velocity. Return the CFL constant
*/
double Pool3D::buildSqeezeField() {
    int mo = this->methodOrd;
    double pnt[3];
    double near[3];
    double maxU = -1;
    double maxV = -1;
    double maxW = -1;

    enumeratePool();
    setUpDomainArray();

    // Make the pool velocities on the interface the distance vector to the MSS
    for (int k = 0; k < this->nz; k++) {
        pnt[2] = simutils::midpoint(z[k], z[k+1]);
        for (int j = 0; j < this->ny; j++) {
            pnt[1] = simutils::midpoint(y[j], y[j+1]);
            for (int i = 0; i < this->nx; i++) {
                pnt[0] = simutils::midpoint(x[i], x[i+1]);

                if (isInterface(objAtIndex(i, j, k))) {
                    int structNum = domainMembership(i, j, k);
                    closestBoundaryPnt(structNum, pnt, near);

                    poolU[mo+k][mo+j][mo+i] = near[0] - pnt[0];
                    poolV[mo+k][mo+j][mo+i] = near[1] - pnt[1];
                    poolW[mo+k][mo+j][mo+i] = near[2] - pnt[2];

                    maxU = (abs(poolU[mo+k][mo+j][mo+i]) > maxU) ? abs(poolU[mo+k][mo+j][mo+i]) : maxU;
                    maxV = (abs(poolV[mo+k][mo+j][mo+i]) > maxV) ? abs(poolV[mo+k][mo+j][mo+i]) : maxV;
                    maxW = (abs(poolW[mo+k][mo+j][mo+i]) > maxW) ? abs(poolW[mo+k][mo+j][mo+i]) : maxW;
                }
            }
        }
    }

    // Extrapolate the velocity field
    fastMarch(true, 1);

    double hx = x[1] - x[0];
    double hy = y[1] - y[0];
    double hz = z[1] - z[0];
    double CFL = 1.0 / ((maxU / hx) + (maxV / hy) + (maxW / hz));

    return CFL;
}

/**
 * Checks if the MSS and the signed distance function remain well asigned.
 * If the interface has drifted too far away from the MSS boundarey edges,
 * we apply the squeeze algorithm to bring it back into place and then reinitialize
 * as a signed distance function.
*/
bool Pool3D::shouldRefitSDF(double tol) {
    // assert(false);
    // For each MSS node, Find the value of the level set function at this point.
    // If it is too far from the level set interface, we require squeeze.
    bool squeeze = false;
    massPoint3D pnt;
    for (auto solid = solids->begin(); solid != solids->end(); ++solid) {
        for (auto bId = solid->boundaryNodeIdList->begin(); bId != solid->boundaryNodeIdList->end(); ++bId) {
            pnt = solid->pntList->at(*bId);
            squeeze = abs(interpolatePhi(pnt.x, pnt.y, pnt.z)) > tol;

            if (!squeeze) {
                return squeeze;
            }
        }
    }

    assert(!squeeze);

    return squeeze;
}

/**
 * Use a simple shrinking algorithm to refit the level set function ot the meshes
 * 
 * TODO: may need to reinitialize before using this
*/
void Pool3D::refitToSolids(int ng) {
    double t = 0;

    double CFL = buildSqeezeField();

    while (t < 1.0) {
        // Update the pool using the third order ENO scheme to test level set in the equation.
        // Using periodic boundary condition.
        tvdRK3HJ(CFL, phi, this, 0, &Pool3D::levelSetRHS_ENO3,
                &Pool3D::applyLevelSetPeriodicBCs);
        t += CFL;

        if (t >= 1.0) {
            break;
        }
        double CFL = buildSqeezeField();

        if (t + CFL > 1.0) {
            CFL = 1.0 - t;
        }
    }

    enumeratePool();
    setUpDomainArray();

    // Make this a signed distance function.
    fastMarch(false, 0);

    simutils::copyVals(nx+2*methodOrd, ny+2*methodOrd, nz+2*methodOrd, phiReInit, phi);
}


void Pool3D::printPool() {
    for (int k = 0; k < this->nz; k++) {
        for (int j = 0; j < this->ny; j++) {
            for (int i = 0; i < this->nx; i++) {
                cout << (int)this->pool[k][j][i] << ", ";
            }
            cout << endl;
        }
        cout << endl;
    }
    cout << endl;
}

void Pool3D::outputPool(const char *fname) {
    double x, y, z;

    ofstream outFile;
    outFile.open(fname);

    int tot = nx * ny * nz;
    int onePercent = (int) (tot / (100.0));

    int Pct = 0;
    int off = 0;


    for (int k = 0; k < this->nz; k++) {
        for (int j = 0; j < this->ny; j++) {
            for (int i = 0; i < this->nx; i++) {
                x = simutils::midpoint(this->x[i], this->x[i+1]);
                y = simutils::midpoint(this->y[j], this->y[j+1]);
                z = simutils::midpoint(this->z[k], this->z[k+1]);
                outFile << x << ", ";
                outFile << y << ", ";
                outFile << z << ", ";

                // double inPnt[3] = {x, y, z};
                // double phi =  this->interpolatePhi(x, y, z);

                outFile << this->interpolatePhi(x, y, z) << endl;
                // if (phi > 0) {
                //     int MAXCANDS = 100;
                //     // outFile << closestBoundaryDist(0, inPnt) << endl;
                //     vector<pair<size_t,double> > ret_matches;
                //     std::vector<size_t> ret_index(MAXCANDS);
                //     std::vector<double> out_dist_sqr(MAXCANDS);

                //     // Candidate list
                //     int numFound = kdTree->knnSearch(&inPnt[0],
                //         MAXCANDS, &ret_index[0], &out_dist_sqr[0]);
                    
                //     outFile << sqrt(out_dist_sqr.at(0)) << endl;

                // } else {
                //     outFile << this->interpolatePhi(x, y, z) << endl;
                    
                // }
            }
        }
    }

    outFile.close();
}

void Pool3D::outputPoolVelocity(const char *fname) {
    double x, y, z;
    int mo = this->methodOrd;

    ofstream outFile;
    outFile.open(fname);

    for (int k = 0; k < this->nz; k++) {
        for (int j = 0; j < this->ny; j++) {
            for (int i = 0; i < this->nx; i++) {
                x = simutils::midpoint(this->x[i], this->x[i+1]);
                y = simutils::midpoint(this->y[j], this->y[j+1]);
                z = simutils::midpoint(this->z[k], this->z[k+1]);

                outFile << x << ", ";
                outFile << y << ", ";
                outFile << z << ", ";
                outFile << this->poolU[mo+k][mo+j][mo+i] << ", ";
                outFile << this->poolV[mo+k][mo+j][mo+i] << ", ";
                outFile << this->poolW[mo+k][mo+j][mo+i] << endl;
            }
        }
    }

    outFile.close();
}

void Pool3D::outputSurfaceCentroids(int structNum, const char *fname) {
    (solids->at(structNum)).outputSurfaceCentroids(fname);
}

void Pool3D::outputStructure(int structNum, const char *fname) {
    (solids->at(structNum)).outputEdges(fname);
}

void Pool3D::outputAllStructures(const char *baseName) {
    string baseString(baseName);

    for (int i = 0; i < solids->size(); i++) {
        string temp = baseString + to_string(i);

        (solids->at(i)).outputEdges(temp.c_str());
    }
}

void Pool3D::outputStructureNodes(int structNum, const char *fname) {
    (solids->at(structNum)).outputNodes(fname);
}

void Pool3D::outputAllStructureNodes(const char *baseName) {
    string baseString(baseName);

    for (int i = 0; i < solids->size(); i++) {
        string temp = baseString + to_string(i);

        (solids->at(i)).outputNodes(temp.c_str());
    }
}

void Pool3D::outputStructureVels(int structNum, const char *fname) {
    (solids->at(structNum)).outputNodeVels(fname);
}

void Pool3D::outputAllStructureVels(const char *baseName) {
    string baseString(baseName);

    for (int i = 0; i < solids->size(); i++) {
        string temp = baseString + to_string(i);

        (solids->at(i)).outputNodeVels(temp.c_str());
    }
}

void Pool3D::outputTracers(const char *fname) {
    double x, y, z;

    ofstream outFile;
    outFile.open(fname);

    for (int i = 0; i < nStructs; i++) {
        x = tracers[i].x;
        y = tracers[i].y;
        z = tracers[i].z;

        outFile << x << ", ";
        outFile << y << ", ";
        outFile << z << endl;
    }
}

void Pool3D::outputMedialAxis(const char *fname) {
    double x, y, z;
    double phiVal;

    ofstream outFile;
    outFile.open(fname);

    for (auto tup = medialAxisPnts->begin(); tup != medialAxisPnts->end(); ++tup) {
        x = get<0>(*tup);
        y = get<1>(*tup);
        z = get<2>(*tup);

        phiVal = interpolatePhi(x, y, z);


        outFile << x << ", ";
        outFile << y << ", ";
        outFile << z << ", ";
        outFile << phiVal << endl;
    }

    outFile.close();
}
/**
 * Return the object type at the given index.
*/
objects::FSIObject Pool3D::objAtIndex(int xInd, int yInd, int zInd) {
    return (this->pool)[zInd][yInd][xInd];
}

/** 
 * Check if this object type has a structure cell in the given direction
*/
bool Pool3D::hasStructInDir(objects::FSIObject obj, objects::FSIObject dir) {
    return (obj % dir) != 0;
}

/**
 * Check if the object type is an interface
 * 
*/
bool Pool3D::isInterface(objects::FSIObject obj) {
    return obj != objects::FLUID_C && obj != objects::STRUCTURE;
}

/**
 * Check if the interface has only a single structure as the neighbour
 *
*/
bool Pool3D::isNormalInterface(objects::FSIObject obj) {
    return obj == objects::FLUID_N || obj == objects::FLUID_S
                || obj == objects::FLUID_E || obj == objects::FLUID_W
                || obj == objects::FLUID_U || obj == objects::FLUID_D;

}

bool Pool3D::enoInRangeX(int val) {
    return simutils::in_range(val, 0, nx);
}

bool Pool3D::enoInRangeY(int val) {
    return simutils::in_range(val, 0, ny);
}

bool Pool3D::enoInRangeZ(int val) {
    return simutils::in_range(val, 0, nz);
}

/**
 * Function to indicates whether a given fluid cell is updatable.
 * TODO: Does not handle boundaries correctly so must be applied only
 *       at internal points
 * 
 * A point is only not updatable in U if its neighbour to the right is an interface.
*/
bool Pool3D::isUpdateableU(int i, int j, int k) {
    return pool[k][j][i+1] == objects::FLUID_C;
}

/**
 * Point only not updatable in v if its top neighbour is an interface
*/
bool Pool3D::isUpdateableV(int i, int j, int k) { 
    return pool[k][j+1][i] == objects::FLUID_C;
}

/**
 * Point only not updatable in w if its up neighbour is an interface
*/
bool Pool3D::isUpdateableW(int i, int j, int k) { 
    return pool[k+1][j][i] == objects::FLUID_C;
}

/* Getters and setters */
int Pool3D::getNx() {
    return this->nx;
}
int Pool3D::getNy() {
    return this->ny;
}
int Pool3D::getNz() {
    return this->nz;
}

double Pool3D::getXMeshVal(int i) {
    return this->x[i];
}

double Pool3D::getYMeshVal(int j) {
    return this->y[j];
}

double Pool3D::getZMeshVal(int k) {
    return this->z[k];
}

double *Pool3D::getXMesh() {
    return x;
}

double *Pool3D::getYMesh() {
    return y;
}

double *Pool3D::getZMesh() {
    return z;
}

double Pool3D::getPhiVal(int i, int j, int k) {
    return this->phi[methodOrd+k][methodOrd+j][methodOrd+i];
}

double Pool3D::getObjU(int i, int j, int k) {
    return this->poolU[methodOrd+k][methodOrd+j][methodOrd+i];
}

double Pool3D::getObjV(int i, int j, int k) {
    return this->poolV[methodOrd+k][methodOrd+j][methodOrd+i];
}

double Pool3D::getObjW(int i, int j, int k) {
    return this->poolW[methodOrd+k][methodOrd+j][methodOrd+i];
}


//-----------------------
/* Tracer utilities */
// -------------------

/**
 * Getter for the domain membership of this index.
*/
int Pool3D::domainMembership(int i, int j, int k) {
    return this->domainTracker[k][j][i];
}

/**
 * Sets up the domain matrix, manipulating the instance vars the correct way.
 * 
 * Uses the domain tracers, so these must be properly updated at all times.
*/
void Pool3D::setUpDomainArray() {
    // Set all of the domain member array to "undescovered"
    simutils::set_constant(nz, ny, nx, DOMAIN_UNDESCOVERED_3D, this->domainTracker);

    // For each immersed structure in the Pool, label subset of the domain they
    // take up.
    for (int l = 0; l < this->nStructs; l++) {
        bfsFromTracer(l);
    }

    // Label all undiscovered domain points as belonging to the fluid domain
    for (int k = 0; k < this->nz; k++) {
        for (int j = 0; j < this->ny; j++) {
            for (int i = 0; i < this->nx; i++) {
                if (this->domainTracker[k][j][i] == DOMAIN_UNDESCOVERED_3D) {
                    this->domainTracker[k][j][i] = DOMAIN_FLUID_3D;
                }
            }
        }
    }
}

/**
 * Performs the bfs to find the structure region for a particular structure in the pool.
 * 
 * Note: where no search has occured, the tracker array is assumed to have value -2. When a node
 *       is discovered, this is changed to -1. Note: using this setup, we actually have a method
 *       to check for when two structures have merged.
 * 
 * A Bredth-first search is done about the pool region that the tracer for this object is
 * currently contained in.
 * 
 * The result is a hashmap (unordered_map) with keys (i, j) coordinates in the pool,
 * with values that are tuples ()
*/
void Pool3D::bfsFromTracer(int structNum) {
    // Find the grid square that the tracer is currently on. This is the root node for the BFS.
    int i = simutils::findLimInfMeshPoint(tracers[structNum].x, this->x, this->nx+1);
    int j = simutils::findLimInfMeshPoint(tracers[structNum].y, this->y, this->ny+1);
    int k = simutils::findLimInfMeshPoint(tracers[structNum].z, this->z, this->nz+1);

    int ri, rj, rk, ni, nj, nk;
    objects::FSIObject FL = objects::FLUID_C; // Make the ternary operator more palatable below

    // Create the queues for each of the indices
    queue<int> queueX;
    queue<int> queueY;
    queue<int> queueZ;

    // Add the first node
    queueX.push(i);
    queueY.push(j);
    queueZ.push(k);

    // Mark this node as a part of the structure. Otherwise, set it to "discovered"
    this->domainTracker[k][j][i] = (objAtIndex(i, j, k) != FL) ? structNum : -1;

    // Start the DFS
    while (!queueX.empty()) {
        ri = queueX.front();
        rj = queueY.front();
        rk = queueZ.front();

        queueX.pop();
        queueY.pop();
        queueZ.pop();

        // Discover the neighbouring points
        for (nk = rk-1; nk <= rk+1; nk++) {
            for (nj = rj-1; nj <= rj+1; nj++) {
                for (ni = ri-1; ni <= ri+1; ni++) {
                    if (ni == ri && nj == rj && nk == rk) {
                        continue;
                    }

                    if (this->domainTracker[nk][nj][ni] == DOMAIN_UNDESCOVERED_3D) { // If this node has not been seen
                        if (objAtIndex(ni, nj, nk) != FL) {
                            // This node belongs to the structure
                            this->domainTracker[nk][nj][ni] = structNum;
                            queueX.push(ni);
                            queueY.push(nj);
                            queueZ.push(nk);
                        } else {
                            // This node belongs to the fluid region
                            this->domainTracker[nk][nj][ni] = DOMAIN_FLUID_3D;
                        }
                    }
                }
            }
        }
    }

    assert(queueX.empty() && queueY.empty() && queueZ.empty());
}

// -------------------

Pool3D::~Pool3D() {
    int mo = this->methodOrd;

    delete solids;
    
    delete[] this->x;
    delete[] this->y;
    delete[] this->z;

    if (this->nStructs > 0) {
        delete[] this->tracers;
    }

    for (int k = 0; k < this->nz; k++) {
        for (int j = 0; j < this->ny; j++) {
            delete[] this->pool[k][j];
        }
        delete[] this->pool[k];
    }
    delete[] this->pool;

    simutils::free_double(this->nz+2*mo, this->ny+2*mo, this->nx+2*mo, this->phi);
    simutils::free_double(this->nz+2*mo, this->ny+2*mo, this->nx+2*mo, this->phiRk1);
    simutils::free_double(this->nz+2*mo, this->ny+2*mo, this->nx+2*mo, this->phiRk2);
    simutils::free_double(this->nz+2*mo, this->ny+2*mo, this->nx+2*mo, this->phiReInit);
    simutils::free_double(this->nz+2*mo, this->ny+2*mo, this->nx+2*mo, this->poolU);
    simutils::free_double(this->nz+2*mo, this->ny+2*mo, this->nx+2*mo, this->poolV);
    simutils::free_double(this->nz+2*mo, this->nz+2*mo, this->nx+2*mo, this->poolW);
    simutils::free_int(this->nz, this->ny, this->nx, this->domainTracker);
}


