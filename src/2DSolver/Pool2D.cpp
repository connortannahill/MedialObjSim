#include "Pool2D.h"
#include "Boundary.h"
#include "SolidObject.h"
#include "../Utils/SimUtilities.h"
#include <math.h>
#include <limits>
#include <fstream>
#include "../Utils/Discretizations.h"
#include <assert.h>
#include <utility>
#include <map>
#include <queue>
#include "MassSpring2D.h"

using namespace mass_spring;
using namespace std;

// Infinity 
const double DINFINITY = numeric_limits<double>::infinity();

// Safety factors for various methods
// const double SFAC_REINIT = 0.1;

int DOMAIN_FLUID = -1;
int DOMAIN_UNDESCOVERED = -2;
int DOMAIN_INTERSECTION = -3;


/* PRIVATE METHODS */

const int CLOSE = 0;
const int ACCEPTED = 1;
const int FAR = 2;

bool Pool2D::indInRange(int i, int j) {
    return !(i < 0 || i > nx-1 || j < 0 || j > ny-1);
}

void Pool2D::fastMarchSetNVal(int i, int j, bool nExtrap, int mode) {
    double hx = x[1] - x[0];
    double hy = y[1] - y[0];
    double hxPerm = x[1] - x[0];
    double hyPerm = y[1] - y[0];
    int mo = methodOrd;

    // Ensure this point is in bounds
    if (!(indInRange(i, j))) {
        return;
    }

    // If this point is in far set it to close, else return
    if (fastMarchingState[j][i] == FAR) {
        fastMarchingState[j][i] = CLOSE;
    } else {
        return;
    }

    // Use upwinding in each axis to determine upwinding direction for this point,
    // as well as the point value to be used in solving the Eikonal equation
    int upwindX = 0;
    int upwindY = 0;

    // x
    if (i+1 > nx-1 || i-1 < 0) {
        // First checks: the point must be in bounds. If one of the points is OOB but the
        // other is in far, we can't use this direction.
        upwindX = (i+1 > nx-1) ? -1 : 1;
        if (fastMarchingState[j][i+upwindX] == FAR) {
            upwindX = 0;
        }
    } else if (fastMarchingState[j][i+1] == FAR && fastMarchingState[j][i-1] == FAR) {
        // Secondly, if both of the points are in FAR, we can't use this direction
        upwindX = 0;
    } else if (fastMarchingState[j][i+1] != FAR && fastMarchingState[j][i-1] != FAR) {
        // Next, if we have both directions, find the proper upwind direction
        upwindX = (phiReInit[mo+j][mo+i-1] < phiReInit[mo+j][mo+i+1]) ? -1 : 1;

        // If the upwind direction is an interface, this will always be the upwind direction
        if (oneGridFromInterfaceStructure(i, j)) {
            if (isInterface(objAtIndex(i-1, j))) {
                upwindX = -1;
            } else if (isInterface(objAtIndex(i+1, j))) {
                upwindX = 1;
            }
        }
    } else {
        // Finally, if there is only one direction, and there are no issues with bounds, we choose this one.
        upwindX = (fastMarchingState[j][i+1] != FAR) ? 1 : -1;
    }

    // now y
    if (j+1 > ny-1 || j-1 < 0) {
        upwindY = (j+1 > ny-1) ? -1 : 1;
        if (fastMarchingState[j+upwindY][i] == FAR) {
            upwindY = 0;
        }
    } else if (fastMarchingState[j+1][i] == FAR && fastMarchingState[j-1][i] == FAR) {
        // Secondly, if both of the points are in FAR, we can't use this direction
        upwindY = 0;
    } else if (fastMarchingState[j+1][i] != FAR && fastMarchingState[j-1][i] != FAR) {
        // Next, if we have both directions, find the proper upwind direction
        upwindY = (phiReInit[mo+j-1][mo+i] < phiReInit[mo+j+1][mo+i]) ? -1 : 1;

        if (oneGridFromInterfaceStructure(i, j)) {
            if (isInterface(objAtIndex(i, j-1))) {
                upwindY = -1;
            } else if (isInterface(objAtIndex(i, j+1))) {
                upwindY = 1;
            }
        }
    } else {
        // Finally, if there is only one direction, and there are no issues with bounds, we choose this one.
        upwindY = (fastMarchingState[j+1][i] != FAR) ? 1 : -1;
    }

    // Now we apply the method outlined in Bridson 2015 to approximate the solution to the
    // reinitialization Eikonal equation
    double phiH, phiV;
    double uH, uV;
    double vH, vV;
    double d;
    const double EPS = 1e-16;

    // Case that should be impossible
    assert(!(upwindX == 0 && upwindY == 0));

    int path;

    if (upwindX == 0)  {
        // Compute signed distance, if difference of sign, approximate the distance to
        // the interface
        hx = hxPerm;
        hy = hyPerm;

        if (isInterface(objAtIndex(i, j+upwindY)) && objAtIndex(i, j) == objects::STRUCTURE) {
            double y0 = simutils::midpoint(y[j], y[j+1]);
            double y1 = simutils::midpoint(y[j+upwindY], y[j+upwindY+1]);
            phiV = 0.0;
            d = abs(y0 - ((phi[mo+j+upwindY][mo+i]*y0 - phi[mo+j][mo+i]*y1)
                / (phi[mo+j+upwindY][mo+i] - phi[mo+j][mo+i]+EPS))); 
            hy = d;
        } else {
            phiV = phiReInit[mo+j+upwindY][mo+i];
            d = phiV + hy;
        }

        // d = phiV + hy;

        // Extrapolate velocity
        if (nExtrap) {
            poolU[mo+j][mo+i] = poolU[mo+j+upwindY][mo+i];
            poolV[mo+j][mo+i] = poolV[mo+j+upwindY][mo+i];
            path = 0;
        }
    } else if (upwindY == 0) {
        hx = hxPerm;
        hy = hyPerm;
        if (isInterface(objAtIndex(i+upwindX, j)) && objAtIndex(j, i) == objects::STRUCTURE) {
            double x0 = simutils::midpoint(x[i], x[i+1]);
            double x1 = simutils::midpoint(x[i+upwindX], x[i+upwindX+1]);
            phiH = 0.0;
            d = abs(x0 - ((phi[mo+j][mo+i+upwindX]*x0 - phi[mo+j][mo+i]*x1)
                / (phi[mo+j][mo+i+upwindX] - phi[mo+j][mo+i]+EPS))); 
            hx = d;
        } else {
            phiH = phiReInit[mo+j][mo+i+upwindX];
            d = phiH + hx;
        }

        // Extrapolate velocity
        if (nExtrap) {
            poolU[mo+j][mo+i] = poolU[mo+j][mo+i+upwindX];
            poolV[mo+j][mo+i] = poolV[mo+j][mo+i+upwindX];
        }
        path = 1;
    } else {
        // Full algorithm
        double phi0, phi1, h;
        char less;
        hx = hxPerm;
        hy = hyPerm;

        // phiV = phiReInit[mo+j+upwindY][mo+i];
        // if (isInterface(objAtIndex(i, j+upwindY)) && objAtIndex(i, j) == objects::STRUCTURE) {
        //     phiV = 0.0;
        // }
        if (isInterface(objAtIndex(i, j+upwindY)) && objAtIndex(i, j) == objects::STRUCTURE) {
            double y0 = simutils::midpoint(y[j], y[j+1]);
            double y1 = simutils::midpoint(y[j+upwindY], y[j+upwindY+1]);
            phiV = 0.0;
            d = abs(y0 - ((phi[mo+j+upwindY][mo+i]*y0 - phi[mo+j][mo+i]*y1)
                / (phi[mo+j+upwindY][mo+i] - phi[mo+j][mo+i]+EPS))); 
            hy = d;
        } else {
            phiV = phiReInit[mo+j+upwindY][mo+i];
            d = phiV + hy;
        }

        if (isInterface(objAtIndex(i+upwindX, j)) && objAtIndex(j, i) == objects::STRUCTURE) {
            double x0 = simutils::midpoint(x[i], x[i+1]);
            double x1 = simutils::midpoint(x[i+upwindX], x[i+upwindX+1]);
            phiH = 0.0;
            d = abs(x0 - ((phi[mo+j][mo+i+upwindX]*x0 - phi[mo+j][mo+i]*x1)
                / (phi[mo+j][mo+i+upwindX] - phi[mo+j][mo+i]+EPS))); 
            hx = d;
        } else {
            phiH = phiReInit[mo+j][mo+i+upwindX];
            d = phiH + hx;
        }

        // phiH = phiReInit[mo+j][mo+i+upwindX];

        // if (isInterface(objAtIndex(i+upwindX, j)) && objAtIndex(j, i) == objects::STRUCTURE) {
        //     phiH = 0.0;
        // }
        // d = phiH + hx;

        if (phiV < phiH) {
            h = hy;
            phi0 = phiV;
            phi1 = phiH;
            less = 'y';
        } else {
            h = hx;
            phi0 = phiH;
            phi1 = phiV;
            less = 'x';
        }

        // Attempting 1 point
        d = phi0 + h;

        double hxsq = simutils::square(hx);
        double hysq = simutils::square(hy);

        if (d > phi1) {
            // Try two closest neighbours
            d = 1.0/((hxsq+hysq)) * ((hysq*phiH+hxsq*phiV)
                        + sqrt(simutils::square(hysq*phiH+hxsq*phiV)
                            - (hxsq+hysq)*(hysq*simutils::square(phiH)
                            + hxsq*simutils::square(phiV) - hysq*hxsq)));

            // Extrapolate the velocities
            // u
            if (nExtrap) {
                uH = poolU[j+mo][i+mo+upwindX];
                uV = poolU[j+mo+upwindY][i+mo];

                // v
                vH = poolV[j+mo][i+mo+upwindX];
                vV = poolV[j+mo+upwindY][i+mo];

                // poolU[j+mo][i+mo] = (hysq*uH*(d - phiH) + hxsq*uV*(d - phiV))/(hxsq*(d - phiH) + hysq*(d - phiV));
                // poolV[j+mo][i+mo] = (hysq*vH*(d - phiH) + hxsq*vV*(d - phiV))/(hxsq*(d - phiH) + hysq*(d - phiV));
                poolU[j+mo][i+mo] = (uH*(d - phiH) + uV*(d - phiV))/((d - phiH) + (d - phiV));
                poolV[j+mo][i+mo] = (vH*(d - phiH) + vV*(d - phiV))/((d - phiH) + (d - phiV));
                path = 2;
            }

        } else {
            if (nExtrap) {
                if (less == 'x') {
                    poolU[mo+j][mo+i] = poolU[mo+j][mo+i+upwindX];
                    poolV[mo+j][mo+i] = poolV[mo+j][mo+i+upwindX];
                } else {
                    poolU[mo+j][mo+i] = poolU[mo+j+upwindY][mo+i];
                    poolV[mo+j][mo+i] = poolV[mo+j+upwindY][mo+i];
                }
            }

            path = 3;
        }
    }

    // Assign the distance
    phiReInit[mo+j][mo+i] = d;

    if (mode == 2) {
        if (d > collisionDist) {
           fastMarchingState[j][i] = ACCEPTED;
           return;
        }
    }

    // Add the point to the heap structure
    heap.push(GridVal2D(i, j, phiReInit[mo+j][mo+i]));
}

/**
 * Assign the marching state of the current node.
*/
void Pool2D::assignDomainMemberships(int i, int j, int mode) {
    // cout << "assigning domain" << endl;
    int curMembership = domainMembership(i, j);
    int neighMembership;

    int iList[4] = {i+1, i-1, i,   i};
    int jList[4] = {j,   j,   j+1, j-1};

    double nIntersections = 0;
    bool intersectionFound = false;
    double medX = 0;
    double medY = 0;

    for (int l = 0; l < 4; l++) {
        int ni = iList[l];
        int nj = jList[l];

        if (indInRange(ni, nj)) {
            neighMembership = domainMembership(ni, nj);

            if (neighMembership == DOMAIN_FLUID) {
                domainTracker[nj][ni] = curMembership;
            } else if (neighMembership != curMembership && domainTracker[nj][ni] != DOMAIN_INTERSECTION) {

                // Keep track of the number if intersections about this point. We accumulate
                // them and take the average.
                nIntersections++;

                int offX = (ni - i > 0) ? 1 : 0;
                int offY = (nj - j > 0) ? 1 : 0;

                medX += x[i+offX];
                medY += y[j+offY];

                domainTracker[j+offY][i+offX] = DOMAIN_INTERSECTION;

                intersectionFound = true;
            }
        }
    }

    if (intersectionFound) {
        medialAxisPnts->push_back(make_tuple(medX/nIntersections, medY/nIntersections));
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
 * TODO: can look for massive performance improvements by (somehow) only doing a single loop over all of the
 *       unknowns, rather than the 2nd and 3rd currently required when extrapolating in the object interior.
 * 
 * mode controls what the fast marching algorithm does
 *  mode = 0: the interpolation based reinititialization
 *  mode = 1: use of the previous isocontour
 *  mode = 2: collision plane detection
*/
void Pool2D::fastMarchPool(bool nExtrap, int mode) {
    // Identify all of the points on the boundaries and initialize the heap with their signed distance
    // function values.
    int mo = methodOrd;

    // Clear the medial axis array
    medialAxisPnts->clear();

    // Initialize the velocity fields and the temp signed distance function to infinity (negative
    // infinity in the interior of the domain.)

    assert(heap.empty());

    double phiVal = 0.0;
    objects::FSIObject obj;
    double pnt[2];
    double vTemp[2];

    for (int j = 0; j < this->ny; j++) {
        pnt[1] = simutils::midpoint(y[j], y[j+1]);
        for (int i = 0; i < this->nx; i++) {
            pnt[0] = simutils::midpoint(x[i], x[i+1]);

            phiVal = INFINITY;
            obj = this->objAtIndex(i, j);

            // Set the default values in the fluid region.
            if (!isInterface(obj)) {
                // If fluid cell or structure cell, set garbage values
                if (nExtrap) {
                    this->poolU[mo+j][mo+i] = DINFINITY;
                    this->poolV[mo+j][mo+i] = DINFINITY;
                }
                this->phiReInit[mo+j][mo+i] = DINFINITY;
            } else {
                // If interface cell, set phi to the distance from the nearest interface pnt
                if (mode == 0) {
                    phiVal = sign(phi[mo+j][mo+i])*solids->at(domainMembership(i, j)).closestBoundaryDist(pnt);
                } else {
                    phiVal = phi[mo+j][mo+i];
                }
                this->phiReInit[mo+j][mo+i] = phiVal;

                // If doing velocity extrapolation, set the velocities at the interface
                if (nExtrap) {
                    if (mode == 0) {
                        solids->at(domainMembership(i, j)).interpFaceVels(pnt, vTemp);

                        this->poolU[mo+j][mo+i] = vTemp[0];
                        this->poolV[mo+j][mo+i] = vTemp[1];
                    }
                }
            }

            // Label each of the points
            fastMarchingState[j][i] = (obj == objects::FLUID_C) ? FAR : ACCEPTED;

            // Add the value of the curent level set (assumed SDF) at the interfaces.
            if (this->isInterface(obj)) {
                fastMarchingState[j][i] = CLOSE;
                heap.push(GridVal2D(i, j, phiVal));
            }
        }
    }

    // Now run the FM algorithm in the fluid region
    int i, j;
    double val;
    int count = 0;
    while (!heap.empty()) {
        // Get and pop the top value 
        GridVal2D pnt = heap.top();
        heap.pop();

        // Extract the point information
        i = pnt.getX();
        j = pnt.getY();
        val = pnt.getVal();

        // Assign the doman memberships of the neighbouring points
        if (domainTracker[j][i] != DOMAIN_INTERSECTION) {
            assignDomainMemberships(i, j, mode);
        }

        // Update all of the neighbours according to the fast marching algorithm
        // and add them to the heap structure
        fastMarchSetNVal(i+1, j, nExtrap, mode);
        fastMarchSetNVal(i-1, j, nExtrap, mode);
        fastMarchSetNVal(i,   j+1, nExtrap, mode);
        fastMarchSetNVal(i,   j-1, nExtrap, mode);

        fastMarchingState[j][i] = ACCEPTED;
        count++;
    }

    // In the case where we are only approximating the medial axis, we don't need to perform
    // the structural domain calculation.
    if (mode == 2) {
        return;
    }


    if (mode != 0) {
        for (int j = 0; j < this->ny; j++) {
            for (int i = 0; i < this->nx; i++) {
                obj = objAtIndex(i, j);

                fastMarchingState[j][i] = (isInterface(obj)) ? CLOSE : FAR;
            }
        }
    }

    // Now: run the same algorithm in the structure domain. Note the interface
    //      and fluid point have been moved into the accepted class by this point
    for (int j = 0; j < this->ny; j++) {
        pnt[1] = simutils::midpoint(y[j], y[j+1]);
        for (int i = 0; i < this->nx; i++) {
            pnt[0] = simutils::midpoint(x[i], x[i+1]);

            // We only act on structure grids here
            obj = this->objAtIndex(i, j);

            if (obj == objects::STRUCTURE) {
                phiVal = INFINITY; // Use the negated phi values.
                assert(phiVal >= 0);

                if (this->oneGridFromInterfaceStructure(i, j)) {
                    // cout << "Getting phi val from vec, domainMem = " << domainMembership(i, j) << endl;
                    if (mode == 0) {

                        phiVal = solids->at(domainMembership(i, j)).closestBoundaryDist(pnt);

                        // Normal extrapolation (very approximate)
                        if (nExtrap) {
                            // Use interpolation of face velocities
                            solids->at(domainMembership(i, j)).interpFaceVels(pnt, vTemp);

                            this->poolU[mo+j][mo+i] = vTemp[0];
                            this->poolV[mo+j][mo+i] = vTemp[1];
                        }

                        // Put this in the close field and push it onto the heap
                        heap.push(GridVal2D(i, j, phiVal));
                    } else {
                        fastMarchSetNVal(i, j, nExtrap, mode);
                     }
                } else {
                    // Otherwise, set the velocity to infinity to make errors easy to spot.
                    if (nExtrap) {
                        this->poolU[mo+j][mo+i] = DINFINITY;
                        this->poolV[mo+j][mo+i] = DINFINITY;
                    }
                }
            }

            // Seperate the logic for assigning the marching state for brain reasons
            if (mode == 0) {
                if (obj == objects::STRUCTURE && this->oneGridFromInterfaceStructure(i, j)) {
                    fastMarchingState[j][i] = CLOSE;
                    phiReInit[mo+j][mo+i] = phiVal;
                } else if (obj == objects::STRUCTURE) {
                    fastMarchingState[j][i] = FAR;
                    phiReInit[mo+j][mo+i] = DINFINITY;
                } else {
                    fastMarchingState[j][i] = ACCEPTED;
                }
            }
        }
    }

    if (mode != 0) {
        for (int j = 0; j < this->ny; j++) {
            for (int i = 0; i < this->nx; i++) {
                obj = objAtIndex(i, j);
                if (obj == objects::STRUCTURE && this->oneGridFromInterfaceStructure(i, j)) {
                    fastMarchingState[j][i] = CLOSE;
                } else if (obj == objects::STRUCTURE) {
                    fastMarchingState[j][i] = FAR;
                    phiReInit[mo+j][mo+i] = DINFINITY;
                } else {
                    fastMarchingState[j][i] = ACCEPTED;
                }
            }
        }
    }

    // Now, run fast marching on the interior data.
    int cnt = 0;
    while (!heap.empty()) {
        GridVal2D pnt = heap.top();
        heap.pop();

        // Extract the point information
        i = pnt.getX();
        j = pnt.getY();
        val = pnt.getVal();

        // Update all of the neighbours according to the fast marching algorithm
        // and add them to the heap structure
        fastMarchSetNVal(i+1, j, nExtrap, mode);
        fastMarchSetNVal(i-1, j, nExtrap, mode);
        fastMarchSetNVal(i,   j+1, nExtrap, mode);
        fastMarchSetNVal(i,   j-1, nExtrap, mode);

        fastMarchingState[j][i] = ACCEPTED;
        cnt++;
    }

    // Finally, flip the sign of the phi value in all of structure points.
    for (int j = 0; j < this->ny; j++) {
        for (int i = 0; i < this->nx; i++) {
            if (objAtIndex(i, j) == objects::STRUCTURE) {
                phiReInit[j+mo][i+mo] *= -1;
            }
        }
    }
}

/**
 * Algorithm to detect collisions and setup data structures for possible force calculations
*/
void Pool2D::detectCollisions() {

    // Wipe the collision memory from the MSS's
    for (auto mss = solids->begin(); mss != solids->end(); ++mss) {

        // For each node, clear the set in the node collisions
        for (auto nodeSet = mss->nodeCols->begin(); nodeSet != mss->nodeCols->end(); ++nodeSet) {
            nodeSet->clear();
        }
    }

    // Now, check the medial axis for any possible collisions using thresholding on phi.
    double medX, medY, medPhi;
    int xCell, yCell;
    set<int> allCollisions;
    set<int> collidingIds;
    vector<pair<double, double>> medialAxisCollisionPnts;

    for (auto tup = medialAxisPnts->begin(); tup != medialAxisPnts->end(); ++tup) {
        // Find the cell location of the medial axis point
        medX = get<0>(*tup);
        medY = get<1>(*tup);

        xCell = simutils::findLimInfMeshPoint(medX, this->x, this->nx+1);
        yCell = simutils::findLimInfMeshPoint(medY, this->y, this->ny+1);

        // Find the value of the level set function in this cell
        medPhi = interpolatePhi(medX, medY);//phiReInit[yCell+mo][xCell+mo];

        // If this value is beneath the threshold for possible collision, look at which objects are colliding
        if (medPhi < repulseDist) {
            // Build up a list of the neighbours around this cell to find which unique colliding nodes
            for (int j = yCell-2; j <= yCell+2; j++) {
                for (int i = xCell-2; i <= xCell+2; i++) {
                    if (indInRange(i, j) && domainMembership(i, j) != DOMAIN_INTERSECTION
                        && domainMembership(i, j) != DOMAIN_FLUID) {
                        collidingIds.insert(domainMembership(i, j));
                        allCollisions.insert(domainMembership(i, j));
                    }
                }
            }

            pair<double, double> medLoc(medX, medY);
            medialAxisCollisionPnts.push_back(medLoc);
        }

        collidingIds.clear();
    }

    // If we are using proper repulsive forces, we now build up the KDTree and
    // all of the node collections
    if (repulseMode == 2 && medialAxisCollisionPnts.size() > 0) {
        // cout << "DETECTED A COLLISION" << endl;
        // Place all of the colliding MSS's into the KDTree
        kdPointCloud->resetCloud();
        for (auto colMss = allCollisions.begin(); colMss != allCollisions.end(); ++colMss) {
            cout << "Adding MSS " << *colMss << endl;
            kdPointCloud->addMSS(solids->at(*colMss));
        }

        // Now rebuild the kdtree for efficient searching.
        kdTree->buildIndex();
        vector<pair<size_t,double> > ret_matches;
        nanoflann::SearchParams params;
        double SCAL_FAC = 1.2;
        int nMatches;
        for (auto pair = medialAxisCollisionPnts.begin(); pair != medialAxisCollisionPnts.end(); ++pair) {
            medX = pair->first;
            medY = pair->second;

            double queryPnt[2] = {medX, medY};

            // Do a radius search around this medial axis cell to discover all of the
            // potentially interacting nodes
            nMatches = kdTree->radiusSearch(&queryPnt[0],
                simutils::square(SCAL_FAC*repulseDist), ret_matches, params);

            // cout << "nMatches = " << nMatches << endl;

            // Find the nearest points on each colliding object to the medial axis point
            vector<int> nearestMss;
            vector<int> nearestPointMatch;
            int id = 0;
            for (auto match = ret_matches.begin(); match != ret_matches.end(); ++match) {
                int mssId = kdPointCloud->points->at(match->first)->structNum;
                // cout << "mssId of match point: " << mssId  << endl;

                // If this is a new Id, since the array is ordered according to proximity to m, we add it
                if (std::find(nearestMss.begin(), nearestMss.end(), mssId) == nearestMss.end()) {
                    nearestMss.push_back(mssId);
                    nearestPointMatch.push_back(id);
                }
                id++;
            }
            
            if (nearestMss.size() < 2) {
                cout << "ERROR IN FINDING MSS NODES AT COLLISION POINT" << endl;
                assert(nearestMss.size() < 2);
            }
            
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
    // assert(false);
}

void Pool2D::create2DPool(Boundary &boundary,
                          vector<SolidObject> &structures,
                          SimParams &params) {
    
    // Get the fixed time step
    if (params.dtFixSet) {
        this->dtFix = params.dtFix;
    } else {
        this->dtFix = -1;
    }
    
    // Set the method order used for the evolution of the level set equation
    this->methodOrd = 3;

    this->nSteps = 0;

    // Set the repulsion mode and distance.
    this->repulseDist = params.repulseDist;
    this->collisionStiffness = params.collisionStiffness;
    this->collisionDist = params.collisionDist;
    this->repulseMode = params.repulseMode;

    // Create the kdTree data object and point cloud
    kdPointCloud = new KDPointCloud();
    kdTree = new KDTree(2, *kdPointCloud, KDTreeSingleIndexAdaptorParams(10));

    // Get the viscosity
    this->mu = params.mu;

    // Create the uniform meshes from the boundary object
    this->nx = params.nx;
    this->ny = params.ny;

    // MSS nx, ny
    int mssNx = params.mssNx;
    int mssNy = params.mssNy;

    // Create the arrays for the externally generated velocity field
    this->poolU = simutils::new_constant(ny+2*methodOrd, nx+2*methodOrd, 0.0);
    this->poolV = simutils::new_constant(ny+2*methodOrd, nx+2*methodOrd, 0.0);

    // this->nStructs = nStructs;
    this->nStructs = structures.size();

    // Create the pool labels
    this->pool = new objects::FSIObject*[ny];
    for (int i = 0; i < ny; i++) {
        this->pool[i] = new objects::FSIObject[nx];
    }

    // Set all cells to fluids to start.
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            this->pool[j][i] = objects::FLUID_C;
        }
    }

    // Generate the 1D meshes
    this->x = boundary.generateXMesh(this->nx);
    this->y = boundary.generateYMesh(this->ny);

    this->hx = x[1] - x[0];
    this->hy = y[1] - y[0];

    // Allocate the phi array and temp array for time stepping
    this->phi = simutils::new_constant(ny+2*methodOrd, nx+2*methodOrd, DINFINITY);
    this->phiReInit = simutils::new_constant(ny+2*methodOrd, nx+2*methodOrd, DINFINITY);
    this->phiRk1 = simutils::new_constant(ny+2*methodOrd, nx+2*methodOrd, 0.0);
    this->phiRk2 = simutils::new_constant(ny+2*methodOrd, nx+2*methodOrd, 0.0);

    // Apply the boundaries using the implicit SolidObject functions.
    // Also places the tracer particals
    if (nStructs > 0) {
        // Initialize the tracer particals
        this->tracers = new simstructs::tracer2D[nStructs];

        if (structures.size() != 0) {
            for (int l = 0; l < nStructs; l++) {
                this->embedShape(structures.at(l), l);
            }
        } else {
            assert(false);
        }
    }

    // Set the initial pool object velocities based on the information from the particles. Important
    // For assigning initial boundary values around the structure.
    if (nStructs == 1) {
        simutils::set_constant(ny+2*methodOrd, nx+2*methodOrd, tracers[0].u, this->poolU);
        simutils::set_constant(ny+2*methodOrd, nx+2*methodOrd, tracers[0].v, this->poolV);
    } else if (nStructs == 0) {
        simutils::set_constant(ny+2*methodOrd, nx+2*methodOrd, 0.0, this->poolU);
        simutils::set_constant(ny+2*methodOrd, nx+2*methodOrd, 0.0, this->poolV);
    }

    // Enumerate the Pool object
    this->enumeratePool();

    // Set up the domain membership matrix, first allocating it
    domainTracker = simutils::new_constant(ny, nx, DOMAIN_UNDESCOVERED);
    this->setUpDomainArray();

    // Create the medial axis point vector
    medialAxisPnts = new vector<tuple<double,double>>();

    // If we are using lower-dimensional MSS's
    solids = new vector<MassSpring2D>();
    if (nx != mssNx || ny != mssNy) {
        // Create different dimension Pool to represent the solid and copy the current Pool.
        SimParams temp = params;
        temp.setNx(mssNx);
        temp.setNy(mssNy);

        // Create new temp Pool
        Pool2D *poolTemp = new Pool2D(boundary, structures, temp);
        
        // Copy the Mass Spring objects from the old array to the new one
        for (int i = 0; i < nStructs; i++) {
            solids->push_back(MassSpring2D((poolTemp->solids)->at(i)));
        }

        // Interpolate the lower dimensional phi value into the higher dimensional space
        double x, y;
        for (int j = 0; j < ny; j++) {
            y = simutils::midpoint(this->y[j], this->y[j+1]);
            for (int i = 0; i < nx; i++) {
                x = simutils::midpoint(this->x[i], this->x[i+1]);
                phi[this->methodOrd+j][this->methodOrd+i] = poolTemp->interpolatePhi(x, y);
            }
        }
        this->enumeratePool();
        this->setUpDomainArray();

        // Delete the temp Pool
        delete poolTemp;
    } else {
        // Create the mass-spring system for each of the objects
        solids = new vector<MassSpring2D>();

        if (params.updateMode == 2 && !params.dtFixSet) {
            cout << "Attempting to use ADMM without fixed time step!" << endl;
            assert(false);
        }

        for (int i = 0; i < nStructs; i++) {
            solids->push_back(MassSpring2D(*this, i, structures.at(i),
                                params.updateMode, params.elementMode));

            if (params.updateMode == 2) {
                solids->at(i).setAdmmTol(params.admmTol);
            }
        }
    }

    // The state of each grid within the fast marching algorithm
    fastMarchingState = simutils::new_constant(ny, nx, FAR);

    // Locate nodes of MSS exactly on isocontour of new MSS
    // for (int i = 0; i < nStructs; i++) {
    //     solids->at(i).interpBoundary(*this, true);
    // }

    // Run the fast marching algorithm
    this->fastMarchPool(false, 0);

    simutils::copyVals(nx+2*methodOrd, ny+2*methodOrd, phiReInit, phi);

    this->enumeratePool();


    // Detect any collisions
    if (repulseMode == 2) {
        detectCollisions();
    }
}

/**
 * Use BFS to update the enumeration based on the current configuration of the MSS objects.
 * 
 * First it labels cells as solids depending on whether they intersect with the solid boundaries.
 * Then it labels the interfaces
*/
void Pool2D::bfsSearchForEdges(int structNum) {
    // Find the grid square that the tracer is currently on. This is the root node for the BFS.
    int i = simutils::findLimInfMeshPoint(tracers[structNum].x, this->x, this->nx+1);
    int j = simutils::findLimInfMeshPoint(tracers[structNum].y, this->y, this->ny+1);


    int ri, rj, ni, nj;
    // objects::FSIObject FL = objects::FLUID_C; // Make the ternary operator more palatable below

    // Create the queues for each of the indices (TODO: find a nicer solution than this)
    queue<int> queueX;
    queue<int> queueY;

    // Add the first node
    queueX.push(i);
    queueY.push(j);

    // Mark this node as a part of the structure. Otherwise, set it to "discovered"
    // this->domainTracker[j][i] = (objAtIndex(i, j) != FL) ? structNum : -1;
    pool[j][i] = objects::STRUCTURE;

    // double hx = x[1] - x[0];
    // double hy = y[1] - y[0];

    // Start the DFS
    double curPnt[2];
    double neighPnt[2];
    while (!queueX.empty()) {
        ri = queueX.front();
        rj = queueY.front();

        queueX.pop();
        queueY.pop();

        curPnt[0] = simutils::midpoint(x[ri], x[ri+1]);
        curPnt[1] = simutils::midpoint(y[rj], y[rj+1]);

        // Mark the current node as a structure
        // objAtIndex(ri, rj) = objects::STRUCTURE;

        // Discover the neighbouring points
        // cout << "ri = " << ri << " rj = " << rj << endl << endl;
        for (nj = rj-1; nj <= rj+1; nj++) {
            for (ni = ri-1; ni <= ri+1; ni++) {
                if (ni == ri && nj == rj) {
                    continue;
                }

                // If a neighbour is a labelled structure, move onto next iteration
                // cout << "ni = " << ni << " nj = " << nj << endl;
                if (objAtIndex(ni, nj) == objects::STRUCTURE || objAtIndex(ni, nj) == objects::UNLABELLED_INTERFACE) {
                    continue;
                }

                neighPnt[0] = simutils::midpoint(x[ni], x[ni+1]);
                neighPnt[1] = simutils::midpoint(y[nj], y[nj+1]);
                // neighPnt[0] = (ni - ri)*hx + curPnt[0];
                // neighPnt[1] = (nj - rj)*hy + curPnt[1];

                // If there is a structure intersection, we don't label the neighbour cell a structure.
                // If there is not an intersection, add this structural cell to the queue
                if (objAtIndex(ni, nj) == objects::UNLABELLED && !solids->at(structNum).intersectsBoundary(curPnt, neighPnt)) {
                    cout << "adding ni = " << ni << " nj = " << nj << endl;
                    queueX.push(ni);
                    queueY.push(nj);

                    this->pool[nj][ni] = objects::STRUCTURE;
                } else {
                    this->pool[nj][ni] = objects::UNLABELLED_INTERFACE;
                    // cout << "intersection found!" << endl;
                    // cout << "ni = " << ni << " nj = " << nj << endl;
                    // cout << "ri = " << ri << " rj = " << rj << endl;
                    // cout << "phi_n " << phi[methodOrd+nj][methodOrd+ni] << endl;
                    // cout << "phi_r " << phi[methodOrd+rj][methodOrd+ri] << endl;

                    // cout << "phi_n interp " << interpolatePhi(curPnt[0], curPnt[1]) << endl;
                    // cout << "phi_r interp " << interpolatePhi(neighPnt[0], neighPnt[1]) << endl;
                    // assert(phi[methodOrd+nj][methodOrd+ni] > 0);
                    // cout << "Intersection detected" << endl;
                }


                // if (this->domainTracker[nj][ni] == DOMAIN_UNDESCOVERED) { // If this node has not been seen
                //     if (objAtIndex(ni, nj) != FL) {
                //         // This node belongs to the structure
                //         queueX.push(ni);
                //         queueY.push(nj);
                //     } else {
                //         // This node belongs to the fluid region
                //         this->domainTracker[nj][ni] = DOMAIN_FLUID;
                //     }
                // }
            }
        }

        printPool();
        cout << "Press enter to continue... ";
        getchar();
        cout << endl;

    }
    cout << "out of loop" << endl;

    assert(queueX.empty() && queueY.empty());

    cout << "FINISHED search for edges" << endl;

}

/**
 * Reset the level set function using a bredth-first search from the tracer particals
*/
void Pool2D::resetLevelSet() {
    int i, j;
    // Change the object enumeration to all unlabelled
    for (j = 0; j < ny; j++) {
        for (i = 0; i < nx; i++) {
            pool[j][i] = objects::UNLABELLED;
        }
    }

    // For each structure, run the BFS
    for (i = 0; i < nStructs; i++) {
        bfsSearchForEdges(i);
    }

    // Using the known structure locations
    // Label the unlabelled interfaces
    for (j = 1; j < ny-1; j++) {
        // Higher level labels: interface and structure-based
        for (i = 1; i < nx-1; i++) {
            if (this->pool[j][i] == objects::UNLABELLED) {
                if (this->isStructConnected(i, j)) {
                    this->pool[j][i] = objects::UNLABELLED_INTERFACE;
                } else {
                    this->pool[j][i] = objects::FLUID_C;
                }
            }
        }
    }

    // Go through and label the unlabelled interfaces
    for (j = 1; j < ny-1; j++) {
        for (i = 1; i < nx-1; i++) {
            if (this->pool[j][i] == objects::UNLABELLED_INTERFACE) {
                this->labelInterface(i, j);
            }
        }
    }

    // Now run FMM to redistance (regenerate) the level set function
    fastMarchPool(false, 0);
    simutils::copyVals(nx+2*methodOrd, ny+2*methodOrd, phiReInit, phi);
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
 * 
 * TODO: refactor this code for better code re-use.
*/
void Pool2D::updateTracer(int structNum, double dt, int mode) {
    // Extract information from the tracer for this structure
    double tracerX = this->tracers[structNum].x;
    double tracerY = this->tracers[structNum].y;

    // Order of the numerical method, used for offsets
    int mo = methodOrd;

    // Find the "lim inf" mesh points for the current tracer location. Needed for all modes
    int i = simutils::findLimInfMeshPoint(tracers[structNum].x, this->x, this->nx+1);
    int j = simutils::findLimInfMeshPoint(tracers[structNum].y, this->y, this->ny+1);

    /* Noting that phi is defined at the midpoints of cell grids compute the stencil for the bilinear inteprolation.
        Use the right (up) shifted stencil unless the rightmost (topmost) cell is where the tracer lays, in this case we
        use the left-most stencil */
    
    int il = (i == nx-1) ? i-1 : i;
    int jl = (j == ny-1) ? j-1 : j;

    // The horizontal (x) stencil
    double xMesh[2] = {simutils::midpoint(this->x[il], this->x[il+1]),
                           simutils::midpoint(this->x[il+1], this->x[il+2])};

    // The vertical (y) stencil
    double yMesh[2] = {simutils::midpoint(this->y[jl], this->y[jl+1]),
                            simutils::midpoint(this->y[jl+1], this->y[jl+2])};

    // Compute the values array, again handling the special case where the tracer is in the topmost or rightmost cell
    double vals[4] = {phi[mo+jl][mo+il], phi[mo+jl][mo+il+1],
                            phi[mo+jl+1][mo+il], phi[mo+jl+1][mo+il+1]};

    // Use bi-linear interpolation to compute the value of phi, as well as its gradient
    double phi_ij = simutils::biLinearInterpolation(tracerX, tracerY, xMesh, yMesh, vals);

    double gradPhi_ij[2];
    simutils::biLinearGradient(tracerX, tracerY, xMesh, yMesh, vals, gradPhi_ij);
    simutils::normalize2D(gradPhi_ij);

    if (mode == 0) {
        /* Jump to interface mode */

        // Normalize the vector
        simutils::normalize2D(gradPhi_ij);

        tracers[structNum].x = tracerX - phi_ij*gradPhi_ij[0];
        tracers[structNum].y = tracerY - phi_ij*gradPhi_ij[1];
    } else if (mode == 1) {
        /* Velocity field jump mode */

        // Use bilinear interpolation to find the values of the velocity field at the tracer partical

        // u:
        double uVals[4] = {poolU[mo+jl][mo+il], poolU[mo+jl][mo+il+1],
                                poolU[mo+jl+1][mo+il], poolU[mo+jl+1][mo+il+1]};
        double uCur = simutils::biLinearInterpolation(tracerX, tracerY, xMesh, yMesh, uVals);

        // v:
        double vVals[4] = {poolV[mo+jl][mo+il], poolV[mo+jl][mo+il+1],
                                poolV[mo+jl+1][mo+il], poolV[mo+jl+1][mo+il+1]};
        double vCur = simutils::biLinearInterpolation(tracerX, tracerY, xMesh, yMesh, vVals);

        // Update the position of the tracer partical
        tracers[structNum].x += dt*uCur;
        tracers[structNum].y += dt*vCur;
    } else {
        /* Use gradient descent with line search and the interior perserving step limit
           to update the position of the tracer partical */
        
        this->interpolatePhiGrad(tracerX, tracerY, gradPhi_ij);
        phi_ij = this->interpolatePhi(tracerX, tracerY);
        
        
        // Use the closest point result to compute initial stepsize, ensuring that the updated point will be
        // within the negative part of the domain
        double alpha = 0.1;
        // double alpha = phi_ij/simutils::eucNorm2D(gradPhi_ij);

        // First gradient update
        double xStep = tracers[structNum].x - alpha*gradPhi_ij[0];
        double yStep = tracers[structNum].y - alpha*gradPhi_ij[1];

        // const double gamma = 0.5;
        const int max_iter = 10;
        int iters = 0;

        while (interpolatePhi(xStep, yStep) > phi_ij) {
            alpha /= 2;
            
            xStep = tracerX - alpha*gradPhi_ij[0];
            yStep = tracerY - alpha*gradPhi_ij[1];

            iters++;

            if (iters == max_iter) {
                break;
            }
        }

       
        if (iters < max_iter) {
            tracers[structNum].x = xStep;
            tracers[structNum].y = yStep;
        }
    }
}

/**
 * Embed the implcitly defined shapes into the pool and initalize the tracer
 * particals.
*/
void Pool2D::embedShape(SolidObject &struc, int structNum) {
    int i, j;
    int min_x, min_y;
    int mo = this->methodOrd;
    double minValue = phi[mo][mo];

    // Compute the signed distance function
    if (this->nStructs >= 0) {
        for (j = 0; j < this->ny; j++) {
            for (i = 0; i < this->nx; i++) {
                this->phiReInit[j+mo][i+mo] = struc.shapeEval(simutils::midpoint(this->x[i],
                    this->x[i+1]), simutils::midpoint(this->y[j], this->y[j+1]));
                
                if (this->phiReInit[j+mo][i+mo] < minValue) {
                    minValue = this->phiReInit[j+mo][i+mo];
                    min_x = i;
                    min_y = j;
                }
            }
        }
    }

    // Set the tracer partical for this level set function to initially be the minimum value found,
    // noting that it will be at the cell center for the ith cell
    if (this->nStructs >= 0) {
        this->tracers[structNum].x = simutils::midpoint(this->x[min_x-1], this->x[min_x]);
        this->tracers[structNum].y = simutils::midpoint(this->y[min_y-1], this->y[min_y]);
        this->tracers[structNum].mass = struc.getMass();
        this->tracers[structNum].u = struc.getU0(); // Velocity of the stucture, useful for update rules
        this->tracers[structNum].v = struc.getV0();
        this->tracers[structNum].isDeformable = struc.getObjType() == SolidObject::ObjectType::DEFORMABLE;
        this->isDeformable = this->isDeformable || struc.getObjType() == SolidObject::ObjectType::DEFORMABLE;
    }


    if (structNum == 0) {
        // If this is the first shape being embedded, simply copy into phi.
        simutils::copyVals(this->nx+2*methodOrd, this->ny+2*methodOrd, phiReInit, phi);
    } else {
        // If this is not the first structure, take the union of the current structure with the existing
        // level set.
        for (j = 0; j < this->ny; j++) {
            for (i = 0; i < this->nx; i++) {
                this->phi[j+mo][i+mo] = min( this->phi[j+mo][i+mo], this->phiReInit[j+mo][i+mo] );
            }
        }
    }
}

/**
 * Label a fluid interface according to the structure. Note that we only allow
 * interfaces that have at leasttwo neighbouring fluid cells. If it does not satisfy this
 * condition, we label it as a structure.
 * 
 * Note interfaces are defined as follows:
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
void Pool2D::labelInterface(int i, int j) {
    int n, e, s, w;
    // int mo = methodOrd;

    // Use the definition of the enum (prime factors for N, E, S, W) to
    // quickly decide the interface
    n = ( this->pool[j+1][i] != objects::STRUCTURE) ? objects::NORTH : 1;
    e = ( this->pool[j][i+1] != objects::STRUCTURE) ? objects::EAST  : 1;
    s = ( this->pool[j-1][i] != objects::STRUCTURE) ? objects::SOUTH : 1;
    w = ( this->pool[j][i-1] != objects::STRUCTURE) ? objects::WEST  : 1;

    int label = n*e*s*w;

    this->pool[j][i] = (objects::FSIObject)label;
}

/**
 * Enumerate the pool based on the current phi.
 * 
 * Boundaries of the rectangular domain are labelled in the obvious way.
 * All others are labelled adaptively using the signed distance function/
*/
void Pool2D::enumeratePool() {
    int i, j;

    // If there is not a signed distance function defined, set all internal points to fluids
    if (this->nStructs != 0) {
        // TODO: refactor this to be more efficient if possible
        for (j = 1; j < ny-1; j++) {
            // Higher level labels: interface and structure-based
            for (i = 1; i < nx-1; i++) {
                // Label the negative signed distance functions that are not being set to 0
                if (this->phi[j+this->methodOrd][i+this->methodOrd] < 0) {
                    this->pool[j][i] = objects::STRUCTURE;
                } else {
                    this->pool[j][i] = objects::FLUID_C;
                }
            }
        }

        // Label the unlabelled interfaces
        for (j = 1; j < ny-1; j++) {
            // Higher level labels: interface and structure-based
            for (i = 1; i < nx-1; i++) {
                if (this->pool[j][i] == objects::FLUID_C) {
                    if (this->isStructConnected(i, j)) {
                        this->pool[j][i] = objects::UNLABELLED_INTERFACE;
                    }
                }
            }
        }

        // Go through and label the unlabelled interfaces
        for (j = 1; j < ny-1; j++) {
            for (i = 1; i < nx-1; i++) {
                if (this->pool[j][i] == objects::UNLABELLED_INTERFACE) {
                    this->labelInterface(i, j);
                }
            }
        }
    }
}

/** 
 * Checks if there is a fluid cell in any of the cardinal directions. On structure cells to discover the
 * boundaries.
*/
bool Pool2D::isFluidConnected(int i, int j) {
    return this->pool[j][i+1] == objects::FLUID_C || this->pool[j][i-1] == objects::FLUID_C
            || this->pool[j+1][i] == objects::FLUID_C || this->pool[j-1][i] == objects::FLUID_C;
}

/** 
 * Checks if there is a structure in any of the cardinal directions from this point (only for boundary calls)
*/
bool Pool2D::isStructConnected(int i, int j) {
    return this->pool[j][i+1] == objects::STRUCTURE || this->pool[j][i-1] == objects::STRUCTURE
            || this->pool[j+1][i] == objects::STRUCTURE || this->pool[j-1][i] == objects::STRUCTURE;
}

/**
 * Get normal direction from interface point.
 * 
 * Think about the unit circle algorithm for this.
*/
void Pool2D::getNormalDir(objects::FSIObject obj, int nDir[2]) {
    nDir[0] = 0;
    nDir[1] = 0;

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
void tvdRK3HJ(double dt, double **arr, Pool2D *pool, int bcType,
                    double (Pool2D::*rhs)(int,int,double*,double*,double*,double*,int),
                    void (Pool2D::*bcFun)(double**)) {
    int i, j;
    int start, endx, endy;

    int mo = pool->methodOrd;

    double hx = pool->x[1] - pool->x[0];
    double hy = pool->y[1] - pool->y[0];

    // Start and end of the loop based on the type of BC
    if (bcType == 2) {
        start = 1; endx = pool->nx-2; endy = pool->ny-2;
    } else {
        start = 0; endx = pool->nx-1; endy = pool->ny-1;
    }

    // Stage 1 - apply BCs to the input array (necissary as boundary conditions
    // will depend on the order of the WENO scheme.
    if (bcType != 1) {
        (pool->*bcFun)(arr);
    }

    // Evaluate the RHS and store the result in the first temp array
    rhsHJENO(hx, hy, arr, pool->phiRk1, bcType, pool, rhs, bcFun);

    // Stage 1 Euler method
    for (j = start; j <= endy; j++) {
        for (i = start; i <= endx; i++) {
            // arr[mo+j][mo+i] = arr[mo+j][mo+i] + dt*pool->phiRk1[mo+j][mo+i];
            pool->phiRk2[mo+j][mo+i] = arr[mo+j][mo+i] + dt*pool->phiRk1[mo+j][mo+i];
        }
    }
    // return;

    // Stage 2
    if (bcType != 1) {
        (pool->*bcFun)(pool->phiRk2);
    }

    rhsHJENO(hx, hy, pool->phiRk2, pool->phiRk1, bcType, pool, rhs, bcFun);

    for (j = start; j <= endy; j++) {
        for (i = start; i <= endx; i++) {
            pool->phiRk2[mo+j][mo+i] = 0.75*arr[mo+j][mo+i]
                + 0.25*pool->phiRk2[mo+j][mo+i]
                + 0.25*dt*pool->phiRk1[mo+j][mo+i];
        }
    }

    // Stage 3 (solution attained)
    if (bcType != 1) {
        (pool->*bcFun)(pool->phiRk2);
    }

    rhsHJENO(hx, hy, pool->phiRk2, pool->phiRk1, bcType, pool, rhs, bcFun);

    for (j = start; j <= endy; j++) {
        for (i = start; i <= endx; i++) {
            arr[mo+j][mo+i] = (1.0/3.0)*arr[mo+j][mo+i]
                + (2.0/3.0)*pool->phiRk2[mo+j][mo+i]
                + (2.0/3.0)*dt*pool->phiRk1[mo+j][mo+i];
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
void rhsHJENO(double hx, double hy, double **in, double **out, int bcType, Pool2D *pool,
            double (Pool2D::*rhs)(int,int,double*,double*,double*,double*,int),
            void (Pool2D::*bcFun)(double**)) {
    int i, j, k;
    int start, endx, endy;

    int mo = pool->methodOrd;

    int nvals = 2*mo+1;
    double mxVals[nvals];
    double myVals[nvals];
    double xVals[nvals];
    double yVals[nvals];

    // Based on the boundary condition, choose how the RHS is updated on the boundary
    if (bcType == 0) {
        start = 0; endx = pool->nx-1; endy = pool->ny-1;
    } else {
        start = 1; endx = pool->nx-2; endy = pool->ny-2;
    }

    for (i = 0; i < nvals; i++) {
        mxVals[i] = ((double)i)*hx;
        myVals[i] = ((double)i)*hy;
    }

    // Compute the RHS array
    for (j = start; j <= endy; j++) {
        for (i = start; i <= endx; i++) {
            for (k = 0; k < nvals; k++) {
                yVals[k] = in[j+k][mo+i];
                xVals[k] = in[mo+j][i+k];
            }
            out[mo+j][mo+i] = (pool->*rhs)(i, j, mxVals, xVals, myVals, yVals, 3);
        }
    }

    // If using Neumann BCs, apply the condition on the spatial discretization.
    if (bcType == 1) {
        // cout << "Updating Nuemann BCs" << endl;
        (pool->*bcFun)(out);
    }
}

/**
 * RHS function for the level set function using the third order ENO scheme
*/
double Pool2D::levelSetRHS_ENO3(int i, int j, double *mxVals, double *xVals,
                                double *myVals, double *yVals, int maxOrd) {
    bool stencil[7] = {true, true, true, true, true, true, true};
    int mo = this->methodOrd;

    // double g[2];
    // this->levelSetGradient(i, j, g);

    int upwindU = (poolU[j+mo][i+mo] > 0) ? -1 : 1;
    int upwindV = (poolV[j+mo][i+mo] > 0) ? -1 : 1;

    double pntX = mxVals[3];
    double pntY = myVals[3];
    return -(poolU[j+mo][i+mo]*discs::thirdOrdENO(pntX, mxVals, xVals, upwindU, stencil, maxOrd)
                +   poolV[j+mo][i+mo]*discs::thirdOrdENO(pntY, myVals, yVals, upwindV, stencil, maxOrd));
}

/**
 * Checks if a point (i, j) is one grid cell away from the interface using the
 * methods from Russo and Smereka 2000.
 * 
 * NOTE: assumes that phi contains the original interface
*/
bool Pool2D::oneGridFromInterface(int i, int j) {
    int ox = methodOrd+i;
    int oy = methodOrd+j;
    return phi[oy][ox]*phi[oy][ox-1] < 0 || phi[oy][ox]*phi[oy][ox+1] < 0
        || phi[oy-1][ox]*phi[oy][ox] < 0 || phi[oy+1][ox]*phi[oy][ox] < 0; 
}

/**
 * Check if point is one away from an interface in the solid direction
 * 
 * NOTE: assumes that phi contains the original interface
*/
bool Pool2D::oneGridFromInterfaceStructure(int i, int j) {
    int ox = methodOrd+i;
    int oy = methodOrd+j;
    return phi[oy][ox] < 0 && (isInterface(objAtIndex(i-1, j)) || isInterface(objAtIndex(i+1, j))
            || isInterface(objAtIndex(i, j+1)) || isInterface(objAtIndex(i, j-1)));
}

/**
 * First-order approximation for the distance to the interface. Taken from Russo and Smereka 2000.
 * 
 * NOTE: assumes that phi contains the original interface before the beginning of time integration.
*/
double Pool2D::distanceFromInterface(int i, int j) {
    double hx, hy;
    double eps = 1e-16;
    int ox = methodOrd+i;
    int oy = methodOrd+j;

    hx = x[1] - x[0];
    hy = y[1] - y[0];

    // Compute centered finite difference discretization of the first derivatives
    double phiGrad[2];
    phiGrad[0] = (phi[oy][ox+1] - phi[oy][ox-1])/(2.0*hx);
    phiGrad[1] = (phi[oy+1][ox] - phi[oy-1][ox])/(2.0*hy);

    // Take the Euclidian norm of the gradient
    double phi_norm = simutils::eucNorm2D(phiGrad);

    // Compute the distance approximation, with a small epsilon in the denom to avoid
    // dividing by 0.
    return phi[oy][ox]/(phi_norm+eps);
}

double Pool2D::signedDistanceReinitializationRHS_CorrectUpwind(int i, int j, double *mxVals,
                                                      double *xVals, double *myVals, 
                                                      double *yVals) {
    
    // Compute the gradient vector
    int s = 3;
    int l = 2;
    int r = 4;
    double mo = methodOrd;
    int ox = mo+i;
    int oy = mo+j;

    double hx, hy, h;

    hx = x[1] - x[0];
    hy = y[1] - y[0];
    h = sqrt(simutils::square(hx)+simutils::square(hy));

    // Compute the gradient using scheme from Russo and Smerka
    double G = 0;
    double a, b, c, d;

    
    if (this->oneGridFromInterface(i, j)) {
        // Use interface-conscious upwind scheme
        return (-this->sign(phi[oy][ox])*abs(xVals[s]) + this->distanceFromInterface(i, j))/h;
    } else {
        // Compute the upwind discretization
        if (phi[oy][ox] > 0) {
            a = simutils::dmax((xVals[s] - xVals[l])/hx, 0.0);
            b = simutils::dmin((xVals[r] - xVals[s])/hx, 0.0);
            c = simutils::dmax((yVals[s] - yVals[l])/hy, 0.0);
            d = simutils::dmin((yVals[r] - yVals[s])/hy, 0.0);

            G = sqrt(simutils::dmax(simutils::square(a), simutils::square(b))
                + simutils::dmax(simutils::square(c), simutils::square(d))) - 1.0;
        } else if (phi[oy][ox] < 0) {
            a = simutils::dmin((xVals[s] - xVals[l])/hx, 0.0);
            b = simutils::dmax((xVals[r] - xVals[s])/hx, 0.0);
            c = simutils::dmin((yVals[s] - yVals[l])/hy, 0.0);
            d = simutils::dmax((yVals[r] - yVals[s])/hy, 0.0);

            G = sqrt(simutils::dmax(simutils::square(a), simutils::square(b))
                + simutils::dmax(simutils::square(c), simutils::square(d))) - 1.0;
        } else {
            G = 0;
        }
        return -this->sign(phi[oy][ox])*G;
    }
}

/**
 * NOTE: depends on phi remaining constant through the whole integration.
 * 
 * TODO: this method does not work, perhaps due to incorrect tracking of the interface.
 *       Additionally, there may be issues with using a non-conservative finite difference 
 *       scheme for this task, as oscillations seem to appear.
 * 
 * TODO: use the first order method near the interface when i, j is available
 * 
*/
double Pool2D::signedDistanceReinitializationRHS_ENO3(int i, int j, double *mxVals,
                                                      double *xVals, double *myVals, 
                                                      double *yVals) {
    
    // Compute the gradient vector
    double h, hx, hy, sign;

    int mo = methodOrd;
    int ox = mo+i;
    int oy = mo+j;
    int s = 3; // Point in (x|y)Vals containing phi[j][i]

    hx = x[1] - x[0];
    hy = y[1] - y[0];
    h = (hx < hy) ? hx : hy;

    if (this->oneGridFromInterface(i, j)) {
        // Use interface-conscious upwind scheme
        return (-this->sign(phi[oy][ox])*abs(xVals[s]) + this->distanceFromInterface(i, j))/h;
    } else {
        bool xStencil[7];
        bool yStencil[7];

        // Create the widest stencils possible, avoiding the ghost cells.
        this->getWidestStencil(0, nx-1, i, xStencil);
        this->getWidestStencil(0, ny-1, j, yStencil);

        int upwindX = 1;
        int upwindY = 1;

        double pnt = xVals[3];
        double phi_x = discs::thirdOrdENO(pnt,mxVals,xVals,upwindX,xStencil,3);
        double phi_y = discs::thirdOrdENO(pnt,myVals,yVals,upwindY,yStencil,3);

        // Compute the smeared signed distance function
        sign = smearSign(h, phi[methodOrd+j][methodOrd+i]);

        // Calculate the RHS of the HJ equation using the third order ENO finite difference scheme.
        return sign*(1.0 - sqrt( simutils::square(phi_x)
                        +  simutils::square(phi_y) ));
    }
}

/**
 * Compute the gradient of the level set using upwind finite differencing.
 * 
 * TODO: add higher order gradient approximation
*/
void Pool2D::levelSetGradient(int i, int j, double g[2]) {
    double hx, hy;
    double phi_min;

    int mo = this->methodOrd;
    
    // Mesh spacing in each dimension
    hx = x[1] - x[0];
    hy = y[1] - y[0];

    bool ip1Less = phi[mo+j][mo+i+1] < phi[mo+j][mo+i-1];
    phi_min = (ip1Less) ? phi[mo+j][mo+i+1] : phi[mo+j][mo+i-1];

    /* Computing the gradient */

    // Take the x partial
    if (ip1Less) {
        g[0] = (phi[mo+j][mo+i+1] - phi[mo+j][mo+i])/hx;
    } else {
        g[0] = (phi[mo+j][mo+i] - phi[mo+j][mo+i-1])/hx;
    }

    // Take the y partial
    ip1Less = phi[mo+j+1][mo+i] < phi[mo+j-1][mo+i];
    phi_min = (ip1Less) ? phi[mo+j+1][mo+i] : phi[mo+j-1][mo+i];

    // if (phi_min < phi[mo+j][mo+i]) {
    if (ip1Less) {
        g[1] = (phi[mo+j+1][mo+i] - phi[mo+j][mo+i])/hy;
    } else {
        g[1] = (phi[mo+j][mo+i] - phi[mo+j-1][mo+i])/hy;
    }
}

/** 
 * Compute the outward normal vector of the level set function at grid point (i, j)
 * (assumed to be an interior point)
*/
void Pool2D::levelSetUnitNormal(int i, int j, double n[2]) {

    // Take the gradient
    this->levelSetGradient(i, j, n);

    // Normalize the gradient
    double scal = sqrt(simutils::square(n[0]) + simutils::square(n[1]));

    if (scal > 0) {
        n[0] /= scal;
        n[1] /= scal;
    }
}

/**
 * Use the two-way velocity extrapolation to extrapolate the velocity into the level set to
 * set constant values along the normal of the level set. Note: two-way formula is more complicated
 * but has the advantage 
 * 
 * Note: this only uses the interior points, i.e., i-1, i+1, j-1, j+1 will never be ghost cells. 
*/
double Pool2D::velocityExtrapolationRHS_ENO3(int i, int j, double *mxVals,
                                             double *xVals, double *myVals, 
                                             double *yVals) {

    // Backup plan
    if (objAtIndex(i, j) != objects::FLUID_C) {
        return 0.0;
    }
    double sign;
    bool stencil[7] = {true, true, true, true, true, true, true};
    

    // bool xStencil[7];
    // bool yStencil[7];
    // double lvlNormal[2];
    int mo = this->methodOrd;

    // Create the widest stencils possible for the ENO schemes, avoiding the ghost cells.
    // this->getWidestStencil(0, nx-1, i, xStencil);
    // this->getWidestStencil(0, ny-1, j, yStencil);

    double hx, hy, h;
    hx = x[1] - x[0];
    hy = y[1] - y[0];

    h = sqrt(simutils::square(hx)+simutils::square(hy));

    // Compute the sign function
    sign = this->sign(phi[mo+j][mo+i]);

    // Compute the normal vector for the level set at the grid point (i, j)
    double phiGrad[2] = {0.0, 0.0};

    if (i == 0) {
        phiGrad[0] = (-(3.0/2.0)*phi[mo+j][mo+i] + 2.0*phi[mo+j][mo+i+1] - 0.5*phi[mo+j][mo+i+2])/hx;
    } else if (i == nx-1) {
        phiGrad[0] = ((3.0/2.0)*phi[mo+j][mo+i] - 2.0*phi[mo+j][mo+i-1] + 0.5*phi[mo+j][mo+i-2])/hx;
    } else {
        phiGrad[0] = (phi[mo+j][mo+i+1] - phi[mo+j][mo+i-1])/(2.0*hx);
    }

    // Compute gradient of phi using upwind finite differences
    // double phi_l = phi[mo+j][mo+i-1];
    // double phi_c = phi[mo+j][mo+i];
    // double phi_r = phi[mo+j][mo+i+1];
    // if (phi_l < phi_c || phi_r < phi_c) {
    //     if (phi_l < phi_r) {
    //         phiGrad[0] = (phi_c - phi_l)/hx;
    //     } else {
    //         phiGrad[0] = (phi_r - phi_c)/hx;
    //     }
    // } else {
    //     // phiGrad[0] = 0.0;
    //     phiGrad[0] = (phi_r - phi_l)/(2.0*hx);
    // }

    // phi_l = phi[mo+j-1][mo+i];
    // phi_c = phi[mo+j][mo+i];
    // phi_r = phi[mo+j+1][mo+i];
    // if (phi_l < phi_c || phi_r < phi_c) {
    //     if (phi_l < phi_r) {
    //         phiGrad[1] = (phi_c - phi_l)/hy;
    //     } else {
    //         phiGrad[1] = (phi_r - phi_c)/hy;
    //     }
    // } else {
    //     // phiGrad[1] = 0.0;
    //     phiGrad[1] = (phi_r - phi_l)/(2.0*hy);
    // }

    if (j == 0) {
        phiGrad[1] = (-(3.0/2.0)*phi[mo+j][mo+i] + 2.0*phi[mo+j+1][mo+i] - 0.5*phi[mo+j+2][mo+i])/hy;
    } else if (j == ny-1) {
        phiGrad[1] = ((3.0/2.0)*phi[mo+j][mo+i] - 2.0*phi[mo+j-1][mo+i] + 0.5*phi[mo+j-2][mo+i])/hy;
    } else {
        phiGrad[1] = (phi[mo+j+1][mo+i] - phi[mo+j-1][mo+i])/(2.0*hy);
    }

    // Normalize the gradient
    simutils::normalize2D(phiGrad);

    int upwindX = (sign*phiGrad[0] > 0)? -1 : 1;
    int upwindY = (sign*phiGrad[1] > 0)? -1 : 1;

    // Calculate the RHS of the HJ equation using the third order ENO finite difference scheme.
    double pnt = xVals[3];
    return -( phiGrad[0]*discs::thirdOrdENO(pnt,mxVals,xVals,upwindX,stencil,3)
                    +  phiGrad[1]*discs::thirdOrdENO(pnt,myVals,yVals,upwindY,stencil,3) );
}

/**
 * Get the largest possible sided stencil that can be used for index i, j, where cell (i,j)
 * is assumed to be a fluid cell (objects::FLUID_C).
 * 
 * NOTE: assumes that there are properly set ghost cells
 * 
 * Different behaviour occurs depending on the axis used:
 *     axis = 0: the x-axis, using u velocity update logic on the staggered grid
 *     axis = 1: the y-axis, using y velocity update logic on the staggered grid
 * 
 * OUT: stencil, loaded with true values for each point on the mesh
 *          X_{k-methodOrd} < ... < X_{k} < ... <  X_{k+methodOrd},
 *      that is a valid stencil point along the specified axis
 * 
 * Returns: the size of the largest possible stencil
*/
int Pool2D::getPossibleStencil(int i, int j, int axis, int methodOrd, bool *stencil) {//, bool ghostSkip) {
    

    for (int k = 0; k < 2*methodOrd+1; k++) {
        if (axis == 0) {
            stencil[k] = isUsableU(i-methodOrd+k, j);///this->objAtIndex(i-methodOrd+k, j) == objects::FLUID_C;
        } else if (axis == 1) {
            // cout << "j = " << j-methodOrd+k << " ";
            stencil[k] = isUsableV(i, j-methodOrd+k);//this->objAtIndex(i, j-methodOrd+k) == objects::FLUID_C;
        }
    }
    // if (axis == 1) {
    //     cout << endl;
    // }
    return simutils::sum(2*methodOrd+1, (int*)stencil);
}


/**
 * Get the widest stencil between min_index and max_index.
 * 
 * Useful for getting a sided stencil when dealing the boundary conditions.
*/
void Pool2D::getWidestStencil(int minIndex, int maxIndex, int curIndex, bool *stencil) {
    int off;

    for (int k = 0; k < 2*methodOrd+1; k++) {
        off = curIndex - methodOrd + k;

        if (off < minIndex || off > maxIndex) {
            stencil[k] = false;
        } else {
            stencil[k] = true;
        }
    }
}

void Pool2D::applyLevelSetPeriodicBCs(double **phi) {
    int i, j;
    int mo = this->methodOrd;

    // Apply a periodic boundary condition using the existing pool
    for (j = 0; j < ny; j++) {
        for (i = 0; i < mo; i++) {
            // LEFT
            phi[j+mo][i] = phi[j+mo][nx+i];

            // RIGHT
            phi[j+mo][mo+nx+i] = phi[j+mo][mo+i];
        }
    }

    for (i = 0; i < this->nx; i++) {
        for (j = 0; j < mo; j++) {
            // BOT
            phi[j][i+mo] = phi[ny+j][i+mo];

            // TOP
            phi[mo+ny+j][i+mo] = phi[mo+j][i+mo];
        }
    }
}

/**
 * Apply the reinitialization boundary conditions mentioned in the Fatori thesis.
*/
void Pool2D::applyReinitializationBCs(double **phi) {
    int i, j;

    int mo = this->methodOrd;

    // Sides
    for (j = 0; j <= ny-1; j++) {
        phi[mo+j][mo] = phi[mo+j][mo+1]
            + this->sign(phi[mo+j][mo+1])*abs(phi[mo+j][mo+1] - phi[mo+j][mo+2]); // LHS
        
        phi[mo+j][mo+(nx-1)] = phi[mo+j][mo+(nx-1)-1]
            + this->sign(phi[mo+j][mo+(nx-1)-1])*abs(phi[mo+j][mo+(nx-1)-1] - phi[mo+j][mo+(nx-1)-2]); // RHS
    }

    for (i = 0; i <= nx-1; i++) {
        phi[mo][mo+i] = phi[mo+1][mo+i]
            + this->sign(phi[mo+1][mo+i])*abs(phi[mo+1][mo+i] - phi[mo+2][mo+i]); // BOT

        phi[mo+(ny-1)][mo+i] = phi[mo+(ny-1)-1][mo+i]
            + this->sign(phi[mo+(ny-1)-1][mo+i])*abs(phi[mo+(ny-1)-1][mo+i] - phi[mo+(ny-1)-2][mo+i]); // TOP

    }
}

/**
 * Apply Homogenous Nuemann conditions to directional derivatives on the boundary. The "no-slip"
 * condition
*/
void Pool2D::applyHomogenousNeumann(double **arr) {
    int i, j;

    int mo = this->methodOrd;

    // Sides
    for (j = 0; j <= ny-1; j++) {
        arr[mo+j][mo] = 0.0;
        arr[mo+j][mo+(nx-1)] = 0.0;
    }

    // Top and bot
    for (i = 0; j <= nx-1; i++) {
        arr[mo][mo+i] = 0.0;
        arr[mo+(ny-1)][mo+i] = 0.0;
    }
}

/**
 * The numerically smeared sign function for the reinitialization equation.
 * 
 * Takes in a grid spacing h and a scalar value of a signed distance function phi.
*/
double Pool2D::phiSmearSign(double h, double phi, double phi_x, double phi_y) {
    return phi/sqrt(simutils::square(phi) + simutils::square(h)*(simutils::square(phi_x)+simutils::square(phi_y)));
}

/**
 * The numerically smeared sign function for the reinitialization equation.
 * 
 * Takes in a grid spacing h and a scalar value of a signed distance function phi.
*/
double Pool2D::smearSign(double h, double phi) {
    return phi/sqrt(simutils::square(phi) + simutils::square(h));
}

/**
 * Sign function.
*/
double Pool2D::sign(double phi) {
    return (phi >= 0.0) ? 1.0 : -1.0;
}

/* PUBLIC METHODS */

/**
 * Note that nx, ny, from the perspective of the pool, is the number of subintervals.
*/
Pool2D::Pool2D(Boundary &boundary, vector<SolidObject> &structures,
               SimParams &params) {
    // Create the pool based on the given implementation
    this->create2DPool(boundary, structures, params);
}

/**
 * Note: fNet on the tracers should be reset in the calling function
 * 
*/
void Pool2D::calculateLocalFNet(int i, int j, objects::FSIObject obj, int ng,
                                double **u, double **v, double **p) {
    // Reset dA
    double dA = 0.0;
    double n[2] = {0.0, 0.0};
    double uGrad[2];
    double vGrad[2];
    int is, js, is_mv, js_mv;
    double hx = x[1] - x[0];
    double hy = y[1] - y[0];
    double pres;

    // Set default stencil
    is = i; is_mv = i+1;
    js = j; js_mv = j+1;

    if (this->hasStructInDir(obj, objects::EAST)) {
        is = i-1;
        is_mv = i-2;

        n[0] -= 1.0;
        dA += simutils::square(hx);
    } else if (this->hasStructInDir(obj, objects::WEST)) {
        n[0] += 1.0;
        dA += simutils::square(hx);
    }

    if (this->hasStructInDir(obj, objects::NORTH))  {
        js = j-1;
        js_mv = j-2;

        n[1] -= 1.0;
        dA += simutils::square(hy);
    } else if (this->hasStructInDir(obj, objects::SOUTH)) {
        n[1] += 1.0;
        dA += simutils::square(hy);
    }

    dA = sqrt(dA);

    // Compute the information for the hydrodynamic stress tensor.
    // By taking the value the unit normal points to.
    pres = p[ng+j+(int)n[1]][ng+i+(int)n[0]];

    // Use either backwards or forwards difference formula to compute the derivatives
    // in each direction

    // x derivatives
    if (is < is_mv) {
        // Forward difference
        uGrad[0] = (u[ng+j][ng+is_mv] - u[ng+j][ng+is])/hx;
        vGrad[0] = (v[ng+js][ng+is_mv] - u[ng+js][ng+is])/hx;
    } else {
        // Backward difference
        uGrad[0] = (u[ng+j][ng+is] - u[ng+j][ng+is_mv])/hx;
        vGrad[0] = (u[ng+js][ng+is] - u[ng+js][ng+is_mv])/hx;
    }

    // y derivatives
    if (js < js_mv) {
        uGrad[1] = (u[ng+js_mv][ng+is] - u[ng+js][ng+is])/hy;
        vGrad[1] = (v[ng+js_mv][ng+i] - v[ng+js][ng+i])/hy;
    } else {
        uGrad[1] = (u[ng+js][ng+is] - u[ng+js_mv][ng+is])/hy;
        vGrad[1] = (v[ng+js][ng+i] - v[ng+js_mv][ng+i])/hy;
    }

    // cout << "n is normalized now" << endl;
    simutils::normalize2D(n);
    // assert(false);

    tracers[domainMembership(i, j)].fNetX += ((-pres + 2.0*this->mu*uGrad[0])*n[0]
        + this->mu*(uGrad[1] + vGrad[0])*n[1])*dA; 
    tracers[domainMembership(i, j)].fNetY += this->mu*((uGrad[1] + vGrad[0])*n[0]
        + (-pres + 2.0*this->mu*vGrad[1])*n[1])*dA;
    // }
}

/** 
 * Method which computes the values of the stresses at each of the boundary points.
 * Note the relationship
 * 
 * Note: i, j has to be on the structure domain.
*/
void Pool2D::computeBoundaryStress(int i, int j, objects::FSIObject obj, int ng,
                                double **u, double **v, double **p) {
    // Reset dA
    double n[2] = {0.0, 0.0};
    double uGrad[2];
    double vGrad[2];
    int is, js, is_mv, js_mv;
    double hx = x[1] - x[0];
    double hy = y[1] - y[0];
    double pres;
    int mo = this->methodOrd;

    // Set default stencil
    is = i; is_mv = i+1;
    js = j; js_mv = j+1;

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

    // Compute the information for the hydrodynamic stress tensor.
    // By taking the value the unit normal points to.
    pres = p[ng+j+(int)n[1]][ng+i+(int)n[0]];

    // Use either backwards or forwards difference formula to compute the derivatives
    // in each direction

    if (is < is_mv) {
        // Forward difference
        uGrad[0] = (u[j][is_mv] - u[j][is])/hx;
        vGrad[0] = (v[js][is_mv] - u[js][is])/hx;
    } else {
        // Backward differenceß
        uGrad[0] = (u[j][is] - u[j][is_mv])/hx;
        vGrad[0] = (u[js][is] - u[js][is_mv])/hx;
    }

    // y derivatives
    if (js < js_mv) {
        uGrad[1] = (u[js_mv][is] - u[js][is])/hy;
        vGrad[1] = (v[js_mv][i] - v[js][i])/hy;
    } else {
        uGrad[1] = (u[js][is] - u[js_mv][is])/hy;
        vGrad[1] = (v[js][i] - v[js_mv][i])/hy;
    }


    // cout << "N is normalized now!!" << endl;
    simutils::normalize2D(n);
    // assert(false);

    // Compute the stresses (assigning them to the pool velocities temporarily for memory saving)
    poolU[mo+j][mo+i] = ((-pres + 2.0*this->mu*uGrad[0])*n[0]
        + this->mu*(uGrad[1] + vGrad[0])*n[1]); 
    poolV[mo+j][mo+i] = this->mu*((uGrad[1] + vGrad[0])*n[0]
        + (-pres + 2.0*this->mu*vGrad[1])*n[1]);
}

/**
 * Set up the Pool velocity field for the case where we have multiple immersed structures.
*/
void Pool2D::multiStructureVelocitySet(double **u, double **v, int ng) {
    int i, j;
    int mo = this->methodOrd;
    objects::FSIObject obj;
    
    cout << "Setting the pool velocity field. Printing tracer vels" << endl;
    for (i = 0; i < this->nStructs; i++) {
        cout << "tracer " << i << " u = " << tracers[i].u << " v = " << tracers[i].v << endl;
    }

    for (j = 0; j < this->ny; j++) {
        for (i = 0; i < this->nx; i++) {
            obj = objAtIndex(i, j);

            if (obj != objects::FLUID_C) {
                this->poolU[j+mo][i+mo] = tracers[domainMembership(i, j)].u;
                this->poolV[j+mo][i+mo] = tracers[domainMembership(i, j)].v;
            } else {
                this->poolU[j+mo][i+mo] = 0.0;
                this->poolV[j+mo][i+mo] = 0.0;
            }
        }
    }

    this->fastMarchPool(true, 0);
}

/**
 * Function to calculate the velocities of the objects in the pool.
 * Note: has only been implemented in the trivial case.
 * 
 * These object velocities are applied to each of the tracer particals.
 * 
 * Using the method from Justin's brain paper. Assuming that the tracking particals are
 * within the respective level set functions.
 * 
 * Note: poolU and poolV are used as temp arrays in this computation
*/
void Pool2D::updatePoolVelocities(double dt, double **u, double **v, double **p, int ng) {
    int i, j;
    int mo = methodOrd;
    objects::FSIObject obj;

    if (this->nStructs >= 1) {

        // Re-initialize the tracer partical forces to 0
        for (i = 0; i < nStructs; i++) {
            tracers[i].fNetX = 0.0;
            tracers[i].fNetY = 0.0;
        }

        for (j = 1; j < ny-1; j++) {
            for (i = 1; i < nx-1; i++) {
                // Get the object at the current cell
                obj = this->objAtIndex(i, j);

                // If this is an interface point, calculate its contribution to the net force being
                // applied to the structure it belongs to
                if (this->isInterface(obj)) {
                    if (!this->isDeformable) {
                        this->calculateLocalFNet(i, j, obj, ng, u, v, p);
                    } else {
                        this->computeBoundaryStress(i, j, obj, ng, u, v, p);
                    }
                } else {
                    this->poolU[mo+j][mo+i] = 0.0;
                    this->poolV[mo+j][mo+i] = 0.0;
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

                for (int k = 0; k < solids->at(i).pntList->size(); k++) {
                    (*solids->at(i).qt)[2*k]   = tracers[i].u;
                    (*solids->at(i).qt)[2*k+1] = tracers[i].v;

                    (*solids->at(i).q)[2*k]   += (*solids->at(i).qt)[2*k]*dt;
                    (*solids->at(i).q)[2*k+1] += (*solids->at(i).qt)[2*k+1]*dt;

                    solids->at(i).pntList->at(k).x = (*solids->at(i).q)[2*k];
                    solids->at(i).pntList->at(k).y = (*solids->at(i).q)[2*k+1];

                    solids->at(i).pntList->at(k).u = (*solids->at(i).qt)[2*k];
                    solids->at(i).pntList->at(k).v = (*solids->at(i).qt)[2*k+1];

                }
            }

            // Move the mass spring system nodes along with the tracer

            // The the existing pool velocity arrays to 0. The velocity values
            // are already encoded into the tracer particals.
            if (nStructs <= 1) {
                simutils::set_constant(ny+2*mo, nx+2*mo, tracers[0].u, this->poolU);
                simutils::set_constant(ny+2*mo, nx+2*mo, tracers[0].v, this->poolV);
            } else {
                this->multiStructureVelocitySet(u, v, ng);
            }
        } else {
            // For each of the structures, update the deformable solids, compute the net force,
            // and update the positions of the tracer particals.
            double **stress[2] = {poolU, poolV};

            if (repulseMode != 0) {
                // TODO: maybe think about a more efficient 
                this->fastMarchPool(false, 2);
                this->detectCollisions();
                this->setUpDomainArray();
            }

            // assert(false);

            double fNet[2] = {0, 0};
            // int elementMode = 0;
            for (int i = 0; i < nStructs; i++) {
                // cout << "OBJ " << i << endl;
                solids->at(i).updateSolidVels(dt, *this, stress, fNet, mo, false);

                // Update the tracer info (not particularly necissary perhaps but will be nice for testing)
                tracers[i].fNetX = fNet[0];
                tracers[i].fNetY = fNet[1];

                tracers[i].u = tracers[i].u + ((tracers[i].fNetX)/(tracers[i].mass))*dt;
                tracers[i].v = tracers[i].v + ((tracers[i].fNetY)/(tracers[i].mass))*dt;

                fNet[0] = 0.0;
                fNet[1] = 0.0;
            }
            // assert(false);
            // cout << "in pool finished updating all structs" << endl;
            // assert(false);

            //

            // Now that we have the velocities, we can apply the fast marching algorithm.
            // if (repulseMode != 0) {
            //     this->fastMarchPool(true);
            // } else {
            //     this->fastMarchPool(true);
            //     this->detectCollisions();
            //     this->fastMarchPool(true);
            // }
            this->fastMarchPool(true, 0);
        }
        
    } else {
        simutils::set_constant(ny+2*mo, nx+2*mo, 0.0, this->poolU);
        simutils::set_constant(ny+2*mo, nx+2*mo, 0.0, this->poolV);
    }
}

/**
 * Update the pool using the third order ENO scheme. Using periodic BC's for simplicity.
 * 
 * ng := number of ghost cells on the boundary.
 * 
*/
void Pool2D::updatePool(double dt, double **u, double **v, double **p, int ng, bool reinitialize) {
    int k;

    // If this is called but there are no immersed structures, simply return
    if (this->nStructs < 1) {
        return;
    }

    // Set up the domain array
    setUpDomainArray();

    // Update the velocity field
    this->updatePoolVelocities(dt, u, v, p, ng);

    if (reinitialize) {
        simutils::copyVals(nx+2*methodOrd, ny+2*methodOrd, phiReInit, phi);

        // Update the enumeration and domain array
        // TODO: are these necissary?
        // enumeratePool();
        // setUpDomainArray();
    }

    // Update the pool using the third order ENO scheme to test level set in the equation.
    // Using periodic boundary condition.
    tvdRK3HJ(dt, phi, this, 0, &Pool2D::levelSetRHS_ENO3,
             &Pool2D::applyLevelSetPeriodicBCs);

    // Update the positions of the tracer particals
    for (k = 0; k < nStructs; k++) {
        updateTracer(k, dt, 2);
    }
    
    // Update the enumeration and domain array for the shifted level set.
    enumeratePool();
    setUpDomainArray();

    // Update the position of the solids
    // cout << "Updating the solid locations" << endl;
    if (isDeformable) {
        for (int i = 0; i < nStructs; i++) {
            solids->at(i).updateSolidLocs(*this, false);
        }
    }

    // If significant drift from the interface detected, correct.
    // if (shouldRefitSDF(min(hx, hy))) {
    //     refitToSolids(ng);
    // }

    nSteps++;
}

/**
 * Checks if the MSS and the signed distance function remain well asigned.
 * If the interface has drifted too far away from the MSS boundarey edges,
 * we apply the squeeze algorithm to bring it back into place and then reinitialize
 * as a signed distance function.
*/
bool Pool2D::shouldRefitSDF(double tol) {
    // For each MSS node, Find the value of the level set function at this point.
    // If it is too far from the level set interface, we require squeeze.
    bool squeeze = false;
    for (auto solid = solids->begin(); solid != solids->end(); ++solid) {
        for (auto pnt = solid->pntList->begin(); pnt != solid->pntList->end(); ++solid) {
            squeeze = abs(interpolatePhi(pnt->x, pnt->y)) < tol;

            if (!squeeze) {
                return squeeze;
            }
        }
    }

    return squeeze;
}

/**
 * Build squeeze field to extrapolate the velocity. Return the CFL constant
*/
double Pool2D::buildSqeezeField() {
    int mo = this->methodOrd;
    double pnt[2];
    double near[2];
    double maxU = -1;
    double maxV = -1;
    // double maxDist = -1;
    // double dist;

    enumeratePool();
    setUpDomainArray();

    // Make the pool velocities on the interface the distance vector to the MSS
    for (int j = 0; j < this->ny; j++) {
        pnt[1] = simutils::midpoint(y[j], y[j+1]);
        for (int i = 0; i < this->nx; i++) {
            pnt[0] = simutils::midpoint(x[i], x[i+1]);

            if (isInterface(objAtIndex(i, j))) {
                solids->at(domainMembership(i, j)).closestBoundaryPnt(pnt, near);

                poolU[mo+j][mo+i] = near[0] - pnt[0];
                poolV[mo+j][mo+i] = near[1] - pnt[1];

                maxU = (abs(poolU[mo+j][mo+i]) > maxU) ? abs(poolU[mo+j][mo+i]) : maxU;
                maxV = (abs(poolV[mo+j][mo+i]) > maxV) ? abs(poolV[mo+j][mo+i]) : maxV;
            }
        }
    }

    // Extrapolate the velocity field
    fastMarchPool(true, 1);

    double hx = x[1] - x[0];
    double hy = y[1] - y[0];
    double CFL = 1.0 / ((maxU / hx) + (maxV / hy));

    return CFL;

}

/**
 * Use a simple shrinking algorithm to refit the level set function ot the meshes
 * 
 * TODO: may need to reinitialize before using this
*/
void Pool2D::refitToSolids(int ng) {
    int mo = this->methodOrd;

    double t = 0;

    double CFL = buildSqeezeField();

    while (t < 1.0) {
        // Update the pool using the third order ENO scheme to test level set in the equation.
        // Using periodic boundary condition.
        tvdRK3HJ(CFL, phi, this, 0, &Pool2D::levelSetRHS_ENO3,
                &Pool2D::applyLevelSetPeriodicBCs);
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
    fastMarchPool(false, 0);

    simutils::copyVals(nx+2*methodOrd, ny+2*methodOrd, phiReInit, phi);
}

/* IO */
//------------------------------

void Pool2D::printPool() {
    for (int j = 0; j < this->ny; j++) {
        for (int i = 0; i < this->nx; i++) {
            cout << (int)this->pool[j][i] << ", ";
        }
        cout << endl;
    }
    cout << endl;
}

void Pool2D::printTracers() {
    for (int i = 0; i < this->nStructs; i++) {
        cout << "Tracer " << i << ": x = " << this->tracers[i].x << " y = " << this->tracers[i].y;
        cout << endl;
    }
}

void Pool2D::outputMedialAxisApprox(const char *fname) {
    double x, y;

    ofstream outFile;
    outFile.open(fname);

    for (auto tup = medialAxisPnts->begin(); tup != medialAxisPnts->end(); ++tup) {
        x = get<0>(*tup);
        y = get<1>(*tup);

        outFile << x << ", ";
        outFile << y << endl;
    }

    outFile.close();

}

void Pool2D::outputPool(const char *fname) {
    double x, y;

    ofstream outFile;
    outFile.open(fname);

    for (int j = 0; j < this->ny; j++) {
        for (int i = 0; i < this->nx; i++) {
            x = simutils::midpoint(this->x[i], this->x[i+1]);
            y = simutils::midpoint(this->y[j], this->y[j+1]);
            outFile << x << ", ";
            outFile << y << ", ";
            outFile << this->phi[this->methodOrd+j][this->methodOrd+i] << endl;
        }
    }

    outFile.close();
}

void Pool2D::outputPoolVelocity(const char *fname) {
    double x, y;
    int mo = this->methodOrd;

    ofstream outFile;
    outFile.open(fname);

    for (int j = 0; j < this->ny; j++) {
        for (int i = 0; i < this->nx; i++) {
            x = simutils::midpoint(this->x[i], this->x[i+1]);
            y = simutils::midpoint(this->y[j], this->y[j+1]);

            outFile << x << ", ";
            outFile << y << ", ";
            outFile << this->poolU[mo+j][mo+i] << ", ";
            outFile << this->poolV[mo+j][mo+i] << endl;
        }
    }
}

void Pool2D::outputStructure(int structNum, const char *fname) {
    solids->at(structNum).outputEdges(fname);
}


void Pool2D::outputAllStructures(const char *baseName) {
    string baseString(baseName);

    for (int i = 0; i < solids->size(); i++) {
        string temp = baseString + to_string(i);

        (solids->at(i)).outputEdges(temp.c_str());
    }
}

void Pool2D::outputStructureNodes(int structNum, const char *fname) {
    solids->at(structNum).outputNodes(fname);
}

void Pool2D::outputAllStructureNodes(const char *baseName) {
    string baseString(baseName);

    for (int i = 0; i < solids->size(); i++) {
        string temp = baseString + to_string(i);

        (solids->at(i)).outputNodes(temp.c_str());
    }
}


void Pool2D::outputStructureVels(int structNum, const char *fname) {
    solids->at(structNum).outputNodeVels(fname);
}

void Pool2D::outputAllStructureVels(const char *baseName) {
    string baseString(baseName);

    for (int i = 0; i < solids->size(); i++) {
        string temp = baseString + to_string(i);

        (solids->at(i)).outputNodeVels(temp.c_str());
    }
}

void Pool2D::outputTracers(const char *fname) {
    double x, y;
    // int mo = this->methodOrd;

    ofstream outFile;
    outFile.open(fname);

    for (int i = 0; i < nStructs; i++) {
        x = tracers[i].x;
        y = tracers[i].y;

        outFile << x << ", ";
        outFile << y << endl;
    }
}

//------------------------------

/* Abstraction helpers for the pool enum array */
// --------------------------------------------


/**
 * Return the object type at the given index.
*/
objects::FSIObject Pool2D::objAtIndex(int xInd, int yInd) {
    return (this->pool)[yInd][xInd];
}

/** 
 * Check if this object type has a structure cell in the given direction
*/
bool Pool2D::hasStructInDir(objects::FSIObject obj, objects::FSIObject dir) {
    return (obj % dir) != 0;
}

/**
 * Check if the object type is an interface
 * 
 * TODO: just check that its not a structure
*/
bool Pool2D::isInterface(objects::FSIObject obj) {
    return obj != objects::FLUID_C && obj != objects::STRUCTURE;
}

/**
 * Check if the interface has only a single structure as the neighbour
 *
*/
bool Pool2D::isNormalInterface(objects::FSIObject obj) {
    return obj == objects::FLUID_N || obj == objects::FLUID_S
                || obj == objects::FLUID_E || objects::FLUID_W;
}

bool Pool2D::enoInRangeX(int val) {
    return simutils::in_range(val, 1, nx-1);
}

bool Pool2D::enoInRangeY(int val) {
    return simutils::in_range(val, 1, ny-1);
}

/**
 * Set of function to indicate if a cell value is usable within an ENO stencil
*/
bool Pool2D::isUsableU(int i, int j) {
    // return pool[j][i] == objects::FLUID_C || pool[j][i+1] == objects::FLUID_C;
    if (enoInRangeX(i) && enoInRangeY(j)) {
        return (pool[j][i] != objects::STRUCTURE || pool[j][i+1] == objects::FLUID_C);
    } else {
        return false;
    }
}

bool Pool2D::isUsableV(int i, int j) {
    // return pool[j][i] == objects::FLUID_C || pool[j+1][i] == objects::FLUID_C;
    if (enoInRangeX(i) && enoInRangeY(j)) {
        return pool[j][i] != objects::STRUCTURE || pool[j+1][i] == objects::FLUID_C;
    } else {
        return false;
    }
}

/** 
 * These functions build the appropriate stencils for the pool 
 * for use with the ENO scheme.
 * 
 * Return the number of stencil points, the valid stencil points
*/

int Pool2D::buildDxStencil(int i, int j, bool sten[7]) {
    int c = 3;

    // Only using second order method
    sten[0] = false;
    sten[6] = false;

    // Go though each of the stencil points and check if they are valid.
    int count = 0;

    // Check the center point
    if (isUsableV(c, j) && isUsableV(c+1, j) && isUsableV(c, j-1)
            && isUsableV(c+1, j-1) && isUsableU(c, j)) {
        sten[c] = true;
        count++;
    }

    // Checking upwards
    for (int l = c+1; l < 7; l++) {
        if (sten[l-1] && isUsableV(l, j) && isUsableV(l+1, j) && isUsableV(l, j-1)
            && isUsableV(l+1, j-1) && isUsableU(l, j)) {
            
            sten[l] = true;
            count++;
        }
    }

    // Checking downwards
    for (int l = c-1; l >= 0; l--) {
        if (sten[l+1] && isUsableV(l, j) && isUsableV(l+1, j) && isUsableV(l, j-1)
            && isUsableV(l+1, j-1) && isUsableU(l, j)) {
            
            sten[l] = true;
            count++;
        }
    }

    return count;
}

int Pool2D::buildDyStencil(int i, int j, bool sten[7]) {
    int c = 3;

    // Only using second order method
    sten[0] = false;
    sten[6] = false;

    // Go though each of the stencil points and check if they are valid.
    int count = 0;

    // Check the center point
    if (isUsableU(i, c) && isUsableU(i, c+1) && isUsableU(i-1, c)
            && isUsableU(i-1, c+1) && isUsableV(i, c)) {
        sten[c] = true;
        count++;
    }

    // Checking upwards
    for (int l = c+1; l < 7; l++) {
        if (sten[l-1] && isUsableU(i, l) && isUsableU(i, l+1) && isUsableU(i-1, l)
            && isUsableV(i-1, l+1) && isUsableV(i, l)) {
            
            sten[l] = true;
            count++;
        }
    }

    // Checking downwards
    for (int l = c-1; l >= 0; l--) {
        if (sten[l+1] && isUsableU(i, l) && isUsableU(i, l+1) && isUsableU(i-1, l)
            && isUsableV(i-1, l+1) && isUsableV(i, l)) {
            
            sten[l] = true;
            count++;
        }
    }
    return count;

}

/**
 * Function to indicates whether a given fluid cell is updatable.
 * TODO: Does not handle boundaries correctly so must be applied only
 *       at internal points
 * 
 * A point is only not updatable in U if its neighbour to the right is an interface.
*/
bool Pool2D::isUpdateableU(int i, int j) {
    return pool[j][i+1] == objects::FLUID_C;
}

/**
 * Point only not updatable in v if its top neighbour is an interface
*/
bool Pool2D::isUpdateableV(int i, int j) { 
    return pool[j+1][i] == objects::FLUID_C;
}

/** 
 * Helper for interpolating general values of phi. Second order accurate. 
*/
double Pool2D::interpolatePhi(double x, double y) {
    int mo = methodOrd;

    // Find the "lim inf" mesh points for the current tracer location. Needed for all modes
    int i = simutils::findLimInfMeshPoint(x, this->x, this->nx+1);
    int j = simutils::findLimInfMeshPoint(y, this->y, this->ny+1);

    // Handle the case where the left-most points are outside of the boundary.
    int il;
    int jl;

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

    double xim = simutils::midpoint(this->x[il-1], this->x[il]);
    double xi = simutils::midpoint(this->x[il], this->x[il+1]);
    double xip = simutils::midpoint(this->x[il+1], this->x[il+2]);

    double yim = simutils::midpoint(this->y[jl-1], this->y[jl]);
    double yi = simutils::midpoint(this->y[jl], this->y[jl+1]);
    double yip = simutils::midpoint(this->y[jl+1], this->y[jl+2]);

    double xMesh[2];
    double yMesh[2];

    int xl, xr, yl, yr;

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

    double vals[4] = {phi[yl][xl], phi[yl][xr],
                            phi[yr][xl], phi[yr][xr]};
    

    // Calcaulte phi with the bilinear interpolation formula
    return simutils::biLinearInterpolation(x, y, xMesh, yMesh, vals);
}

/** 
 * Helper for interpolating general values of phi. Second order accurate.
*/
void Pool2D::interpolatePhiGrad(double x, double y, double phiGrad[2]) {
    int mo = methodOrd;

    // Find the "lim inf" mesh points for the current tracer location. Needed for all modes
    int i = simutils::findLimInfMeshPoint(x, this->x, this->nx+1);
    int j = simutils::findLimInfMeshPoint(y, this->y, this->ny+1);

    // Handle the case where the left-most points are outside of the boundary.
    int il;
    int jl;

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

    double xim = simutils::midpoint(this->x[il-1], this->x[il]);
    double xi = simutils::midpoint(this->x[il], this->x[il+1]);
    double xip = simutils::midpoint(this->x[il+1], this->x[il+2]);

    double yim = simutils::midpoint(this->y[jl-1], this->y[jl]);
    double yi = simutils::midpoint(this->y[jl], this->y[jl+1]);
    double yip = simutils::midpoint(this->y[jl+1], this->y[jl+2]);

    double xMesh[2];
    double yMesh[2];

    int xl, xr, yl, yr;

    if (xim <= x && x <= xi) {
        xMesh[0] = xim;
        xMesh[1] = xi;

        // xl = mo+il;
        // xr = mo+il+1;
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

    double vals[4] = {phi[yl][xl], phi[yl][xr],
                            phi[yr][xl], phi[yr][xr]};

    // Calcaulte phi with the bilinear interpolation formula
    simutils::biLinearGradient(x, y, xMesh, yMesh, vals, phiGrad);
}

double Pool2D::getLevelSetVal(int i, int j) {
    int mo = methodOrd;

    return phi[mo+j][mo+i];
}

/**
 * Get the closest point (xc, yc) on the interface to input point (x, y)
*/
void Pool2D::getClosestInterfacePoint(double P[2], double XC[2]) {
    // Compute the gradient of the level set function using bilinear interpolation
    // formula.
    double phiGrad[2];
    double phiVal;

    // Calcaulte phi
    phiVal = this->interpolatePhi(P[0], P[1]);

    // Calculate phi's gradient
    this->interpolatePhiGrad(P[0], P[1], phiGrad);

    // Normalize the gradient to get the outward normal vector of the signed distance function
    simutils::normalize2D(phiGrad);

    // cout << "phi = " << phiVal << endl;
    // cout << "phiGrad[0] " << phiGrad[0] << " phiGrad[1] = " << phiGrad[1] << endl;
    // cout << "||grad|| " << simutils::eucNorm2D(phiGrad) << endl;

    // Now, use the gradient values to calculate the closest interface point
    XC[0] = P[0] - phiVal*phiGrad[0];
    XC[1] = P[1] - phiVal*phiGrad[1];
}

// --------------------------------------------


/* Getters and setters */
//-----------------------

int Pool2D::getNx() {
    return this->nx;
}
int Pool2D::getNy() {
    return this->ny;
}

double Pool2D::getXMeshVal(int i) {
    return this->x[i];
}

double Pool2D::getYMeshVal(int j) {
    return this->y[j];
}

double Pool2D::getPhiVal(int i, int j) {
    return this->phi[methodOrd+j][methodOrd+i];
}

double Pool2D::getObjU(int i, int j) {
    return this->poolU[methodOrd+j][methodOrd+i];
}

double Pool2D::getObjV(int i, int j) {
    return this->poolV[methodOrd+j][methodOrd+i];
}

//-----------------------
/* Tracer utilities */
// -------------------

/**
 * Getter for the domain membership of this index.
*/
int Pool2D::domainMembership(int i, int j) {
    return this->domainTracker[j][i];
}

/**
 * Print out the member matrix for debugging
*/
void Pool2D::printDomainMatrix() {
    for (int j = 0; j < this->ny; j++) {
        for (int i = 0; i < this->nx; i++) {
            cout << this->domainTracker[j][i] << ", ";
        }
        cout << endl;
    }
    cout << endl;
}

/**
 * Sets up the domain matrix, manipulating the instance vars the correct way.
 * 
 * Uses the domain tracers, so these must be properly updated at all times.
*/
void Pool2D::setUpDomainArray() {
    // Set all of the domain member array to "undescovered"
    simutils::set_constant(ny, nx, DOMAIN_UNDESCOVERED, this->domainTracker);

    // For each immersed structure in the Pool, label subset of the domain they
    // take up.
    for (int l = 0; l < this->nStructs; l++) {
        bfsFromTracer(l);
    }

    // Label all undiscovered domain points as belonging to the fluid domain
    for (int j = 0; j < this->ny; j++) {
        for (int i = 0; i < this->nx; i++) {
            if (this->domainTracker[j][i] == DOMAIN_UNDESCOVERED) {
                this->domainTracker[j][i] = DOMAIN_FLUID;
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
 * TODO: should not allow any kind of array bound violation. In general usage this should be impossible but still.
 * 
 * A Bredth-first search is done about the pool region that the tracer for this object is
 * currently contained in.
 * 
 * The result is a hashmap (unordered_map) with keys (i, j) coordinates in the pool,
 * with values that are tuples ()
 * 
*/
void Pool2D::bfsFromTracer(int structNum) {
    // Find the grid square that the tracer is currently on. This is the root node for the BFS.
    int i = simutils::findLimInfMeshPoint(tracers[structNum].x, this->x, this->nx+1);
    int j = simutils::findLimInfMeshPoint(tracers[structNum].y, this->y, this->ny+1);

    cout << "For stuct num " << structNum << " tracker grid location is x = " << i << " y = " << j << endl;
    int ri, rj, ni, nj;
    objects::FSIObject FL = objects::FLUID_C; // Make the ternary operator more palatable below

    // Create the queues for each of the indices
    queue<int> queueX;
    queue<int> queueY;

    // Add the first node
    queueX.push(i);
    queueY.push(j);

    // Mark this node as a part of the structure. Otherwise, set it to "discovered"
    this->domainTracker[j][i] = (objAtIndex(i, j) != FL) ? structNum : -1;

    // Start the DFS
    while (!queueX.empty()) {
        ri = queueX.front();
        rj = queueY.front();

        queueX.pop();
        queueY.pop();

        // Discover the neighbouring points
        for (nj = rj-1; nj <= rj+1; nj++) {
            for (ni = ri-1; ni <= ri+1; ni++) {
                if (ni == ri && nj == rj) {
                    continue;
                }

                if (this->domainTracker[nj][ni] == DOMAIN_UNDESCOVERED) { // If this node has not been seen
                    if (objAtIndex(ni, nj) != FL) {
                        // This node belongs to the structure
                        this->domainTracker[nj][ni] = structNum;
                        queueX.push(ni);
                        queueY.push(nj);
                    } else {
                        // This node belongs to the fluid region
                        this->domainTracker[nj][ni] = DOMAIN_FLUID;
                    }
                }
            }
        }

    }

    assert(queueX.empty() && queueY.empty());
}

double *Pool2D::getXMesh() {
    return this->x;
}

double *Pool2D::getYMesh() {
    return this->y;
}

// -------------------

Pool2D::~Pool2D() {
    delete kdPointCloud;
    delete kdTree;

    delete[] this->x;
    delete[] this->y;
    delete solids;

    for (int j = 0; j < this->ny; j++) {
        delete[] this->pool[j];
    }
    delete[] this->pool;

    delete[] this->tracers;

    simutils::free_double(this->ny+2*this->methodOrd, this->nx+2*this->methodOrd, this->phi);
    simutils::free_double(this->ny+2*this->methodOrd, this->nx+2*this->methodOrd, this->phiRk1);
    simutils::free_double(this->ny+2*this->methodOrd, this->nx+2*this->methodOrd, this->phiRk2);
    simutils::free_double(this->ny+2*this->methodOrd, this->nx+2*this->methodOrd, this->phiReInit);
    simutils::free_double(this->ny+2*this->methodOrd, this->nx+2*this->methodOrd, this->poolU);
    simutils::free_double(this->ny+2*this->methodOrd, this->nx+2*this->methodOrd, this->poolV);
    simutils::free_int(this->ny, this->nx, this->domainTracker);
    simutils::free_int(this->ny, this->nx, this->fastMarchingState);

    // delete this->changedBoundaries;
}
