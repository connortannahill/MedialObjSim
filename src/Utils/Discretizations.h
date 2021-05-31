#ifndef DISCRETIZATIONS_H
#define DISCRETIZATIONS_H 

#include "SimUtilities.h"

/**
 * A collection of inline discretization functions. The point of these
 * is to ease programming, as many possible discretizations may be employed
 * for each term. Note that the offets must be correct, we do not check the bounds.
*/
namespace discs {
    /**
     * Discretization schemes used in the book. Use strange averageing procedures
     * that we will throw out ASAP
    */

    /**
     * Discretization schemes used in the book. Use strange averaging procedures
     * that we will throw out ASAP
    */
    inline double firstOrder_conv_usqx(int i, int j, double dx, double **u) {
        return ( simutils::square(u[j][i] + u[j][i+1])
            - simutils::square(u[j][i-1] + u[j][i]))/(4.0*dx);
    }

    inline double firstOrder_conv_uvy(int i, int j, double dy, double **u, double **v) {
        return ( (v[j][i] + v[j][i+1]) * (u[j][i] + u[j+1][i]) 
            -  (v[j-1][i] + v[j-1][i+1]) * (u[j-1][i] + u[j][i]) ) / (4.0*dy);
    }

    inline double firstOrder_conv_uvx(int i, int j, double dx, double **u, double **v) {
        return ( (u[j][i] + u[j+1][i]) * (v[j][i] + v[j][i+1])
            - (u[j][i-1] + u[j+1][i-1]) * (v[j][i-1] + v[j][i]) ) / (4.0*dx);
    }

    inline double firstOrder_conv_vsqy(int i, int j, double dy, double **v) {
        return ( simutils::square(v[j][i] + v[j+1][i]) 
            - simutils::square(v[j-1][i] + v[j][i])) / (4.0*dy);
    }

    inline double firstOrder_lap_uxx(int i, int j, double dx, double **u) {
        return (u[j][i+1] - 2.0*u[j][i] + u[j][i-1]) / (simutils::square(dx));
    }

    inline double firstOrder_lap_uyy(int i, int j, double dy, double **u) {
        return (u[j+1][i] - 2.0*u[j][i] + u[j-1][i]) / (simutils::square(dy));
    }

    inline double firstOrder_lap_vxx(int i, int j, double dx, double **v) {
        return (v[j][i+1] - 2.0*v[j][i] + v[j][i-1]) / (simutils::square(dx));
    }
    inline double firstOrder_lap_vyy(int i, int j, double dy, double **v) {
        return (v[j+1][i] - 2.0*v[j][i] + v[j-1][i]) / (simutils::square(dy));
    }

    /**
     * Third order discretizations
    */
    // First equation
    inline double firstOrder3D_conv_usqx(int i, int j, int k, double dx, double ***u) {
        return ( simutils::square(u[k][j][i] + u[k][j][i+1])
            - simutils::square(u[k][j][i-1] + u[k][j][i]))/(4.0*dx);
    }

    inline double firstOrder3D_conv_uvy(int i, int j, int k, double dy, double ***u, double ***v) {
        return ( (v[k][j][i] + v[k][j][i+1]) * (u[k][j][i] + u[k][j+1][i]) 
            -  (v[k][j-1][i] + v[k][j-1][i+1]) * (u[k][j-1][i] + u[k][j][i]) ) / (4.0*dy);
    }

    inline double firstOrder3D_conv_uwz(int i, int j, int k, double dz, double ***u, double ***w) {
        return ( (w[k][j][i] + w[k][j][i+1]) * (u[k][j][i] + u[k+1][j][i]) 
            -  (w[k-1][j][i] + w[k-1][j][i+1]) * (u[k-1][j][i] + u[k][j][i]) ) / (4.0*dz);
    }

    inline double firstOrder3D_lap_uxx(int i, int j, int k, double dx, double ***u) {
        return (u[k][j][i+1] - 2.0*u[k][j][i] + u[k][j][i-1]) / (simutils::square(dx));
    }

    inline double firstOrder3D_lap_uyy(int i, int j, int k, double dy, double ***u) {
        return (u[k][j+1][i] - 2.0*u[k][j][i] + u[k][j-1][i]) / (simutils::square(dy));
    }
    
    inline double firstOrder3D_lap_uzz(int i, int j, int k, double dz, double ***u) {
        return (u[k+1][j][i] - 2.0*u[k][j][i] + u[k-1][j][i]) / (simutils::square(dz));
    }

    // Second equation
    inline double firstOrder3D_conv_vux(int i, int j, int k, double dx, double ***v, double ***u) {
        return ( (u[k][j][i] + u[k][j+1][i]) * (v[k][j][i] + v[k][j][i+1])
            - (u[k][j][i-1] + u[k][j+1][i-1]) * (v[k][j][i-1] + v[k][j][i]) ) / (4.0*dx);
    }

    inline double firstOrder3D_conv_vsqy(int i, int j, int k, double dy, double ***v) {
        return ( simutils::square(v[k][j][i] + v[k][j+1][i]) 
            - simutils::square(v[k][j-1][i] + v[k][j][i])) / (4.0*dy);
    }

    inline double firstOrder3D_conv_vwz(int i, int j, int k, double dz, double ***v, double ***w) {
        return ( (w[k][j][i] + w[k][j+1][i]) * (v[k][j][i] + v[k+1][j][i])
            - (w[k-1][j][i] + w[k-1][j+1][i]) * (v[k-1][j][i] + v[k][j][i]) ) / (4.0*dz);
    }

    inline double firstOrder3D_lap_vxx(int i, int j, int k, double dx, double ***v) {
        return (v[k][j][i+1] - 2.0*v[k][j][i] + v[k][j][i-1]) / (simutils::square(dx));
    }

    inline double firstOrder3D_lap_vyy(int i, int j, int k, double dy, double ***v) {
        return (v[k][j+1][i] - 2.0*v[k][j][i] + v[k][j-1][i]) / (simutils::square(dy));
    }

    inline double firstOrder3D_lap_vzz(int i, int j, int k, double dz, double ***v) {
        return (v[k+1][j][i] - 2.0*v[k][j][i] + v[k-1][j][i]) / (simutils::square(dz));
    }

    // Third equation
    inline double firstOrder3D_conv_wux(int i, int j, int k, double dx, double ***w, double ***u) {
        return ( (u[k][j][i] + u[k+1][j][i]) * (w[k][j][i] + w[k][j][i+1])
            - (u[k][j][i-1] + u[k+1][j][i-1]) * (w[k][j][i-1] + w[k][j][i]) )  / (4.0*dx);
    }

    inline double firstOrder3D_conv_wvy(int i, int j, int k, double dy, double ***w, double ***v) {
        return ( (v[k][j][i] + v[k+1][j][i]) * (w[k][j][i] + w[k][j+1][i])
            - (v[k][j-1][i] + v[k+1][j-1][i]) * (w[k][j-1][i] + w[k][j][i] ) )  / (4.0*dy);
    }
    
    inline double firstOrder3D_conv_wsqz(int i, int j, int k, double dz, double ***w) {
        return ( simutils::square(w[k][j][i] + w[k+1][j][i]) 
            - simutils::square(w[k-1][j][i] + w[k][j][i])) / (4.0*dz);
    }

    inline double firstOrder3D_lap_wxx(int i, int j, int k, double dx, double ***w) {
        return (w[k][j][i+1] - 2.0*w[k][j][i] + w[k][j][i-1]) / (simutils::square(dx));
    }

    inline double firstOrder3D_lap_wyy(int i, int j, int k, double dy, double ***w) {
        return (w[k][j+1][i] - 2.0*w[k][j][i] + w[k][j-1][i]) / (simutils::square(dy));
    }

    inline double firstOrder3D_lap_wzz(int i, int j, int k, double dz, double ***w) {
        return (w[k+1][j][i] - 2.0*w[k][j][i] + w[k-1][j][i]) / (simutils::square(dz));
    }

    /**
     * Third order ENO scheme with stencil
     *     x_{i-3} < x_{i-2} < x_{i-1} < x_{i} < x_{i+1} < x_{i+2} < x_{i+3}
     * Evaluated at x_i
     * 
     * prefStencil is an array of boolean values for each stencil point. If index j
     * is 1, then it can be used within the ENO scheme. If index j is set to 0, it can not
     * be used. Note index i and 3 points on either side must be set to 1.
     * 
     * NOTE: assumes a uniform mesh
     * 
     * TODO: must be a more ellegant and efficient implementation
     * 
     * Modify so that if stencil does not have enough points, we use the less accurate method
    */
    inline double thirdOrdENO(double evalPnt, double *sP, double *sY, int upDir, bool *prefStencil) {
        // Current leftmost and rightmost stencil points when building up the scheme
        // int l = 3;
        // int r = 3;
        int i = 3;
        int k = i; // k tracks the left-most point in the stencil

        // The value fo the interpolant
        double n = 0;

        int count = 0;
        for (int l = 0; l < 7; l++) 
            count += prefStencil[l];
        
        if (count < 2) {
            assert(false);
        }

        // Mesh spacing
        // double h = sP[1] - sP[0];

        // Undivided differences for comparisons
        double ul, ur, c;

        // First term chosen based on the upwind direction if the information is available
        // ul = (sY[k] - sY[k-1])/h;
        // ur = (sY[k+1] - sY[k])/h;
        ul = (sY[k] - sY[k-1])/(sP[k] - sP[k-1]);
        ur = (sY[k+1] - sY[k])/(sP[k+1] - sP[k]);

        if (prefStencil[k-1] && prefStencil[k+1]) {
            if (upDir == -1) {
                k--;
                n += ul;
            } else {
                n += ur;
            }
        } else {
            if (prefStencil[k-1]) {
                k--;
                n += ul;
            } else {
                n += ur;
            }
        }

        if (count < 3) {
            return n;
        }

        // Second term comparison
        // ul = sY[k+1] - 2.0*sY[k] + sY[k-1];
        // ur = sY[k+2] - 2.0*sY[k+1] + sY[k];
        ul = ( (sY[k+1] - sY[k])*(sP[k] - sP[k-1]) - (sY[k] - sY[k-1])*(sP[k+1] - sP[k]) ) 
            / ( (sP[k+1] - sP[k-1])*(sP[k+1]-sP[k])*(sP[k]-sP[k-1]) );
        ur = ( (sY[k+2] - sY[k+1])*(sP[k+1] - sP[k]) - (sY[k+1] - sY[k])*(sP[k+2] - sP[k+1]) ) 
            / ( (sP[k+2] - sP[k])*(sP[k+2]-sP[k+1])*(sP[k+1]-sP[k]) );

        if (prefStencil[k-1] && prefStencil[k+2]) {
            if (abs(ul) < abs(ur)) {
                c = ul; /// (2.0*simutils::square(h));
                n += c*((evalPnt - sP[k]) + (evalPnt - sP[k+1]));
                k--;
            } else {
                c = ur;// / (2.0*simutils::square(h));
                n += c*((evalPnt - sP[k]) + (evalPnt - sP[k+1]));
            }
        } else {
            if (prefStencil[k-1]) {
                c = ul;// / (2.0*simutils::square(h));
                n += c*((evalPnt - sP[k]) + (evalPnt - sP[k+1]));
                k--;
            } else {
                c = ur;// / (2.0*simutils::square(h));
                n += c*((evalPnt - sP[k]) + (evalPnt - sP[k+1]));
            }
        }

        if (count < 4) {
            return n;
        }

        // Third term comparison
        // ul = sY[k+2] - 3.0*sY[k+1] + 3.0*sY[k] - sY[k-1];
        // ur = sY[k+3] - 3.0*sY[k+2] + 3.0*sY[k+1] - sY[k];
        ul = ( ( (sY[k+2] - sY[k+1])*(sP[k+1] - sP[k]) - (sY[k+1] - sY[k])*(sP[k+2] - sP[k+1]) ) / ( (sP[k+2] - sP[k])*(sP[k+2]-sP[k+1])*(sP[k+1]-sP[k]) )
            - ( (sY[k+1] - sY[k])*(sP[k] - sP[k-1]) - (sY[k] - sY[k-1])*(sP[k+1] - sP[k]) ) / ( (sP[k+1] - sP[k-1])*(sP[k+1]-sP[k])*(sP[k]-sP[k-1]) ) )
            / (sP[k+2] - sP[k-1]);
        ur = ( ( (sY[k+3] - sY[k+2])*(sP[k+2] - sP[k+1]) - (sY[k+2] - sY[k+1])*(sP[k+3] - sP[k+2]) ) / ( (sP[k+3] - sP[k+1])*(sP[k+3]-sP[k+2])*(sP[k+2]-sP[k+1]) )
            - ( (sY[k+2] - sY[k+1])*(sP[k+1] - sP[k]) - (sY[k+1] - sY[k])*(sP[k+2] - sP[k+1]) ) / ( (sP[k+2] - sP[k])*(sP[k+2]-sP[k+1])*(sP[k+1]-sP[k]) ) )
            / (sP[k+3] - sP[k]);
        
        if (prefStencil[k-1] && prefStencil[k+3]) {
            if (abs(ul) < abs(ur)) {
                c = ul;// / (6.0*simutils::cube(h));
                n += c*((evalPnt-sP[k+1])*(evalPnt-sP[k+2]) 
                    + (evalPnt-sP[k])*(evalPnt-sP[k+2])
                    + (evalPnt-sP[k])*(evalPnt-sP[k+1]));
                k--;
            } else {
                c = ur;// / (6.0*simutils::cube(h));
                n += c*((evalPnt-sP[k+1])*(evalPnt-sP[k+2]) 
                    + (evalPnt-sP[k])*(evalPnt-sP[k+2])
                    + (evalPnt-sP[k])*(evalPnt-sP[k+1]));
            }
        } else {
            if (prefStencil[k-1]) {
                c = ul;// / (6.0*simutils::cube(h));
                n += c*((evalPnt-sP[k+1])*(evalPnt-sP[k+2]) 
                    + (evalPnt-sP[k])*(evalPnt-sP[k+2])
                    + (evalPnt-sP[k])*(evalPnt-sP[k+1]));
                k--;
            } else {
                c = ur;// / (6.0*simutils::cube(h));
                n += c*((evalPnt-sP[k+1])*(evalPnt-sP[k+2]) 
                    + (evalPnt-sP[k])*(evalPnt-sP[k+2])
                    + (evalPnt-sP[k])*(evalPnt-sP[k+1]));
            }

        }

        return n;
    }

    /**
     * Fifth order WENO scheme.
    */
   inline double fifthOrdWENO(double *sPoints, double *sY) {
       return 0;
   }

    inline double varOrdENO(int methodOrd, double *sPoints, double *sY) {
        std::cout << "Hi, you shouldn't be using variable order ENO yet." << std::endl;
        return 0;
    }
}

#endif