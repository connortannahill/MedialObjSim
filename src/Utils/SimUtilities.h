#ifndef SIM_UTILITIES_H
#define SIM_UTILITIES_H

#include <iostream>
#include <stdio.h>
#include <math.h>

/**
 * TODO: refactor this namespace so that not everything is inline. Check out why there is a compiler error
 *       in some cases.
*/

namespace simstructs {
    struct tracer2D {
        double x;
        double y;
        double mass;
        double volume;
        double density;
        double fNetX;
        double fNetY;
        double u; // x component of velocity
        double v; // y component of velocity
        bool isDeformable;
    };

    struct tracer3D {
        double x;
        double y;
        double z;
        double mass;
        double volume;
        double density;
        double fNetX;
        double fNetY;
        double fNetZ;
        double u; // x component of velocity
        double v; // y component of velocity
        double w; // w component of velocity
        bool isDeformable;
    };

    struct surfaceDiff {
        int pointId;
        double val;
    };
}

namespace simutils {
    /**
     * Set all values in array of size n to constant value.
     */
    inline void set_constant(int n, double val, double *x) {
        int i;
        for (i = 0; i < n; i++) {
            x[i] = val;
        }
    }

    /**
     * Set all values in array of size n to constant value.
     */
    inline void set_constant(int n, int val, int *x) {
        int i;
        for (i = 0; i < n; i++) {
            x[i] = val;
        }
    }

    inline double dmax(double a, double b) {
        return (a > b) ? a : b;
    }

    inline double dmin(double a, double b) {
        return (a < b) ? a : b;
    }

    /**
     * Epsilon equality
    */
    inline bool eps_equal(double a, double b, double eps) {
        return abs(a - b) < eps;
    }

    /**
     * Set all values in array of size n x m to constant value
     */
    inline void set_constant(int n, int m, double val, double **x) {
        int i, j;
        for (i = 0; i < n; i++) {
            for (j = 0; j < m; j++) {
                x[i][j] = val;
            }
        }
    }

    /**
     * Compute the 3D cross product
    */
   inline void cross_product_3D(double u[3], double v[3], double uCrossV[3]) {
       uCrossV[0] = u[1]*v[2] - u[2]*v[1];
       uCrossV[1] = u[2]*v[0] - u[0]*v[2];
       uCrossV[2] = u[0]*v[1] - u[1]*v[0];
   }

    /**
     * Set all values in array of size n x m to constant value
     */
    inline void set_constant(int n, int m, int val, int **x) {
        int i, j;
        for (i = 0; i < n; i++) {
            for (j = 0; j < m; j++) {
                x[i][j] = val;
            }
        }
    }

    /**
     * Set all values in array of size n x m x l to constant value
     */
    inline void set_constant(int n, int m, int l, double val, double ***x) {
        int i, j, k;
        for (i = 0; i < n; i++) {
            for (j = 0; j < m; j++) {
                for (k = 0; k < l; k++) {
                    x[i][j][k] = val;
                }
            }
        }
    }

    /**
     * Set all values in array of size n x m x l to constant value
     */
    inline void set_constant(int n, int m, int l, int val, int ***x) {
        int i, j, k;
        for (i = 0; i < n; i++) {
            for (j = 0; j < m; j++) {
                for (k = 0; k < l; k++) {
                    x[i][j][k] = val;
                }
            }
        }
    }

    /**
     * Create new array of size n and set to constant value.
     */
    inline double *new_constant(int n, double val) {
        double *x = new double[n]; 

        int i;
        for (i = 0; i < n; i++) {
            x[i] = val;
        }

        return x;
    }

    /**
     * Create new array of size n and set to constant value.
     */
    inline int *new_constant(int n, int val) {
        int *x = new int[n]; 

        int i;
        for (i = 0; i < n; i++) {
            x[i] = val;
        }

        return x;
    }

    // Given three colinear points p, q, r, the function checks if 
    // point q lies on line segment 'pr' 
    inline bool onSegment(double p[2], double q[2], double r[2]) 
    { 
        if (q[0] <= dmax(p[0], r[0]) && q[0] >= dmin(p[0], r[0]) && 
            q[1] <= dmax(p[1], r[1]) && q[1] >= dmin(p[1], r[1])) 
        return true; 
    
        return false; 
    } 
    
    // To find orientation of ordered triplet (p, q, r). 
    // The function returns following values 
    // 0 --> p, q and r are colinear 
    // 1 --> Clockwise 
    // 2 --> Counterclockwise 
    inline int orientation(double p[2], double q[2], double r[2]) 
    { 
        // See https://www.geeksforgeeks.org/orientation-3-ordered-points/ 
        // for details of below formula. 
        int val = (q[1] - p[1]) * (r[0] - q[0]) - 
                (q[0] - p[0]) * (r[1] - q[1]); 
    
        if (val == 0) return 0;  // colinear 
    
        return (val > 0) ? 1: 2; // clock or counterclock wise 
    } 
    
    // The main function that returns true if line segment 'p1q1' 
    // and 'p2q2' intersect. 
    inline bool line_intersect(double p1[2], double q1[2], double p2[2], double q2[2]) 
    { 
        // Find the four orientations needed for general and 
        // special cases 
        int o1 = orientation(p1, q1, p2); 
        int o2 = orientation(p1, q1, q2); 
        int o3 = orientation(p2, q2, p1); 
        int o4 = orientation(p2, q2, q1); 
    
        // General case 
        if (o1 != o2 && o3 != o4) 
            return true; 
    
        // Special Cases 
        // p1, q1 and p2 are colinear and p2 lies on segment p1q1 
        if (o1 == 0 && onSegment(p1, p2, q1)) return true; 
    
        // p1, q1 and q2 are colinear and q2 lies on segment p1q1 
        if (o2 == 0 && onSegment(p1, q2, q1)) return true; 
    
        // p2, q2 and p1 are colinear and p1 lies on segment p2q2 
        if (o3 == 0 && onSegment(p2, p1, q2)) return true; 
    
        // p2, q2 and q1 are colinear and q1 lies on segment p2q2 
        if (o4 == 0 && onSegment(p2, q1, q2)) return true; 
    
        return false; // Doesn't fall in any of the above cases 
    } 

    /**
     * Create new array of size n x m and set to constant value
     */
    inline double **new_constant(int n, int m, double val) {
        int i, j;
        double **x = new double*[n];
        for (i = 0; i < n; i++)
            x[i] = new double[m];

        for (i = 0; i < n; i++) {
            for (j = 0; j < m; j++) {
                x[i][j] = val;
            }
        }

        return x;
    }

    /**
     * Create new array of size n x m and set to constant value
     */
    inline int **new_constant(int n, int m, int val) {
        int i, j;
        int **x = new int*[n];
        for (i = 0; i < n; i++)
            x[i] = new int[m];

        for (i = 0; i < n; i++) {
            for (j = 0; j < m; j++) {
                x[i][j] = val;
            }
        }

        return x;
    }

    /**
     * Create new array of size n x m x l and set to constant value
     */
    inline double ***new_constant(int n, int m, int l, double val) {
        int i, j, k;

        double ***x = new double**[n];

        for (i = 0; i < n; i++) {
            x[i] = new double*[m];
            for (j = 0; j < m; j++) {
                x[i][j] = new double[l];
            }
        }

        for (i = 0; i < n; i++) {
            for (j = 0; j < m; j++) {
                for (k = 0; k < l; k++) {
                    x[i][j][k] = val;
                }
            }
        }

        return x;
    }

    /**
     * Create new array of size n x m x l and set to constant value
     */
    inline int ***new_constant(int n, int m, int l, int val) {
        int i, j, k;

        int ***x = new int**[n];

        for (i = 0; i < n; i++) {
            x[i] = new int*[m];
            for (j = 0; j < m; j++) {
                x[i][j] = new int[l];
            }
        }

        for (i = 0; i < n; i++) {
            for (j = 0; j < m; j++) {
                for (k = 0; k < l; k++) {
                    x[i][j][k] = val;
                }
            }
        }

        return x;
    }

    /**
     * Easily free 2D array
    */
   inline void free_double(int n, int m, double **x) {
       for (int i = 0; i < n; i++) {
           delete[] x[i];
       }
       delete[] x;
   }

    /**
     * Easily free 2D array
    */
   inline void free_int(int n, int m, int **x) {
       for (int i = 0; i < n; i++) {
           delete[] x[i];
       }
       delete[] x;
   }
   
    /**
     * Easily free 3D array
    */
   inline void free_double(int n, int m, int l, double ***x) {
       int i, j;
       for (i = 0; i < n; i++) {
           for (j = 0; j < m; j++) {
               delete[] x[i][j];
           }
           delete[] x[i];
       }
       delete[] x;
   }

    /**
     * Easily free 3D array
    */
   inline void free_int(int n, int m, int l, int ***x) {
       int i, j;
       for (i = 0; i < n; i++) {
           for (j = 0; j < m; j++) {
               delete[] x[i][j];
           }
           delete[] x[i];
       }
       delete[] x;
   }

   inline double mean(int n, double *x) {
       double sum = 0;
       for (int i = 0; i < n; i++) {
           sum += x[i];
       }
       return sum /( (double) n);
   }

   /**
    * Take the mean of a 2D matrix
   */
    inline double mean(int n, int m, double **x) {
        double sum = 0;
        double count = (double) n * m;

        for (int j = 0; j < n; j++) {
            for (int i = 0; i < m; i++) {
                sum += x[j][i];
            }
        }

        return sum / count;
    }

    /**
     * Uniformly partition [xa, xb] into ns subintervals.
     * Result is stored in x.
    */
    inline void linspace(double xa, double xb, int ns, double *x){
        for (int i = 0; i < ns+1; i++) {
            x[i] = xa + ((double)i)*(xb - xa)/ns;
        }
    }

    /**
     * Get the maximum value of n x m 2D array.
    */
   inline double max(int n, int m, double **u) {
       double max = u[0][0];
       int i, j;
       for (i = 0; i < n; i++) {
           for (j = 0; j < m; j++) {
               max = (u[i][j] > max) ? u[i][j] : max;
           }
       }

       return max;
   }

   /**
    * Compute the dot product of two 2D vecs
   */
  inline double ddot2d(double a[2], double b[2]) {
      return a[0]*b[0] + a[1]*b[1];
  }

   /**
    * Compute the dot product of two 3D vecs
   */
  inline double ddot3d(double a[3], double b[3]) {
      return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  }

  /**
   * Check if value is in range.
  */
  inline bool in_range(double val, double low, double high) {
      return (val >= low) && (val <= high);
  }

  /**
   * Check if value is in range.
  */
  inline bool int_in_range(int val, int low, int high) {
      return (val >= low) && (val <= high);
  }

   /**
    * Compute the 1 norm of a large vector
   */
   inline double l1(double *q, int len) {
       double acc = 0.0;
       for (int i = 0; i < len; i++) {
           acc += abs(q[i]);
       }

       return acc;
   }

   /**
    * Compute the 2 norm of a large vector
   */
   inline double l2(double *q, int len) {
       double acc = 0.0;
       for (int i = 0; i < len; i++) {
           acc += q[i]*q[i];
       }

       return sqrt(acc);
   }

   /**
    * Get the maximum value of an n x m x l array
   */
   inline double max(int n, int m, int l, double ***u) {
        double max = u[0][0][0];
        for (int k = 0; k < n; k++) {
            for (int j = 0; j < m; j++) {
                for (int i = 0; i < l+1; i++) {
                    max = (u[k][j][i] > max) ? u[k][j][i] : max;
                }
            }
        }

       return max;
   }

    /**
    * Print u.
    */
    inline void printU(int nx, int ny, int nGhost, double **u) {
        for (int j = 0; j < ny+2; j++) {
            for (int i = 0; i < nx+1; i++) {
                // std::cout << u[j][i] << ", ";
                printf("%.6f, ", u[j][i]);
            }
            printf("\n");
            // std::cout << std::endl;
        }
        printf("\n");
        // std::cout << std::endl;
    }

    inline void printU(int nx, int ny, int nz, int nGhost, double ***u) {
        for (int k = 0; k < nz+2; k++) {
            for (int j = 0; j < ny+2; j++) {
                for (int i = 0; i < nx+1; i++) {
                    printf("%.6f, ", u[k][j][i]);
                }
                printf("\n");
            }
            printf("\n\n");
        }
    }

    /**
    * Print v.
    */
    inline void printV(int nx, int ny, int nGhost, double **v) {
        for (int j = 0; j < ny+1; j++) {
            for (int i = 0; i < nx+2; i++) {
                // std::cout << v[j][i] << ", ";
                printf("%.6f, ", v[j][i]);
            }
            printf("\n");
            // std::cout << std::endl;
        }
        printf("\n");
        // std::cout << std::endl;
    }

    inline void printMat(int nx, int ny, double **m) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                std::cout << m[j][i] << ", ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    /**
     * Construct the Barycentric weights on the interval [0, 1]
     * for the uniform mesh [0, h, 2h, ..., nh], where h = 1/n
     * 
     * Args:
     * - n is the degree of the interpolating polynomial.
     * - bw is the Barycentric weights, array of length n+1.
    */
   inline void barycentricInterp(int n, double bw[]) {
       // Calculate the mesh spacing
       double h = 1.0/((double) n);

       // Compute the Barycentric weights
       double temp;
       for (int j = 0; j < n+1; j++) {
           temp = 1.0;
           for (int k = 0; k < n+1; k++) {
               if (k == j) {
                   continue;
               } else {
                   temp *= h * (double)(j - k);
               }
           }

           bw[j] = 1/temp;
       }
   }

   /**
    * Evaluates the Lagrange interpolant lb at points x. 
    * Number of points at which the solution will be evaluated
    * Let n be the degree of the desired polynomial interpolant.
    * 
    * Args:
    * - x: the point at which the interpolant will be evaluated.
    * - xa: the right-most interpolation point
    * - xb: the left-most interpolation point
    * - intWidth: absolute size of the interval on which the function being interpolated is defined.
    * - n: the degree of the interpolating polynomial
    * - bw: the Barycentric weights for a uniform mesh on the interval [0, 1].
    * - y: the n+1 points being interpolated.
    */
   inline double barycentricInterv(double x, double xa, double xb, int n, double bw[], double y[]) {
       double h = 1.0/((double)n);

       // Theta value: x mapped onto [0, 1]
       double theta = (x - xa) / (xb - xa);

    //    intWidth = xb - xa;

       double p = 0;
       double l = 1;
       double pdiff;
       for (int i = 0; i < n+1; i++) {
           pdiff = theta - ((double)i)*h;
           p += (bw[i]/pdiff)*y[i];
           l *= pdiff;
       }

    //    return (l/intWidth)*p;
       return l*p;
   }

   /**
    * Evaluates the derivatives of the second Barycentric form at an interpolation
    * point xi. Let n be the degree of the desired polynomial interpolant (note its
    * degree will be of one degree less and one order of accuracy less).
    * 
    * Args:
    * - xi: the interpolation point at which we want to evaluate the derivative.
    * - xa: the right-most interpolation point
    * - xb: the left-most interpolation point
    * - n: the degree of the interpolating polynomial
    * - bw: the Barycentric weights for a uniform mesh on the interval [0, 1].
    * - y: the n+1 points being interpolated.
    * Note: untested
   */
  inline double barycentricDeriv(int xi, double xa, double xb, int n, double bw[], double y[]) {
      double h = 1.0/((double)n);

      double t;
      double p = 0;
      double li = 0;
      for (int j = 0; j < n+1; j++) {
          if (j == xi)
              continue;

          t = (bw[j]/bw[xi]) / (h*(xi - j));

          p += t*y[j];
          li += t;
      }

      p += -li*y[xi];
      p /= (xb - xa);

      return p;
    }

    inline double square(double u) {
        return u * u;
    }

    inline double cube(double u) {
        return u * u * u;
    }

    inline int sign(double x) {
        return (x > 0) ? 1 : -1;
    }

    inline void copyVals(int nx, double *src, double *tar) {
        for (int i = 0; i < nx; i++) {
            tar[i] = src[i];
        }
    }

    inline void copyVals(int nx, int ny, double **src, double **tar) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                tar[j][i] = src[j][i];
            }
        }
    }

    inline void copyVals(int nx, int ny, int **src, int **tar) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                tar[j][i] = src[j][i];
            }
        }
    }

    inline void copyVals(int nx, int ny, int nz, double ***src, double ***tar) {
        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    tar[k][j][i] = src[k][j][i];
                }
            }
        }
    }

    inline void copyVals(int nx, int ny, int nz, int ***src, int ***tar) {
        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    tar[k][j][i] = src[k][j][i];
                }
            }
        }
    }

    inline double eucNorm2D(double n[2]) {
        return sqrt(simutils::square(n[0]) + simutils::square(n[1]));
    }

    inline double eucDiff2D(double a[2], double b[2]) {
        return sqrt(simutils::square(a[0] - b[0]) + simutils::square(a[1] - b[1]));
    }

    inline double eucNorm3D(double n[3]) {
        return sqrt(simutils::square(n[0]) + simutils::square(n[1])
                    + simutils::square(n[2]));
    }

    inline double eucDiff3D(double a[3], double b[3]) {
        return sqrt(simutils::square(a[0] - b[0]) + simutils::square(a[1] - b[1]) + simutils::square(a[2] - b[2]));
    }

    inline void ghostCopy(int nx, int ny, int nGhost, double **arr, double **copy) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                copy[j+nGhost][i+nGhost] = arr[j+nGhost][i+nGhost];
            }
        }
    }

    inline int sum(int n, int *arr) {
        int sum = 0;
        for (int i = 0; i < n; i++) {
            sum += arr[i];
        }
        return sum;
    }

    inline double midpoint(double xa, double xb) {
        return xa + (xb - xa)/2.0;
    }

    inline int imax(int a, int b) {
        return (a > b) ? a : b;
    }

    inline int imin(int a, int b) {
        return (a < b) ? a : b;
    }


    /**
     * Bi-linear interpolation on a tensor product grid
     * 
     * The vals array is set up such like:
     * 
     * (vals[6]=c011)------------(vals[7]=c111)
     * |                            |
     * |                            |
     * |                            |
     * |                            |
     * |                            |
     * (vals[4]=c001)------------(vals[5]=c101)
     *                ^ 
     *                |
     *                |
     *               "up"
     * 
     * (vals[2]=c010)------------(vals[3]=c110)
     * |                            |
     * |                            |
     * |                            |
     * |                            |
     * |                            |
     * (vals[0]=c000)------------(vals[1]=c100)
     * 
     * 
    */
    inline double triLinearInterpolation(double x, double y, double z, double xMesh[2],
                                        double yMesh[2], double zMesh[2], double f[8]) {
        double xd = (x - xMesh[0])/(xMesh[1] - xMesh[0]);
        double yd = (y - yMesh[0])/(yMesh[1] - yMesh[0]);
        double zd = (z - zMesh[0])/(zMesh[1] - zMesh[0]);

        double c00 = f[0]*(1 - xd) + f[1]*xd;
        double c01 = f[4]*(1 - xd) + f[5]*xd;
        double c10 = f[2]*(1 - xd) + f[3]*xd;
        double c11 = f[6]*(1 - xd) + f[7]*xd;

        double c0 = c00*(1 - yd) + c10*yd;
        double c1 = c01*(1 - yd) + c11*yd;

        // double c0 = (f[0]*(1 - xd) + f[1]*xd)*(1 - yd) + (f[2]*(1 - xd) + f[3]*xd)*yd;
        // double c1 = (f[4]*(1 - xd) + f[5]*xd)*(1 - yd) + (f[6]*(1 - xd) + f[7]*xd)*yd;

        return c0*(1 - zd) + c1*zd;
        // return f[0]*(1 - xd)*(1 - yd)*(1 - zd) + f[1]*xd*(1 - yd)*(1 - zd) + f[2]*(1 - xd)*yd*(1 - zd) + f[3]*xd*yd*(1 - zd) 
        // + f[4]*(1 - xd)*(1 - yd)*zd + f[5]*xd*(1 - yd)*zd + f[6]*(1 - xd)*yd*zd + f[7]*xd*yd*zd;
    }

    inline void triLinearGradient(double x, double y, double z, double xMesh[2],
                                        double yMesh[2], double zMesh[2], double f[8],
                                        double phiGrad[3]) {
        double xd = (x - xMesh[0])/(xMesh[1] - xMesh[0]);
        double yd = (y - yMesh[0])/(yMesh[1] - yMesh[0]);
        double zd = (z - zMesh[0])/(zMesh[1] - zMesh[0]);

        double xdp = 1.0/(xMesh[1] - xMesh[0]);
        double ydp = 1.0/(yMesh[1] - yMesh[0]);
        double zdp = 1.0/(zMesh[1] - zMesh[0]);

        double c00 = f[0]*(1 - xd) + f[1]*xd;
        double c01 = f[4]*(1 - xd) + f[5]*xd;
        double c10 = f[2]*(1 - xd) + f[3]*xd;
        double c11 = f[6]*(1 - xd) + f[7]*xd;

        double c0 = c00*(1 - yd) + c10*yd;
        double c1 = c01*(1 - yd) + c11*yd;

        phiGrad[0] = ( (-f[0]*xdp + f[1]*xdp)*(1.0-yd) + (-f[2]*xdp + f[3]*xdp)*yd )*(1 - zd)
            + ( (-f[4]*xdp + f[5]*xdp)*(1-yd) + (-f[6]*xdp + f[7]*xdp)*yd )*zd;
        phiGrad[1] = (c00*(-ydp) + c10*ydp)*(1.0 - zd) + (c01*(-ydp) + c11*ydp)*zd;
        phiGrad[2] = c0*(-zdp) + c1*(zdp);


    }

    /**
     * Bi-linear interpolation on a tensor product grid
     * 
     * The vals array is set up such like:
     * 
     * (vals[2])------------(vals[3])
     * |                            |
     * |                            |
     * |                            |
     * |                            |
     * |                            |
     * (vals[0])------------(vals[1])
     * 
     * 
    */
    inline double biLinearInterpolation(double x, double y, double xMesh[2],
                                        double yMesh[2], double f[4]) {
        return (1/((xMesh[1]-xMesh[0])*(yMesh[1] - yMesh[0])))*(
                    f[0]*(xMesh[1] - x)*(yMesh[1] - y)
                        + f[1]*(x - xMesh[0])*(yMesh[1] - y)
                        + f[2]*(xMesh[1] - x)*(y - yMesh[0])
                        + f[3]*(x - xMesh[0])*(y - yMesh[0]) );
    }

    /**
     * Use bi-linear interpolation formule to get first order gradient approximation
    */
    inline void biLinearGradient(double x, double y, double xMesh[2],
                                    double yMesh[2], double f[4], double grad[2]) {
        double hx = xMesh[1] - xMesh[0];
        double hy = yMesh[1] - yMesh[0];

        grad[0] = (1/(hx*hy))*(-f[0]*(yMesh[1] - y) + f[1]*(yMesh[1] - y)
                        - f[2]*(y - yMesh[0]) + f[3]*(y - yMesh[0]));
        grad[1] = (1/(hx*hy))*(-f[0]*(xMesh[1] - x) - f[1]*(x - xMesh[0])
                        + f[2]*(xMesh[1] - x) + f[3]*(x - xMesh[0]));
    }

    /**
     * Perform a bisection-esque search on a 1D array w_mesh to with indices [0, ..., nw-1] for the maximum mesh
     * point w_i for which w_i < w. Return the i that satisfies this.
    */
    inline int findLimInfMeshPoint(double w, double* w_mesh, int nw)
    {
        return int((w - w_mesh[0])/(w_mesh[1] - w_mesh[0]));
    }
    // inline int findLimInfMeshPoint(double w, double* w_mesh, int nw)
    // {
    //     int lb = 0;
    //     int ub = nw-1;
    //     int mid;

    //     while (ub - lb > 1)
    //     {
    //         mid = (ub + lb)/2;
    //         // cout << "mid " << mid << endl;
    //         // cout << "wmesh mid " << w_mesh[mid] << endl;

    //         if (w_mesh[mid] > w)
    //         {
    //             ub = mid;
    //         }
    //         else
    //         {
    //             lb = mid;
    //         }
    //     }

    //     return lb;
    // }

    inline void normalize2D(double x[2]) {
        double norm = simutils::eucNorm2D(x);

        if (norm == 0.0) {
            norm += 1e-16;
        }

        x[0] = x[0] / norm;
        x[1] = x[1] / norm;
    }

    inline void normalize3D(double x[3]) {
        double norm = simutils::eucNorm3D(x);
        x[0] = x[0] / norm;
        x[1] = x[1] / norm;
        x[2] = x[2] / norm;
    }
}


#endif