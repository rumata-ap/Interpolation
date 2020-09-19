using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static Mathematics.alglib;

namespace Mathematics
{
    class Cspline : Spline1dCalcs
    {
        /*************************************************************************
        This subroutine builds cubic spline interpolant.

        INPUT PARAMETERS:
            X           -   spline nodes, array[0..N-1].
            Y           -   function values, array[0..N-1].

        OPTIONAL PARAMETERS:
            N           -   points count:
                            * N>=2
                            * if given, only first N points are used to build spline
                            * if not given, automatically detected from X/Y sizes
                              (len(X) must be equal to len(Y))
            BoundLType  -   boundary condition type for the left boundary
            BoundL      -   left boundary condition (first or second derivative,
                            depending on the BoundLType)
            BoundRType  -   boundary condition type for the right boundary
            BoundR      -   right boundary condition (first or second derivative,
                            depending on the BoundRType)

        OUTPUT PARAMETERS:
            C           -   spline interpolant

        ORDER OF POINTS

        Subroutine automatically sorts points, so caller may pass unsorted array.

        SETTING BOUNDARY VALUES:

        The BoundLType/BoundRType parameters can have the following values:
            * -1, which corresonds to the periodic (cyclic) boundary conditions.
                  In this case:
                  * both BoundLType and BoundRType must be equal to -1.
                  * BoundL/BoundR are ignored
                  * Y[last] is ignored (it is assumed to be equal to Y[first]).
            *  0, which  corresponds  to  the  parabolically   terminated  spline
                  (BoundL and/or BoundR are ignored).
            *  1, which corresponds to the first derivative boundary condition
            *  2, which corresponds to the second derivative boundary condition
            *  by default, BoundType=0 is used

        PROBLEMS WITH PERIODIC BOUNDARY CONDITIONS:

        Problems with periodic boundary conditions have Y[first_point]=Y[last_point].
        However, this subroutine doesn't require you to specify equal  values  for
        the first and last points - it automatically forces them  to  be  equal by
        copying  Y[first_point]  (corresponds  to the leftmost,  minimal  X[])  to
        Y[last_point]. However it is recommended to pass consistent values of Y[],
        i.e. to make Y[first_point]=Y[last_point].

          -- ALGLIB PROJECT --
             Copyright 23.06.2007 by Bochkanov Sergey
        *************************************************************************/
        public static void spline1dbuildcubic(double[] x, double[] y, int n, int boundltype, double boundl, int boundrtype, double boundr, out Spline1dinterpolant c)
        {
            c = new Spline1dinterpolant();
            Spline1d.spline1dbuildcubic(x, y, n, boundltype, boundl, boundrtype, boundr, c.innerobj);
            return;
        }
        public static void spline1dbuildcubic(double[] x, double[] y, out Spline1dinterpolant c)
        {
            int n;
            int boundltype;
            double boundl;
            int boundrtype;
            double boundr;
            if ((ap.len(x) != ap.len(y)))
                throw new alglibexception("Error while calling 'spline1dbuildcubic': looks like one of arguments has wrong size");
            c = new Spline1dinterpolant();
            n = ap.len(x);
            boundltype = 0;
            boundl = 0;
            boundrtype = 0;
            boundr = 0;
            Spline1d.spline1dbuildcubic(x, y, n, boundltype, boundl, boundrtype, boundr, c.innerobj);

            return;
        }
    }
}
