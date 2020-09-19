using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static Mathematics.alglib;

namespace Mathematics
{
    class CRspline : Spline1dCalcs
    {
        /*************************************************************************
        This subroutine builds Catmull-Rom spline interpolant.

        INPUT PARAMETERS:
            X           -   spline nodes, array[0..N-1].
            Y           -   function values, array[0..N-1].

        OPTIONAL PARAMETERS:
            N           -   points count:
                            * N>=2
                            * if given, only first N points are used to build spline
                            * if not given, automatically detected from X/Y sizes
                              (len(X) must be equal to len(Y))
            BoundType   -   boundary condition type:
                            * -1 for periodic boundary condition
                            *  0 for parabolically terminated spline (default)
            Tension     -   tension parameter:
                            * tension=0   corresponds to classic Catmull-Rom spline (default)
                            * 0<tension<1 corresponds to more general form - cardinal spline

        OUTPUT PARAMETERS:
            C           -   spline interpolant


        ORDER OF POINTS

        Subroutine automatically sorts points, so caller may pass unsorted array.

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
        public static void spline1dbuildcatmullrom(double[] x, double[] y, int n, int boundtype, double tension, out Spline1dinterpolant c)
        {
            c = new Spline1dinterpolant();
            Spline1d.spline1dbuildcatmullrom(x, y, n, boundtype, tension, c.innerobj);
            return;
        }
        public static void spline1dbuildcatmullrom(double[] x, double[] y, out Spline1dinterpolant c)
        {
            int n;
            int boundtype;
            double tension;
            if ((ap.len(x) != ap.len(y)))
                throw new alglibexception("Error while calling 'spline1dbuildcatmullrom': looks like one of arguments has wrong size");
            c = new Spline1dinterpolant();
            n = ap.len(x);
            boundtype = 0;
            tension = 0;
            Spline1d.spline1dbuildcatmullrom(x, y, n, boundtype, tension, c.innerobj);

            return;
        }
    }
}
