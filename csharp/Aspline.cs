using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static Mathematics.alglib;

namespace Mathematics
{
    class Aspline : Spline1dCalcs
    {
        /*************************************************************************
        This subroutine builds Akima spline interpolant

        INPUT PARAMETERS:
            X           -   spline nodes, array[0..N-1]
            Y           -   function values, array[0..N-1]
            N           -   points count (optional):
                            * N>=2
                            * if given, only first N points are used to build spline
                            * if not given, automatically detected from X/Y sizes
                              (len(X) must be equal to len(Y))

        OUTPUT PARAMETERS:
            C           -   spline interpolant


        ORDER OF POINTS

        Subroutine automatically sorts points, so caller may pass unsorted array.

          -- ALGLIB PROJECT --
             Copyright 24.06.2007 by Bochkanov Sergey
        *************************************************************************/
        public static void spline1dbuildakima(double[] x, double[] y, int n, out Spline1dinterpolant c)
        {
            c = new Spline1dinterpolant();
            Spline1d.spline1dbuildakima(x, y, n, c.innerobj);
            return;
        }
        public static void spline1dbuildakima(double[] x, double[] y, out Spline1dinterpolant c)
        {
            int n;
            if ((ap.len(x) != ap.len(y)))
                throw new alglibexception("Error while calling 'spline1dbuildakima': looks like one of arguments has wrong size");
            c = new Spline1dinterpolant();
            n = ap.len(x);
            Spline1d.spline1dbuildakima(x, y, n, c.innerobj);

            return;
        }

    }
}
