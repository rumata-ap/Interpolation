using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Mathematics
{
    public class Spline1dCalcs
    {
        /*************************************************************************
        This subroutine calculates the value of the spline at the given point X.

        INPUT PARAMETERS:
            C   -   spline interpolant
            X   -   point

        Result:
            S(x)

          -- ALGLIB PROJECT --
             Copyright 23.06.2007 by Bochkanov Sergey
        *************************************************************************/
        public static double spline1dcalc(Spline1dinterpolant c, double x)
        {

            double result = Spline1d.spline1dcalc(c.innerobj, x);
            return result;
        }

        /*************************************************************************
        This subroutine differentiates the spline.

        INPUT PARAMETERS:
            C   -   spline interpolant.
            X   -   point

        Result:
            S   -   S(x)
            DS  -   S'(x)
            D2S -   S''(x)

          -- ALGLIB PROJECT --
             Copyright 24.06.2007 by Bochkanov Sergey
        *************************************************************************/
        public static void spline1ddiff(Spline1dinterpolant c, double x, out double s, out double ds, out double d2s)
        {
            s = 0;
            ds = 0;
            d2s = 0;
            Spline1d.spline1ddiff(c.innerobj, x, ref s, ref ds, ref d2s);
            return;
        }

        /*************************************************************************
        This subroutine integrates the spline.

        INPUT PARAMETERS:
            C   -   spline interpolant.
            X   -   right bound of the integration interval [a, x],
                    here 'a' denotes min(x[])
        Result:
            integral(S(t)dt,a,x)

          -- ALGLIB PROJECT --
             Copyright 23.06.2007 by Bochkanov Sergey
        *************************************************************************/
        public static double spline1dintegrate(Spline1dinterpolant c, double x)
        {

            double result = Spline1d.spline1dintegrate(c.innerobj, x);
            return result;
        }
    }
}
