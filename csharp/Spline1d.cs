using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static Mathematics.alglib;

namespace Mathematics
{
    public class Spline1d
    {
        /*************************************************************************
        1-мерный интерполянт сплайна
        *************************************************************************/
        public class spline1dinterpolant : apobject
        {
            public bool periodic;
            public int n;
            public int k;
            public int continuity;
            public double[] x;
            public double[] c;
            public spline1dinterpolant()
            {
                init();
            }
            public override void init()
            {
                x = new double[0];
                c = new double[0];
            }
            public override alglib.apobject make_copy()
            {
                spline1dinterpolant _result = new spline1dinterpolant();
                _result.periodic = periodic;
                _result.n = n;
                _result.k = k;
                _result.continuity = continuity;
                _result.x = (double[])x.Clone();
                _result.c = (double[])c.Clone();
                return _result;
            }
        };

        /*************************************************************************
        Эта подпрограмма строит интеполянт линейного сплайна

        ВХОДНЫЕ ПАРАМЕТРЫ:
            X   -   узлы сплайна, массив[0..N-1]
            Y   -   значения функци, массив[0..N-1]
            N   -   число точек (опционально):
                    * N>=2
                    * если задано, для построения сплайна используются только первые N точек
                    * если не указано, автоматически определяется по размерам X / Y
                      (length(X) должно быть равно length(Y))

        ВЫХОДНЫЕ ПАРАМЕТРЫ:
            C   -   интерполянт сплайна


        ПОРЯДОК ТОЧЕК:

        Подпрограмма автоматически сортирует точки, поэтому вызывающий может передать несортированный массив.

          -- ALGLIB PROJECT --
             Copyright 24.06.2007 by Bochkanov Sergey
        *************************************************************************/
        public static void spline1dbuildlinear(double[] x, double[] y, int n, spline1dinterpolant c)
        {
            int i = 0;

            x = (double[])x.Clone();
            y = (double[])y.Clone();

            alglib.ap.assert(n > 1, "Spline1DBuildLinear: N<2!");
            alglib.ap.assert(alglib.ap.len(x) >= n, "Spline1DBuildLinear: Length(X)<N!");
            alglib.ap.assert(alglib.ap.len(y) >= n, "Spline1DBuildLinear: Length(Y)<N!");

            //
            // check and sort points
            //
            alglib.ap.assert(apserv.isfinitevector(x, n), "Spline1DBuildLinear: X contains infinite or NAN values!");
            alglib.ap.assert(apserv.isfinitevector(y, n), "Spline1DBuildLinear: Y contains infinite or NAN values!");
            Heapsortpoints(ref x, ref y, n);
            alglib.ap.assert(apserv.aredistinct(x, n), "Spline1DBuildLinear: at least two consequent points are too close!");

            //
            // Build
            //
            c.periodic = false;
            c.n = n;
            c.k = 3;
            c.continuity = 0;
            c.x = new double[n];
            c.c = new double[4 * (n - 1) + 2];
            for (i = 0; i <= n - 1; i++)
            {
                c.x[i] = x[i];
            }
            for (i = 0; i <= n - 2; i++)
            {
                c.c[4 * i + 0] = y[i];
                c.c[4 * i + 1] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
                c.c[4 * i + 2] = 0;
                c.c[4 * i + 3] = 0;
            }
            c.c[4 * (n - 1) + 0] = y[n - 1];
            c.c[4 * (n - 1) + 1] = c.c[4 * (n - 2) + 1];
        }

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
        public static void spline1dbuildcubic(double[] x,
            double[] y,
            int n,
            int boundltype,
            double boundl,
            int boundrtype,
            double boundr,
            spline1dinterpolant c)
        {
            double[] a1 = new double[0];
            double[] a2 = new double[0];
            double[] a3 = new double[0];
            double[] b = new double[0];
            double[] dt = new double[0];
            double[] d = new double[0];
            int[] p = new int[0];
            int ylen = 0;

            x = (double[])x.Clone();
            y = (double[])y.Clone();


            //
            // check correctness of boundary conditions
            //
            alglib.ap.assert(((boundltype == -1 || boundltype == 0) || boundltype == 1) || boundltype == 2, "Spline1DBuildCubic: incorrect BoundLType!");
            alglib.ap.assert(((boundrtype == -1 || boundrtype == 0) || boundrtype == 1) || boundrtype == 2, "Spline1DBuildCubic: incorrect BoundRType!");
            alglib.ap.assert((boundrtype == -1 && boundltype == -1) || (boundrtype != -1 && boundltype != -1), "Spline1DBuildCubic: incorrect BoundLType/BoundRType!");
            if (boundltype == 1 || boundltype == 2)
            {
                alglib.ap.assert(math.isfinite(boundl), "Spline1DBuildCubic: BoundL is infinite or NAN!");
            }
            if (boundrtype == 1 || boundrtype == 2)
            {
                alglib.ap.assert(math.isfinite(boundr), "Spline1DBuildCubic: BoundR is infinite or NAN!");
            }

            //
            // check lengths of arguments
            //
            alglib.ap.assert(n >= 2, "Spline1DBuildCubic: N<2!");
            alglib.ap.assert(alglib.ap.len(x) >= n, "Spline1DBuildCubic: Length(X)<N!");
            alglib.ap.assert(alglib.ap.len(y) >= n, "Spline1DBuildCubic: Length(Y)<N!");

            //
            // check and sort points
            //
            ylen = n;
            if (boundltype == -1)
            {
                ylen = n - 1;
            }
            alglib.ap.assert(apserv.isfinitevector(x, n), "Spline1DBuildCubic: X contains infinite or NAN values!");
            alglib.ap.assert(apserv.isfinitevector(y, ylen), "Spline1DBuildCubic: Y contains infinite or NAN values!");
            heapsortppoints(ref x, ref y, ref p, n);
            alglib.ap.assert(apserv.aredistinct(x, n), "Spline1DBuildCubic: at least two consequent points are too close!");

            //
            // Now we've checked and preordered everything,
            // so we can call internal function to calculate derivatives,
            // and then build Hermite spline using these derivatives
            //
            if (boundltype == -1 || boundrtype == -1)
            {
                y[n - 1] = y[0];
            }
            spline1dgriddiffcubicinternal(x, ref y, n, boundltype, boundl, boundrtype, boundr, ref d, ref a1, ref a2, ref a3, ref b, ref dt);
            spline1dbuildhermite(x, y, d, n, c);
            c.periodic = boundltype == -1 || boundrtype == -1;
            c.continuity = 2;
        }

        /*************************************************************************
        This subroutine builds Hermite spline interpolant.

        INPUT PARAMETERS:
            X           -   spline nodes, array[0..N-1]
            Y           -   function values, array[0..N-1]
            D           -   derivatives, array[0..N-1]
            N           -   points count (optional):
                            * N>=2
                            * if given, only first N points are used to build spline
                            * if not given, automatically detected from X/Y sizes
                              (len(X) must be equal to len(Y))

        OUTPUT PARAMETERS:
            C           -   spline interpolant.


        ORDER OF POINTS

        Subroutine automatically sorts points, so caller may pass unsorted array.

          -- ALGLIB PROJECT --
             Copyright 23.06.2007 by Bochkanov Sergey
        *************************************************************************/
        public static void spline1dbuildhermite(double[] x,
            double[] y,
            double[] d,
            int n,
            spline1dinterpolant c)
        {
            int i = 0;
            double delta = 0;
            double delta2 = 0;
            double delta3 = 0;

            x = (double[])x.Clone();
            y = (double[])y.Clone();
            d = (double[])d.Clone();

            alglib.ap.assert(n >= 2, "Spline1DBuildHermite: N<2!");
            alglib.ap.assert(alglib.ap.len(x) >= n, "Spline1DBuildHermite: Length(X)<N!");
            alglib.ap.assert(alglib.ap.len(y) >= n, "Spline1DBuildHermite: Length(Y)<N!");
            alglib.ap.assert(alglib.ap.len(d) >= n, "Spline1DBuildHermite: Length(D)<N!");

            //
            // check and sort points
            //
            alglib.ap.assert(apserv.isfinitevector(x, n), "Spline1DBuildHermite: X contains infinite or NAN values!");
            alglib.ap.assert(apserv.isfinitevector(y, n), "Spline1DBuildHermite: Y contains infinite or NAN values!");
            alglib.ap.assert(apserv.isfinitevector(d, n), "Spline1DBuildHermite: D contains infinite or NAN values!");
            heapsortdpoints(ref x, ref y, ref d, n);
            alglib.ap.assert(apserv.aredistinct(x, n), "Spline1DBuildHermite: at least two consequent points are too close!");

            //
            // Build
            //
            c.x = new double[n];
            c.c = new double[4 * (n - 1) + 2];
            c.periodic = false;
            c.k = 3;
            c.n = n;
            c.continuity = 1;
            for (i = 0; i <= n - 1; i++)
            {
                c.x[i] = x[i];
            }
            for (i = 0; i <= n - 2; i++)
            {
                delta = x[i + 1] - x[i];
                delta2 = math.sqr(delta);
                delta3 = delta * delta2;
                c.c[4 * i + 0] = y[i];
                c.c[4 * i + 1] = d[i];
                c.c[4 * i + 2] = (3 * (y[i + 1] - y[i]) - 2 * d[i] * delta - d[i + 1] * delta) / delta2;
                c.c[4 * i + 3] = (2 * (y[i] - y[i + 1]) + d[i] * delta + d[i + 1] * delta) / delta3;
            }
            c.c[4 * (n - 1) + 0] = y[n - 1];
            c.c[4 * (n - 1) + 1] = d[n - 1];
        }

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
        public static void spline1dbuildcatmullrom(double[] x,
            double[] y,
            int n,
            int boundtype,
            double tension,
            spline1dinterpolant c)
        {
            double[] d = new double[0];
            int i = 0;

            x = (double[])x.Clone();
            y = (double[])y.Clone();

            alglib.ap.assert(n >= 2, "Spline1DBuildCatmullRom: N<2!");
            alglib.ap.assert(boundtype == -1 || boundtype == 0, "Spline1DBuildCatmullRom: incorrect BoundType!");
            alglib.ap.assert((double)(tension) >= (double)(0), "Spline1DBuildCatmullRom: Tension<0!");
            alglib.ap.assert((double)(tension) <= (double)(1), "Spline1DBuildCatmullRom: Tension>1!");
            alglib.ap.assert(alglib.ap.len(x) >= n, "Spline1DBuildCatmullRom: Length(X)<N!");
            alglib.ap.assert(alglib.ap.len(y) >= n, "Spline1DBuildCatmullRom: Length(Y)<N!");

            //
            // check and sort points
            //
            alglib.ap.assert(apserv.isfinitevector(x, n), "Spline1DBuildCatmullRom: X contains infinite or NAN values!");
            alglib.ap.assert(apserv.isfinitevector(y, n), "Spline1DBuildCatmullRom: Y contains infinite or NAN values!");
            heapsortpoints(ref x, ref y, n);
            alglib.ap.assert(apserv.aredistinct(x, n), "Spline1DBuildCatmullRom: at least two consequent points are too close!");

            //
            // Special cases:
            // * N=2, parabolic terminated boundary condition on both ends
            // * N=2, periodic boundary condition
            //
            if (n == 2 && boundtype == 0)
            {

                //
                // Just linear spline
                //
                spline1dbuildlinear(x, y, n, c);
                return;
            }
            if (n == 2 && boundtype == -1)
            {

                //
                // Same as cubic spline with periodic conditions
                //
                spline1dbuildcubic(x, y, n, -1, 0.0, -1, 0.0, c);
                return;
            }

            //
            // Periodic or non-periodic boundary conditions
            //
            if (boundtype == -1)
            {

                //
                // Periodic boundary conditions
                //
                y[n - 1] = y[0];
                d = new double[n];
                d[0] = (y[1] - y[n - 2]) / (2 * (x[1] - x[0] + x[n - 1] - x[n - 2]));
                for (i = 1; i <= n - 2; i++)
                {
                    d[i] = (1 - tension) * (y[i + 1] - y[i - 1]) / (x[i + 1] - x[i - 1]);
                }
                d[n - 1] = d[0];

                //
                // Now problem is reduced to the cubic Hermite spline
                //
                spline1dbuildhermite(x, y, d, n, c);
                c.periodic = true;
            }
            else
            {

                //
                // Non-periodic boundary conditions
                //
                d = new double[n];
                for (i = 1; i <= n - 2; i++)
                {
                    d[i] = (1 - tension) * (y[i + 1] - y[i - 1]) / (x[i + 1] - x[i - 1]);
                }
                d[0] = 2 * (y[1] - y[0]) / (x[1] - x[0]) - d[1];
                d[n - 1] = 2 * (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]) - d[n - 2];

                //
                // Now problem is reduced to the cubic Hermite spline
                //
                spline1dbuildhermite(x, y, d, n, c);
            }
        }


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
        public static void spline1dbuildakima(double[] x, double[] y, int n, spline1dinterpolant c)
        {
            int i = 0;
            double[] d = new double[0];
            double[] w = new double[0];
            double[] diff = new double[0];

            x = (double[])x.Clone();
            y = (double[])y.Clone();

            alglib.ap.assert(n >= 2, "Spline1DBuildAkima: N<2!");
            alglib.ap.assert(alglib.ap.len(x) >= n, "Spline1DBuildAkima: Length(X)<N!");
            alglib.ap.assert(alglib.ap.len(y) >= n, "Spline1DBuildAkima: Length(Y)<N!");

            //
            // check and sort points
            //
            alglib.ap.assert(apserv.isfinitevector(x, n), "Spline1DBuildAkima: X contains infinite or NAN values!");
            alglib.ap.assert(apserv.isfinitevector(y, n), "Spline1DBuildAkima: Y contains infinite or NAN values!");
            heapsortpoints(ref x, ref y, n);
            alglib.ap.assert(apserv.aredistinct(x, n), "Spline1DBuildAkima: at least two consequent points are too close!");

            //
            // Handle special cases: N=2, N=3, N=4
            //
            if (n <= 4)
            {
                spline1dbuildcubic(x, y, n, 0, 0.0, 0, 0.0, c);
                return;
            }

            //
            // Prepare W (weights), Diff (divided differences)
            //
            w = new double[n - 1];
            diff = new double[n - 1];
            for (i = 0; i <= n - 2; i++)
            {
                diff[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
            }
            for (i = 1; i <= n - 2; i++)
            {
                w[i] = Math.Abs(diff[i] - diff[i - 1]);
            }

            //
            // Prepare Hermite interpolation scheme
            //
            d = new double[n];
            for (i = 2; i <= n - 3; i++)
            {
                if ((double)(Math.Abs(w[i - 1]) + Math.Abs(w[i + 1])) != (double)(0))
                {
                    d[i] = (w[i + 1] * diff[i - 1] + w[i - 1] * diff[i]) / (w[i + 1] + w[i - 1]);
                }
                else
                {
                    d[i] = ((x[i + 1] - x[i]) * diff[i - 1] + (x[i] - x[i - 1]) * diff[i]) / (x[i + 1] - x[i - 1]);
                }
            }
            d[0] = diffthreepoint(x[0], x[0], y[0], x[1], y[1], x[2], y[2]);
            d[1] = diffthreepoint(x[1], x[0], y[0], x[1], y[1], x[2], y[2]);
            d[n - 2] = diffthreepoint(x[n - 2], x[n - 3], y[n - 3], x[n - 2], y[n - 2], x[n - 1], y[n - 1]);
            d[n - 1] = diffthreepoint(x[n - 1], x[n - 3], y[n - 3], x[n - 2], y[n - 2], x[n - 1], y[n - 1]);

            //
            // Build Akima spline using Hermite interpolation scheme
            //
            spline1dbuildhermite(x, y, d, n, c);
        }


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
        public static double spline1dcalc(spline1dinterpolant c, double x)
        {
            double result = 0;
            int l = 0;
            int r = 0;
            int m = 0;
            double t = 0;

            alglib.ap.assert(c.k == 3, "Spline1DCalc: internal error");
            alglib.ap.assert(!Double.IsInfinity(x), "Spline1DCalc: infinite X!");

            //
            // special case: NaN
            //
            if (Double.IsNaN(x))
            {
                result = Double.NaN;
                return result;
            }

            //
            // correct if periodic
            //
            if (c.periodic)
            {
                apserv.apperiodicmap(ref x, c.x[0], c.x[c.n - 1], ref t);
            }

            //
            // Binary search in the [ x[0], ..., x[n-2] ] (x[n-1] is not included)
            //
            l = 0;
            r = c.n - 2 + 1;
            while (l != r - 1)
            {
                m = (l + r) / 2;
                if (c.x[m] >= x)
                {
                    r = m;
                }
                else
                {
                    l = m;
                }
            }

            //
            // Interpolation
            //
            x = x - c.x[l];
            m = 4 * l;
            result = c.c[m] + x * (c.c[m + 1] + x * (c.c[m + 2] + x * c.c[m + 3]));
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
        public static void spline1ddiff(spline1dinterpolant c,
            double x,
            ref double s,
            ref double ds,
            ref double d2s)
        {
            int l = 0;
            int r = 0;
            int m = 0;
            double t = 0;

            s = 0;
            ds = 0;
            d2s = 0;

            alglib.ap.assert(c.k == 3, "Spline1DDiff: internal error");
            alglib.ap.assert(!Double.IsInfinity(x), "Spline1DDiff: infinite X!");

            //
            // special case: NaN
            //
            if (Double.IsNaN(x))
            {
                s = Double.NaN;
                ds = Double.NaN;
                d2s = Double.NaN;
                return;
            }

            //
            // correct if periodic
            //
            if (c.periodic)
            {
                apserv.apperiodicmap(ref x, c.x[0], c.x[c.n - 1], ref t);
            }

            //
            // Binary search
            //
            l = 0;
            r = c.n - 2 + 1;
            while (l != r - 1)
            {
                m = (l + r) / 2;
                if (c.x[m] >= x)
                {
                    r = m;
                }
                else
                {
                    l = m;
                }
            }

            //
            // Differentiation
            //
            x = x - c.x[l];
            m = 4 * l;
            s = c.c[m] + x * (c.c[m + 1] + x * (c.c[m + 2] + x * c.c[m + 3]));
            ds = c.c[m + 1] + 2 * x * c.c[m + 2] + 3 * math.sqr(x) * c.c[m + 3];
            d2s = 2 * c.c[m + 2] + 6 * x * c.c[m + 3];
        }

        /*************************************************************************
        This subroutine integrates the spline.

        INPUT PARAMETERS:
            C   -   spline interpolant.
            X   -   правая граница интервала интегрирования [a, x],
                    здесь 'a' обозначает min(x[])
        Result:
            integral(S(t)dt,a,x)

          -- ALGLIB PROJECT --
             Copyright 23.06.2007 by Bochkanov Sergey
        *************************************************************************/
        public static double spline1dintegrate(spline1dinterpolant c, double x)
        {
            double result = 0;
            int n = 0;
            int i = 0;
            int j = 0;
            int l = 0;
            int r = 0;
            int m = 0;
            double w = 0;
            double v = 0;
            double t = 0;
            double intab = 0;
            double additionalterm = 0;

            n = c.n;

            //
            // Periodic splines require special treatment. We make
            // following transformation:
            //
            //     integral(S(t)dt,A,X) = integral(S(t)dt,A,Z)+AdditionalTerm
            //
            // here X may lie outside of [A,B], Z lies strictly in [A,B],
            // AdditionalTerm is equals to integral(S(t)dt,A,B) times some
            // integer number (may be zero).
            //
            if (c.periodic && ((double)(x) < (double)(c.x[0]) || (double)(x) > (double)(c.x[c.n - 1])))
            {

                //
                // compute integral(S(x)dx,A,B)
                //
                intab = 0;
                for (i = 0; i <= c.n - 2; i++)
                {
                    w = c.x[i + 1] - c.x[i];
                    m = (c.k + 1) * i;
                    intab = intab + c.c[m] * w;
                    v = w;
                    for (j = 1; j <= c.k; j++)
                    {
                        v = v * w;
                        intab = intab + c.c[m + j] * v / (j + 1);
                    }
                }

                //
                // map X into [A,B]
                //
                apserv.apperiodicmap(ref x, c.x[0], c.x[c.n - 1], ref t);
                additionalterm = t * intab;
            }
            else
            {
                additionalterm = 0;
            }

            //
            // Binary search in the [ x[0], ..., x[n-2] ] (x[n-1] is not included)
            //
            l = 0;
            r = n - 2 + 1;
            while (l != r - 1)
            {
                m = (l + r) / 2;
                if ((double)(c.x[m]) >= (double)(x))
                {
                    r = m;
                }
                else
                {
                    l = m;
                }
            }

            //
            // Integration
            //
            result = 0;
            for (i = 0; i <= l - 1; i++)
            {
                w = c.x[i + 1] - c.x[i];
                m = (c.k + 1) * i;
                result = result + c.c[m] * w;
                v = w;
                for (j = 1; j <= c.k; j++)
                {
                    v = v * w;
                    result = result + c.c[m + j] * v / (j + 1);
                }
            }
            w = x - c.x[l];
            m = (c.k + 1) * l;
            v = w;
            result = result + c.c[m] * w;
            for (j = 1; j <= c.k; j++)
            {
                v = v * w;
                result = result + c.c[m + j] * v / (j + 1);
            }
            result = result + additionalterm;
            return result;
        }

        /*************************************************************************
        Internal subroutine. Heap sort.
        *************************************************************************/
        private static void Heapsortpoints(ref double[] x,
            ref double[] y,
            int n)
        {
            double[] bufx = new double[0];
            double[] bufy = new double[0];

            tsort.tagsortfastr(ref x, ref y, ref bufx, ref bufy, n);
        }

        /*************************************************************************
        Internal subroutine. Heap sort.
        *************************************************************************/
        private static void heapsortpoints(ref double[] x,
            ref double[] y,
            int n)
        {
            double[] bufx = new double[0];
            double[] bufy = new double[0];

            tsort.tagsortfastr(ref x, ref y, ref bufx, ref bufy, n);
        }

        /*************************************************************************
          Внутренняя подпрограмма. Сортировка кучи.
          *************************************************************************/
        public static void heapsortdpoints(ref double[] x,
            ref double[] y,
            ref double[] d,
            int n)
        {
            double[] rbuf = new double[0];
            int[] ibuf = new int[0];
            double[] rbuf2 = new double[0];
            int[] ibuf2 = new int[0];
            int i = 0;
            int i_ = 0;

            ibuf = new int[n];
            rbuf = new double[n];
            for (i = 0; i <= n - 1; i++)
            {
                ibuf[i] = i;
            }
            tsort.tagsortfasti(ref x, ref ibuf, ref rbuf2, ref ibuf2, n);
            for (i = 0; i <= n - 1; i++)
            {
                rbuf[i] = y[ibuf[i]];
            }
            for (i_ = 0; i_ <= n - 1; i_++)
            {
                y[i_] = rbuf[i_];
            }
            for (i = 0; i <= n - 1; i++)
            {
                rbuf[i] = d[ibuf[i]];
            }
            for (i_ = 0; i_ <= n - 1; i_++)
            {
                d[i_] = rbuf[i_];
            }
        }

        /*************************************************************************
        Внутренняя подпрограмма. Сортировка кучи.

        Принимает:
            X, Y    -   точки
            P       -   пустой или предварительно выделенный массив

        Возвращает:
            X, Y    -   сортированные X
            P       -   массив перестановок; I-е позиции выходных массивов 
                        X/Y содержит (X[P[I]],Y[P[I]])
        *************************************************************************/
        private static void heapsortppoints(ref double[] x,
            ref double[] y,
            ref int[] p,
            int n)
        {
            double[] rbuf = new double[0];
            int[] ibuf = new int[0];
            int i = 0;
            int i_ = 0;

            if (alglib.ap.len(p) < n)
            {
                p = new int[n];
            }
            rbuf = new double[n];
            for (i = 0; i <= n - 1; i++)
            {
                p[i] = i;
            }
            tsort.tagsortfasti(ref x, ref p, ref rbuf, ref ibuf, n);
            for (i = 0; i <= n - 1; i++)
            {
                rbuf[i] = y[p[i]];
            }
            for (i_ = 0; i_ <= n - 1; i_++)
            {
                y[i_] = rbuf[i_];
            }
        }

        /*************************************************************************
        Internal version of Spline1DGridDiffCubic.

        Accepts pre-ordered X/Y, temporary arrays (which may be  preallocated,  if
        you want to save time, or not) and output array (which may be preallocated
        too).

        Y is passed as var-parameter because we may need to force last element  to
        be equal to the first one (if periodic boundary conditions are specified).

          -- ALGLIB PROJECT --
             Copyright 03.09.2010 by Bochkanov Sergey
        *************************************************************************/
        private static void spline1dgriddiffcubicinternal(double[] x,
            ref double[] y,
            int n,
            int boundltype,
            double boundl,
            int boundrtype,
            double boundr,
            ref double[] d,
            ref double[] a1,
            ref double[] a2,
            ref double[] a3,
            ref double[] b,
            ref double[] dt)
        {
            int i = 0;
            int i_ = 0;


            //
            // allocate arrays
            //
            if (alglib.ap.len(d) < n)
            {
                d = new double[n];
            }
            if (alglib.ap.len(a1) < n)
            {
                a1 = new double[n];
            }
            if (alglib.ap.len(a2) < n)
            {
                a2 = new double[n];
            }
            if (alglib.ap.len(a3) < n)
            {
                a3 = new double[n];
            }
            if (alglib.ap.len(b) < n)
            {
                b = new double[n];
            }
            if (alglib.ap.len(dt) < n)
            {
                dt = new double[n];
            }

            //
            // Special cases:
            // * N=2, parabolic terminated boundary condition on both ends
            // * N=2, periodic boundary condition
            //
            if ((n == 2 && boundltype == 0) && boundrtype == 0)
            {
                d[0] = (y[1] - y[0]) / (x[1] - x[0]);
                d[1] = d[0];
                return;
            }
            if ((n == 2 && boundltype == -1) && boundrtype == -1)
            {
                d[0] = 0;
                d[1] = 0;
                return;
            }

            //
            // Periodic and non-periodic boundary conditions are
            // two separate classes
            //
            if (boundrtype == -1 && boundltype == -1)
            {

                //
                // Periodic boundary conditions
                //
                y[n - 1] = y[0];

                //
                // Boundary conditions at N-1 points
                // (one point less because last point is the same as first point).
                //
                a1[0] = x[1] - x[0];
                a2[0] = 2 * (x[1] - x[0] + x[n - 1] - x[n - 2]);
                a3[0] = x[n - 1] - x[n - 2];
                b[0] = 3 * (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]) * (x[1] - x[0]) + 3 * (y[1] - y[0]) / (x[1] - x[0]) * (x[n - 1] - x[n - 2]);
                for (i = 1; i <= n - 2; i++)
                {

                    //
                    // Altough last point is [N-2], we use X[N-1] and Y[N-1]
                    // (because of periodicity)
                    //
                    a1[i] = x[i + 1] - x[i];
                    a2[i] = 2 * (x[i + 1] - x[i - 1]);
                    a3[i] = x[i] - x[i - 1];
                    b[i] = 3 * (y[i] - y[i - 1]) / (x[i] - x[i - 1]) * (x[i + 1] - x[i]) + 3 * (y[i + 1] - y[i]) / (x[i + 1] - x[i]) * (x[i] - x[i - 1]);
                }

                //
                // Solve, add last point (with index N-1)
                //
                solvecyclictridiagonal(a1, a2, a3, b, n - 1, ref dt);
                for (i_ = 0; i_ <= n - 2; i_++)
                {
                    d[i_] = dt[i_];
                }
                d[n - 1] = d[0];
            }
            else
            {

                //
                // Non-periodic boundary condition.
                // Left boundary conditions.
                //
                if (boundltype == 0)
                {
                    a1[0] = 0;
                    a2[0] = 1;
                    a3[0] = 1;
                    b[0] = 2 * (y[1] - y[0]) / (x[1] - x[0]);
                }
                if (boundltype == 1)
                {
                    a1[0] = 0;
                    a2[0] = 1;
                    a3[0] = 0;
                    b[0] = boundl;
                }
                if (boundltype == 2)
                {
                    a1[0] = 0;
                    a2[0] = 2;
                    a3[0] = 1;
                    b[0] = 3 * (y[1] - y[0]) / (x[1] - x[0]) - 0.5 * boundl * (x[1] - x[0]);
                }

                //
                // Central conditions
                //
                for (i = 1; i <= n - 2; i++)
                {
                    a1[i] = x[i + 1] - x[i];
                    a2[i] = 2 * (x[i + 1] - x[i - 1]);
                    a3[i] = x[i] - x[i - 1];
                    b[i] = 3 * (y[i] - y[i - 1]) / (x[i] - x[i - 1]) * (x[i + 1] - x[i]) + 3 * (y[i + 1] - y[i]) / (x[i + 1] - x[i]) * (x[i] - x[i - 1]);
                }

                //
                // Right boundary conditions
                //
                if (boundrtype == 0)
                {
                    a1[n - 1] = 1;
                    a2[n - 1] = 1;
                    a3[n - 1] = 0;
                    b[n - 1] = 2 * (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]);
                }
                if (boundrtype == 1)
                {
                    a1[n - 1] = 0;
                    a2[n - 1] = 1;
                    a3[n - 1] = 0;
                    b[n - 1] = boundr;
                }
                if (boundrtype == 2)
                {
                    a1[n - 1] = 1;
                    a2[n - 1] = 2;
                    a3[n - 1] = 0;
                    b[n - 1] = 3 * (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]) + 0.5 * boundr * (x[n - 1] - x[n - 2]);
                }

                //
                // Solve
                //
                solvetridiagonal(a1, a2, a3, b, n, ref d);
            }


        }

        /*************************************************************************
        Internal subroutine. Tridiagonal solver. Solves

        ( B[0] C[0]                      
        ( A[1] B[1] C[1]                 )
        (      A[2] B[2] C[2]            )
        (            ..........          ) * X = D
        (            ..........          )
        (           A[N-2] B[N-2] C[N-2] )
        (                  A[N-1] B[N-1] )

        *************************************************************************/
        private static void solvetridiagonal(double[] a,
            double[] b,
            double[] c,
            double[] d,
            int n,
            ref double[] x)
        {
            int k = 0;
            double t = 0;

            b = (double[])b.Clone();
            d = (double[])d.Clone();

            if (alglib.ap.len(x) < n)
            {
                x = new double[n];
            }
            for (k = 1; k <= n - 1; k++)
            {
                t = a[k] / b[k - 1];
                b[k] = b[k] - t * c[k - 1];
                d[k] = d[k] - t * d[k - 1];
            }
            x[n - 1] = d[n - 1] / b[n - 1];
            for (k = n - 2; k >= 0; k--)
            {
                x[k] = (d[k] - c[k] * x[k + 1]) / b[k];
            }
        }

        /*************************************************************************
        Internal subroutine. Cyclic tridiagonal solver. Solves

        ( B[0] C[0]                 A[0] )
        ( A[1] B[1] C[1]                 )
        (      A[2] B[2] C[2]            )
        (            ..........          ) * X = D
        (            ..........          )
        (           A[N-2] B[N-2] C[N-2] )
        ( C[N-1]           A[N-1] B[N-1] )
        *************************************************************************/
        private static void solvecyclictridiagonal(double[] a,
            double[] b,
            double[] c,
            double[] d,
            int n,
            ref double[] x)
        {
            int k = 0;
            double alpha = 0;
            double beta = 0;
            double gamma = 0;
            double[] y = new double[0];
            double[] z = new double[0];
            double[] u = new double[0];

            b = (double[])b.Clone();

            if (alglib.ap.len(x) < n)
            {
                x = new double[n];
            }
            beta = a[0];
            alpha = c[n - 1];
            gamma = -b[0];
            b[0] = 2 * b[0];
            b[n - 1] = b[n - 1] - alpha * beta / gamma;
            u = new double[n];
            for (k = 0; k <= n - 1; k++)
            {
                u[k] = 0;
            }
            u[0] = gamma;
            u[n - 1] = alpha;
            solvetridiagonal(a, b, c, d, n, ref y);
            solvetridiagonal(a, b, c, u, n, ref z);
            for (k = 0; k <= n - 1; k++)
            {
                x[k] = y[k] - (y[0] + beta / gamma * y[n - 1]) / (1 + z[0] + beta / gamma * z[n - 1]) * z[k];
            }
        }

        /*************************************************************************
        Internal subroutine. Three-point differentiation
        *************************************************************************/
        private static double diffthreepoint(double t,
            double x0,
            double f0,
            double x1,
            double f1,
            double x2,
            double f2)
        {
            double result = 0;
            double a = 0;
            double b = 0;

            t = t - x0;
            x1 = x1 - x0;
            x2 = x2 - x0;
            a = (f2 - f0 - x2 / x1 * (f1 - f0)) / (math.sqr(x2) - x1 * x2);
            b = (f1 - f0 - a * math.sqr(x1)) / x1;
            result = 2 * a * t + b;
            return result;
        }
    }
}
